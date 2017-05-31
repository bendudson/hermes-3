
#include "mixed.hxx"

#include <options.hxx>
#include <boutexception.hxx>

#include "div_ops.hxx"

NeutralMixed::NeutralMixed(Solver *solver, Mesh *mesh, Options *options) : NeutralModel(options) {
  if(options == NULL) 
    options = Options::getRoot()->getSection("neutral");
  
  solver->add(Nn, "Nn");
  solver->add(Pn, "Pn");
  solver->add(NVn, "NVn");
  
  OPTION(options, sheath_ydown, true);
  OPTION(options, sheath_yup,   true);
  
  OPTION(options, neutral_gamma, 5./4);
  
  OPTION(options, numdiff, 1.0);
  
  // Optionally output time derivatives
  bool output_ddt;
  OPTION(options, output_ddt, false);
  if(output_ddt) {
    SAVE_REPEAT3(ddt(Nn), ddt(Pn), ddt(NVn));
  }
  
  Dnn = 0.0; // Neutral gas diffusion
  
  S = 0;
  F = 0;
  Qi = 0;
  Rp = 0;
  
  SAVE_REPEAT(Dnn);
}

void NeutralMixed::update(const Field3D &Ne, const Field3D &Te, const Field3D &Ti, const Field3D &Vi) {
  TRACE("NeutralMixed::update");
  
  mesh->communicate(Nn, Pn, NVn);
  
  //////////////////////////////////////////////////////
  // 3D model, diffusive in X-Z, fluid in Y
  // 
  // Evolves neutral density Nn, pressure Pn, and
  // parallel momentum NVn
  // 
  
  Nn = floor(Nn, 1e-8);
  Field3D Nnlim = floor(Nn, 1e-5); // Used where division by neutral density is needed
  Field3D Tn = Pn / Nn; 
  Tn = floor(Tn, 0.01/Tnorm);
  
  Field3D Vn = NVn / Nnlim; // Neutral parallel velocity
  
  Field3D Pnlim = Nn * Tn;
  Pnlim.applyBoundary("neumann");
  
  /////////////////////////////////////////////////////
  // Boundary conditions
  TRACE("Neutral boundary conditions");
  
  if(sheath_ydown) {
    for(RangeIterator r=mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      for(int jz=0; jz<mesh->ngz; jz++) {
        // Free boundary (constant gradient) density
        BoutReal nnwall = 0.5*( 3.*Nn(r.ind, mesh->ystart, jz) - Nn(r.ind, mesh->ystart+1, jz) );
        if(nnwall < 0.0)
          nnwall = 0.0;
        
        BoutReal tnwall = Tn(r.ind, mesh->ystart,jz);
        
        Nn(r.ind, mesh->ystart-1,jz) = 2*nnwall - Nn(r.ind, mesh->ystart, jz);
        
        // Zero gradient temperature, heat flux added later
        Tn(r.ind, mesh->ystart-1,jz) = tnwall;
        
        // Set pressure consistent at the boundary
        Pn(r.ind, mesh->ystart-1,jz) = 2.*nnwall*tnwall - Pn(r.ind, mesh->ystart,jz);
        // No flow into wall
        Vn(r.ind, mesh->ystart-1,jz) = -Vn(r.ind, mesh->ystart,jz);
        NVn(r.ind, mesh->ystart-1,jz) = -NVn(r.ind, mesh->ystart,jz);
      }
    }
  }
  
  if(sheath_yup) {
    for(RangeIterator r=mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      for(int jz=0; jz<mesh->ngz; jz++) {
        // Free boundary (constant gradient) density
        BoutReal nnwall = 0.5*( 3.*Nn(r.ind, mesh->yend, jz) - Nn(r.ind, mesh->yend-1, jz) );
        if(nnwall < 0.0)
          nnwall = 0.0;
        
        BoutReal tnwall = Tn(r.ind, mesh->yend,jz);
        
        Nn(r.ind, mesh->yend+1,jz) = 2*nnwall - Nn(r.ind, mesh->yend, jz);
        
        // Zero gradient temperature, heat flux added later
        Tn(r.ind, mesh->yend+1,jz) = tnwall;
        
        // Set pressure consistent at the boundary
        Pn(r.ind, mesh->yend+1,jz) = 2.*nnwall*tnwall - Pn(r.ind, mesh->yend,jz);
        // No flow into wall
        Vn(r.ind, mesh->yend+1,jz) = -Vn(r.ind, mesh->yend,jz);
        NVn(r.ind, mesh->yend+1,jz) = -NVn(r.ind, mesh->yend,jz);
      }
    }
  }
  
  /////////////////////////////////////////////////////
  // Atomic processes
  TRACE("Atomic processes");
  
  Field3D Riz, Rrc, Rcx;
  neutral_rates(Ne, Te, Ti, Vi, 
                Nn, Tn, Vn,
                S, F, Qi, Rp,
                Riz, Rrc, Rcx);
  
  // Neutral cross-field diffusion coefficient
  BoutReal neutral_lmax = 0.1 / Lnorm;
  Field3D Rnn = Nn*sqrt(Tn)/neutral_lmax; // Neutral-neutral collisions
  Dnn = Pnlim / (Riz + Rcx + Rnn);
  mesh->communicate(Dnn);
  Dnn.applyBoundary("dirichlet_o2");
  
  // Logarithms used to calculate perpendicular velocity
  // V_perp = -Dnn * ( Grad_perp(Nn)/Nn + Grad_perp(Tn)/Tn )
  //
  // Grad(Pn) / Pn = Grad(Tn)/Tn + Grad(Nn)/Nn
  //               = Grad(logTn + logNn)
  //Field3D logNn = log(Nn);
  //Field3D logTn = log(Tn);
  
  Field3D logPnlim = log(Pnlim);
  logPnlim.applyBoundary("neumann");
      
  /////////////////////////////////////////////////////
  // Neutral density
  TRACE("Neutral density");
  
  ddt(Nn) = 
    - Div_par_FV(Nn, Vn)        // Advection
    + S                         // Source from recombining plasma
    + Div_Perp_Lap_XYZ(Dnn * Nn, logPnlim, false) // Perpendicular diffusion
    ;
  
  /////////////////////////////////////////////////////
  // Neutral momentum
  TRACE("Neutral momentum");
  
  ddt(NVn) = 
    - Div_par_FV(NVn, Vn)        // Momentum flow
    + F                          // Friction with plasma
    - Grad_par(Pnlim)            // Pressure gradient
    + Div_Perp_Lap_XYZ(Dnn * NVn, logPnlim, false) // Perpendicular diffusion
    
    + Div_par_diffusion(Dnn*Nn, Vn, false) // Parallel viscosity
    ;
  
  if(numdiff > 0.0) {
    ddt(NVn) += Div_par_diffusion(SQ(mesh->dy)*mesh->g_22*numdiff, Vn);
  }
  
  Fperp = Rcx + Rrc;
  
  /////////////////////////////////////////////////////
  // Neutral pressure
  TRACE("Neutral pressure");
  
  ddt(Pn) =
    - Div_par_FV(Pn, Vn)                   // Advection
    - (2./3)*Pnlim*Div_par(Vn)             // Compression
    + (2./3) * Qi
    + Div_Perp_Lap_XYZ(Dnn * Pn, logPnlim, false) // Perpendicular diffusion
    + Div_Perp_Lap_XYZ(Dnn * Nn, Tn, false) // Conduction
    + Div_par_diffusion(Dnn * Nn, Tn)       // Parallel conduction
    ;
  
  //////////////////////////////////////////////////////
  // Boundary condition on fluxes
  TRACE("Neutral boundary fluxes");
    
  if(sheath_ydown) {
    for(RangeIterator r=mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      for(int jz=0;jz<mesh->ngz-1;jz++) {
        
        // Loss of thermal energy to the target.
        // This depends on the reflection coefficient
        // and is controlled by the option neutral_gamma
        //         q = neutral_gamma * n * T * cs
        
        // Density at the target
        BoutReal Nnout = 0.5*(Nn(r.ind, mesh->ystart,jz) + Nn(r.ind, mesh->ystart-1,jz));
        if(Nnout < 0.0)
          Nnout = 0.0;
        // Temperature at the target
        BoutReal Tnout = 0.5*(Tn(r.ind, mesh->ystart,jz) + Tn(r.ind, mesh->ystart-1,jz));
        if(Tnout < 0.0)
          Tnout = 0.0;
        
        // gamma * n * T * cs
        BoutReal q = neutral_gamma * Nnout * Tnout * sqrt(Tnout);
        // Multiply by cell area to get power
        BoutReal heatflux = q * (mesh->J(r.ind, mesh->ystart)+mesh->J(r.ind, mesh->ystart-1))/(sqrt(mesh->g_22(r.ind, mesh->ystart)) + sqrt(mesh->g_22(r.ind, mesh->ystart-11)));
        
        // Divide by volume of cell, and multiply by 2/3 to get pressure
        ddt(Pn)(r.ind, mesh->ystart,jz) -= (2./3)*heatflux / (mesh->dy(r.ind, mesh->ystart)*mesh->J(r.ind, mesh->ystart));
      }
    }
  }
      
  if(sheath_yup) {
    for(RangeIterator r=mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      for(int jz=0;jz<mesh->ngz-1;jz++) {
        
        // Loss of thermal energy to the target.
        // This depends on the reflection coefficient
        // and is controlled by the option neutral_gamma
        //         q = neutral_gamma * n * T * cs
        
        // Density at the target
        BoutReal Nnout = 0.5*(Nn(r.ind, mesh->yend,jz) + Nn(r.ind, mesh->yend+1,jz));
        if(Nnout < 0.0)
          Nnout = 0.0;
        // Temperature at the target
        BoutReal Tnout = 0.5*(Tn(r.ind, mesh->yend,jz) + Tn(r.ind, mesh->yend+1,jz));
        if(Tnout < 0.0)
          Tnout = 0.0;
        
        // gamma * n * T * cs
        BoutReal q = neutral_gamma * Nnout * Tnout * sqrt(Tnout);
        // Multiply by cell area to get power
        BoutReal heatflux = q * (mesh->J(r.ind, mesh->yend)+mesh->J(r.ind, mesh->yend+1))/(sqrt(mesh->g_22(r.ind, mesh->yend)) + sqrt(mesh->g_22(r.ind, mesh->yend+1)));
        
        // Divide by volume of cell, and multiply by 2/3 to get pressure
        ddt(Pn)(r.ind, mesh->yend,jz) -= (2./3)*heatflux / (mesh->dy(r.ind, mesh->yend)*mesh->J(r.ind, mesh->yend));
      }
    }
  }
  
}

