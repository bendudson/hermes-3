#include "recycling.hxx"

#include <options.hxx>
#include <boutexception.hxx>

NeutralRecycling::NeutralRecycling(Solver *solver, Mesh *mesh, Options *options) : NeutralModel(options) {
  if(options == NULL) 
    options = Options::getRoot()->getSection("neutral");
  
  OPTION(options, Lmax, 1.0); // Maximum mean free path [m]
  OPTION(options, frecycle, 0.9); // Recycling fraction

  Nn = 0.0; // initialise to 0
  Tn = 3.5 / Tnorm; // 3.5eV  for Franck-Condon energy
  
  S = 0;
  F = 0;
  Qi = 0;
  Rp = 0;
  SAVE_REPEAT4(Nn0, lambda, lambda_int, Nn);
  
  GRID_LOAD(hthe);
}

void NeutralRecycling::update(const Field3D &Ne, const Field3D &Te, const Field3D &Ti, const Field3D &Vi) {
  TRACE("NeutralRecycling::update");
  
  // Lower limit for neutral density (mainly for first time through when Nn = 0)
  Nn = floor(Nn, 1e-8);
  
  Field3D Nelim = floor(Ne, 1e-19); // Smaller limit for rate coefficients
  
  // Calculate scaling of Nn exponential (such that n_lost = n_recy on each field line)
  BoutReal nnexp, nlost, nnexp_tot;
  BoutReal vth_n, sigma_cx, sigma_iz, fluxout;
  BoutReal Nn0max = 10.0;//max(Nelim,true);
  lambda = 0.0;
  Nn0 = 0.0;
  // Approximate neutral density at t=0 to be exponential away from plate with max density equal to ion max density
  static bool first_time = true;
  if(first_time) {
    Field2D ll;
    ll = CumSumY2D(hthe*mesh->dy/Lmax, true);
    Nn = max(Nelim) * exp(-ll);
    first_time = false;
  }
  // calculate iz and cx mean free paths
  for(int i=mesh->xstart;i<=mesh->xend;i++) {
    for(int j=mesh->ystart;j<=mesh->yend;j++) {
      for(int k=0;k<mesh->ngz;k++) {
        vth_n = sqrt(Tn(i,j,k)) * Lnorm * Fnorm; // in m/s
        sigma_cx = Nelim(i,j,k) * Nnorm * hydrogen.chargeExchange(Te(i,j,k)*Tnorm); // in seconds
        sigma_iz = Nelim(i,j,k) * Nnorm * Nn(i,j,k) * hydrogen.ionisation(Te(i,j,k)*Tnorm);  // in seconds
        // mean-free path for cx and iz is sqrt(l_cx*l_iz), and each are vth/<sig-v>
        lambda(i,j,k) = vth_n / sqrt(sigma_cx*sigma_iz);//min( vth_n / sqrt(sigma_cx*sigma_iz), Lmax); // in meters
      }
    }
  }
  
  // Sum in y dy/lambda (in meters) for exponential (if lambda is constant, results in exp(-y/lambda) )
  lambda_int = CumSumY3D(hthe*mesh->dy/lambda, true); // length in poloidal plane - appropriate for neutrals
  
  // int(S dV) = fr * N_lost (ie. volume integral of ionization density source = recycling fraction * particle loss rate)
  for(int k=0;k<mesh->ngz;k++) {
    for(int i=mesh->xstart;i<=mesh->xend;i++) {
      // plasma density flux out to divertor plate
      fluxout = 0.25 * (Ne(i,mesh->yend,k) + Ne(i,mesh->yend+1,k))*(Vi(i,mesh->yend,k) + Vi(i,mesh->yend+1,k)); // [d^-2][t^-1]
      
      // number of particles lost per dt (flux out times perp volume) [t^-1]
      nlost = bcast_lasty( fluxout * 0.5 * (mesh->J(i,mesh->yend) * mesh->dx(i,mesh->yend) * mesh->dz / sqrt(mesh->g_22(i,mesh->yend)) 
                                            + mesh->J(i,mesh->yend+1) * mesh->dx(i,mesh->yend+1) * mesh->dz / sqrt(mesh->g_22(i,mesh->yend+1))) );
      
      // Integrate ionization rate over volume to get volume loss rate (simple integration using trap rule)
      nnexp = 0.0;
      for(int j=mesh->ystart;j<=mesh->yend;j++) {
        sigma_iz = hydrogen.ionisation(Te(i,j,k)*Tnorm) * Nnorm / Fnorm; // ionization rate [d^3]/[t]
        BoutReal dV = mesh->J(i,j) * mesh->dx(i,j) * mesh->dy(i,j) * mesh->dz; // volume element
        nnexp += Nelim(i,j,k) * sigma_iz * exp(-lambda_int(i,j,k)) * dV; // full integral of density source [d^3]/[t]
      }
      MPI_Allreduce(&nnexp, &nnexp_tot, 1, MPI_DOUBLE, MPI_SUM, mesh->getYcomm(i)); // add up all y-procs
      
      // neutral density factor
      for(int j=0;j<mesh->ngy;j++) {
        Nn0(i,j,k) = frecycle * nlost / nnexp_tot; // ( [d^-3] ) Max neutral density can't exceed max plasma density
      }
    }
  }
  Nn0 = -floor(-Nn0,-max(Ne,true));
  
  // Approximate neutral density as exponential
  Nn = Nn0 * exp(-lambda_int);
  
  // Calculate neutral density and source terms
  for(int i=mesh->xstart;i<=mesh->xend;i++)
    for(int j=mesh->ystart;j<=mesh->yend;j++)
      for(int k=0;k<mesh->ngz;k++) {
        // Rates
        BoutReal R_rc = SQ(Nelim(i,j,k)) * hydrogen.recombination(Nelim(i,j,k)*Nnorm, Te(i,j,k)*Tnorm) * Nnorm / Fnorm; // Time normalisation
        BoutReal R_iz = Nelim(i,j,k)*Nn(i,j,k) * hydrogen.ionisation(Te(i,j,k)*Tnorm) * Nnorm / Fnorm; // Time normalisation
        BoutReal R_cx = Nelim(i,j,k)*Nn(i,j,k) * hydrogen.chargeExchange(Te(i,j,k)*Tnorm) * Nnorm / Fnorm;
        
        // Ionisation plasma source, recombination sink
        S(i,j,k) = R_rc - R_iz;
        
        // Ionisation power transfer to plasma from neutrals
        Qi(i,j,k) = -Tn(i,j,k) * R_iz;
        // Recombination plasma power sink
        Qi(i,j,k) += Te(i,j,k) * R_rc;
        // Power transfer from plasma to neutrals due to charge exchange
        Qi(i,j,k) += R_cx * (Te(i,j,k) - Tn(i,j,k));
        
        // Ion-neutral friction due to charge exchange
        Fperp(i,j,k) = R_cx; 
        // Friction due to recombination
        Fperp(i,j,k) += R_rc;
        
        // Radiated power from plasma
        // Factor of 1.09 so that recombination becomes an energy source at 5.25eV
        Rp(i,j,k) = (1.09*Te(i,j,k) - 13.6/Tnorm)*R_rc
          + (Eionize/Tnorm) * R_iz; // Ionisation energy
      }
}


const Field2D NeutralRecycling::CumSumY2D(const Field2D &f, bool reverse) {
  // Cumulative sum in Y one xz-slice at a time
  //    -- reverse is option to sum from the end of Y
  Field2D result = 0.0;

  if(reverse) {
    for(int i=mesh->xstart;i<=mesh->xend;i++) {
      // All but last processor receive
      if(!mesh->lastY()) {
        mesh->wait(mesh->irecvYOutOutdest(&result(i,mesh->yend+1),1,mesh->ngx*i));
      }
      // Calculate sum (reversed)
      for(int j=mesh->yend;j>=mesh->ystart;j--) {
        result[i][j] = result[i][j+1] + f[i][j];
      }
      // Send the value at yend to the next processor.
      mesh->sendYInOutdest(&result(i,mesh->ystart),1,mesh->ngx*i);
    }
  } else{
    for(int i=mesh->xstart;i<=mesh->xend;i++) {
      // All but first processor receive
      if(!mesh->firstY()) {
        mesh->wait(mesh->irecvYInIndest(&result(i,mesh->ystart-1),1,mesh->ngx*i));
      }
      // Calculate sum
      for(int j=mesh->ystart;j<=mesh->yend;j++) {
        result[i][j] = result[i][j-1] + f[i][j];
      }
      // Send the value at yend to the next processor.
      mesh->sendYOutIndest(&result(i,mesh->yend),1,mesh->ngx*i);
    }
  }

  return result;
}

const Field3D NeutralRecycling::CumSumY3D(const Field3D &f, bool reverse) {
  // Cumulative sum in Y one xz-slice at a time
  //    -- reverse is option to sum from the end of Y
  Field3D result = 0.0;

  if(reverse) {
    for(int k=0;k<mesh->ngz;k++) {
      for(int i=mesh->xstart;i<=mesh->xend;i++) {
        // All but last processor receive
        if(!mesh->lastY()) {
          mesh->wait(mesh->irecvYOutOutdest(&result(i,mesh->yend+1,k),1,mesh->ngx*i+mesh->ngz*k));
        }
        // Calculate sum (reversed)
        for(int j=mesh->yend;j>=mesh->ystart;j--) {
          result[i][j][k] = result[i][j+1][k] + f[i][j][k];
        }
        // Send the value at yend to the next processor.
        mesh->sendYInOutdest(&result(i,mesh->ystart,k),1,mesh->ngx*i+mesh->ngz*k);
      }
    }
  } else{
    for(int k=0;k<mesh->ngz;k++) {
      for(int i=mesh->xstart;i<=mesh->xend;i++) {
        // All but first processor receive
        if(!mesh->firstY()) {
          mesh->wait(mesh->irecvYInIndest(&result(i,mesh->ystart-1,k),1,mesh->ngx*i+mesh->ngz*k));
        }
        // Calculate sum
        for(int j=mesh->ystart;j<=mesh->yend;j++) {
          result[i][j][k] = result[i][j-1][k] + f[i][j][k];
        }
        // Send the value at yend to the next processor.
        mesh->sendYOutIndest(&result(i,mesh->yend,k),1,mesh->ngx*i+mesh->ngz*k);
      }
    }
  }

  return result;
}

const BoutReal NeutralRecycling::bcast_lasty(const BoutReal f) {
  BoutReal myf;
  // All but last processor receive
  if(!mesh->lastY()) {
    mesh->wait(mesh->irecvYOutOutdest(&myf,1,0));
  } else {
    myf = f;
  }
  // Send the value at yend to the next processor.
  mesh->sendYInOutdest(&myf,1,0);

  return myf;
}
