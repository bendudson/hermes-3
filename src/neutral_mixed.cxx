
#include <bout/constants.hxx>
#include <bout/derivs.hxx>
#include <bout/difops.hxx>
#include <bout/fv_ops.hxx>
#include <bout/output_bout_types.hxx>

#include "../include/hermes_utils.hxx"
#include "../include/div_ops.hxx"
#include "../include/hermes_build_config.hxx"
#include "../include/neutral_mixed.hxx"

using bout::globals::mesh;

using ParLimiter = FV::Upwind;

NeutralMixed::NeutralMixed(const std::string& name, Options& alloptions, Solver* solver)
    : name(name) {
  AUTO_TRACE();

  // Normalisations
  const Options& units = alloptions["units"];
  const BoutReal meters = units["meters"];
  const BoutReal seconds = units["seconds"];
  const BoutReal Nnorm = units["inv_meters_cubed"];
  const BoutReal Tnorm = units["eV"];
  const BoutReal Omega_ci = 1. / units["seconds"].as<BoutReal>();

  // Need to take derivatives in X for cross-field diffusion terms
  ASSERT0(mesh->xstart > 0);

  auto& options = alloptions[name];

  // Evolving variables e.g name is "h" or "h+"
  solver->add(Nn, std::string("N") + name);
  solver->add(Pn, std::string("P") + name);

  evolve_momentum = options["evolve_momentum"]
                        .doc("Evolve parallel neutral momentum?")
                        .withDefault<bool>(true);

  if (evolve_momentum) {
    solver->add(NVn, std::string("NV") + name);
  } else {
    output_warn.write(
        "WARNING: Not evolving neutral parallel momentum. NVn and Vn set to zero\n");
    NVn = 0.0;
    Vn = 0.0;
  }

  sheath_ydown = options["sheath_ydown"]
                     .doc("Enable wall boundary conditions at ydown")
                     .withDefault<bool>(true);

  sheath_yup = options["sheath_yup"]
                   .doc("Enable wall boundary conditions at yup")
                   .withDefault<bool>(true);

  nn_floor = options["nn_floor"]
                 .doc("A minimum density used when dividing NVn by Nn. "
                      "Normalised units.")
                 .withDefault(1e-8);

  pn_floor = nn_floor * (1./get<BoutReal>(alloptions["units"]["eV"]));

  precondition = options["precondition"]
                     .doc("Enable preconditioning in neutral model?")
                     .withDefault<bool>(true);

  lax_flux = options["lax_flux"]
                     .doc("Enable stabilising lax flux?")
                     .withDefault<bool>(true);

  // flux_limit =
  //     options["flux_limit"]
  //         .doc("Limit diffusive fluxes to fraction of thermal speed. <0 means off.")
  //         .withDefault(0.2);

  advection_limit_alpha =
      options["advection_limit_alpha"]
          .doc("Limit perpendicular advection fluxes to fraction of thermal speed. <0 means off.")
          .withDefault(1.0);

  conduction_limit_alpha = options["conduction_limit_alpha"]
    .doc("Scale heat flux limiter")
    .withDefault(1.0);

  viscosity_limit_alpha = options["viscosity_limit_alpha"]
    .doc("Scale momentum flux limiter")
    .withDefault(1.0);

  flux_limit_gamma =
      options["flux_limit_gamma"]
          .doc("Sharpness of flux limiter. 1 is very loose, 2 loose and 5 reasonably tight")
          .withDefault(2);

  override_limiter = options["override_limiter"]
                     .doc("Force conduction and viscosity limiters to use the advection limiter?")
                     .withDefault<bool>(false);

  legacy_limiter = options["legacy_limiter"]
                     .doc("Use old form of the flux limiter?")
                     .withDefault<bool>(false);

  legacy_separate_conduction = options["legacy_separate_conduction"]
                     .doc("Feature separate conduction limiter when legacy form enabled? Requires legacy_limiter=True.")
                     .withDefault<bool>(false);

  maximum_mfp =
      options["maximum_mfp"]
          .doc("Add a pseudo-collisionality representing physical MFP limit to pressure diffusion model")
          .withDefault(0.1) 
                    / get<BoutReal>(alloptions["units"]["meters"]);

  diffusion_limit = options["diffusion_limit"]
                        .doc("Upper limit on diffusion coefficient [m^2/s]. <0 means off")
                        .withDefault(-1.0)
                    / (meters * meters / seconds); // Normalise

  neutral_viscosity = options["neutral_viscosity"]
                          .doc("Include neutral gas viscosity?")
                          .withDefault<bool>(true);

  neutral_conduction = options["neutral_conduction"]
                          .doc("Include neutral gas heat conduction?")
                          .withDefault<bool>(true);

  legacy_vth_limiter = options["legacy_vth_limiter"]
                          .doc("Use old formulation for Vth in limiter?")
                          .withDefault<bool>(true);

  diffusion_collisions_mode = options["diffusion_collisions_mode"]
      .doc("Can be legacy: all enabled collisions excl. IZ, or afn: CX, IZ and NN collisions")
      .withDefault<std::string>("legacy");

  if (precondition) {
    inv = std::unique_ptr<Laplacian>(Laplacian::create(&options["precon_laplace"]));

    inv->setInnerBoundaryFlags(INVERT_DC_GRAD | INVERT_AC_GRAD);
    inv->setOuterBoundaryFlags(INVERT_DC_GRAD | INVERT_AC_GRAD);

    inv->setCoefA(1.0);
  }

  // Optionally output time derivatives
  output_ddt =
      options["output_ddt"].doc("Save derivatives to output?").withDefault<bool>(false);

  diagnose =
      options["diagnose"].doc("Save additional diagnostics?").withDefault<bool>(false);

  AA = options["AA"].doc("Particle atomic mass. Proton = 1").withDefault(1.0);

  // Try to read the density source from the mesh
  // Units of particles per cubic meter per second
  density_source = 0.0;
  mesh->get(density_source, std::string("N") + name + "_src");
  // Allow the user to override the source
  density_source =
      alloptions[std::string("N") + name]["source"]
          .doc("Source term in ddt(N" + name + std::string("). Units [m^-3/s]"))
          .withDefault(density_source)
      / (Nnorm * Omega_ci);

  // Try to read the pressure source from the mesh
  // Units of Pascals per second
  pressure_source = 0.0;
  mesh->get(pressure_source, std::string("P") + name + "_src");
  // Allow the user to override the source
  pressure_source = alloptions[std::string("P") + name]["source"]
                        .doc(std::string("Source term in ddt(P") + name
                             + std::string("). Units [N/m^2/s]"))
                        .withDefault(pressure_source)
                    / (SI::qe * Nnorm * Tnorm * Omega_ci);

  // Set boundary condition defaults: Neumann for all but the diffusivity.
  // The dirichlet on diffusivity ensures no radial flux.
  // NV and V are ignored as they are hardcoded in the parallel BC code.
  alloptions[std::string("Dnn") + name]["bndry_all"] =
      alloptions[std::string("Dnn") + name]["bndry_all"].withDefault("dirichlet");
  alloptions[std::string("T") + name]["bndry_all"] =
      alloptions[std::string("T") + name]["bndry_all"].withDefault("neumann");
  alloptions[std::string("P") + name]["bndry_all"] =
      alloptions[std::string("P") + name]["bndry_all"].withDefault("neumann");
  alloptions[std::string("N") + name]["bndry_all"] =
      alloptions[std::string("N") + name]["bndry_all"].withDefault("neumann");

  // Pick up BCs from input file
  Dnn.setBoundary(std::string("Dnn") + name);
  Tn.setBoundary(std::string("T") + name);
  Pn.setBoundary(std::string("P") + name);
  Nn.setBoundary(std::string("N") + name);

  // All floored versions of variables get the same boundary as the original
  Tnlim.setBoundary(std::string("T") + name);
  Pnlim.setBoundary(std::string("P") + name);
  logPnlim.setBoundary(std::string("P") + name);
  Nnlim.setBoundary(std::string("N") + name);

  // Product of Dnn and another parameter has same BC as Dnn - see eqns to see why this is
  // necessary
  DnnNn.setBoundary(std::string("Dnn") + name);
  DnnPn.setBoundary(std::string("Dnn") + name);
  DnnNVn.setBoundary(std::string("Dnn") + name);
}

void NeutralMixed::transform(Options& state) {
  AUTO_TRACE();

  mesh->communicate(Nn, Pn, NVn);

  Nn.clearParallelSlices();
  Pn.clearParallelSlices();
  NVn.clearParallelSlices();

  Nn = floor(Nn, 0.0);
  Pn = floor(Pn, 0.0);

  // Nnlim Used where division by neutral density is needed
  Nnlim = floor(Nn, nn_floor);
  Tn = Pn / Nnlim;
  Tn.applyBoundary();

  Vn = NVn / (AA * Nnlim);
  Vnlim = Vn;

  Vn.applyBoundary("neumann");
  Vnlim.applyBoundary("neumann");

  Pnlim = floor(Pn, pn_floor);
  Pnlim.applyBoundary();

  Tnlim = Pnlim / Nnlim;

  /////////////////////////////////////////////////////
  // Parallel boundary conditions
  TRACE("Neutral boundary conditions");

  if (sheath_ydown) {
    for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        // Free boundary (constant gradient) density
        BoutReal nnwall =
            0.5 * (3. * Nn(r.ind, mesh->ystart, jz) - Nn(r.ind, mesh->ystart + 1, jz));
        if (nnwall < 0.0)
          nnwall = 0.0;

        BoutReal tnwall = Tn(r.ind, mesh->ystart, jz);

        Nn(r.ind, mesh->ystart - 1, jz) = 2 * nnwall - Nn(r.ind, mesh->ystart, jz);

        // Zero gradient temperature, heat flux added later
        Tn(r.ind, mesh->ystart - 1, jz) = tnwall;

        // Set pressure consistent at the boundary
        // Pn(r.ind, mesh->ystart - 1, jz) =
        //     2. * nnwall * tnwall - Pn(r.ind, mesh->ystart, jz);

        // Zero-gradient pressure
        Pn(r.ind, mesh->ystart - 1, jz) = Pn(r.ind, mesh->ystart, jz);
        Pnlim(r.ind, mesh->ystart - 1, jz) = Pnlim(r.ind, mesh->ystart, jz);

        // No flow into wall
        Vn(r.ind, mesh->ystart - 1, jz) = -Vn(r.ind, mesh->ystart, jz);
        Vnlim(r.ind, mesh->ystart - 1, jz) = -Vnlim(r.ind, mesh->ystart, jz);
        NVn(r.ind, mesh->ystart - 1, jz) = -NVn(r.ind, mesh->ystart, jz);
      }
    }
  }

  if (sheath_yup) {
    for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        // Free boundary (constant gradient) density
        BoutReal nnwall =
            0.5 * (3. * Nn(r.ind, mesh->yend, jz) - Nn(r.ind, mesh->yend - 1, jz));
        if (nnwall < 0.0)
          nnwall = 0.0;

        BoutReal tnwall = Tn(r.ind, mesh->yend, jz);

        Nn(r.ind, mesh->yend + 1, jz) = 2 * nnwall - Nn(r.ind, mesh->yend, jz);

        // Zero gradient temperature, heat flux added later
        Tn(r.ind, mesh->yend + 1, jz) = tnwall;

        // Zero-gradient pressure
        Pn(r.ind, mesh->yend + 1, jz) = Pn(r.ind, mesh->yend, jz);
        Pnlim(r.ind, mesh->yend + 1, jz) = Pnlim(r.ind, mesh->yend, jz);

        // No flow into wall
        Vn(r.ind, mesh->yend + 1, jz) = -Vn(r.ind, mesh->yend, jz);
        Vnlim(r.ind, mesh->yend + 1, jz) = -Vnlim(r.ind, mesh->yend, jz);
        NVn(r.ind, mesh->yend + 1, jz) = -NVn(r.ind, mesh->yend, jz);
      }
    }
  }

  // Set values in the state
  auto& localstate = state["species"][name];
  set(localstate["density"], Nn);
  set(localstate["AA"], AA); // Atomic mass
  set(localstate["pressure"], Pn);
  set(localstate["momentum"], NVn);
  set(localstate["velocity"], Vn);
  set(localstate["temperature"], Tn);
}

void NeutralMixed::finally(const Options& state) {
  AUTO_TRACE();
  auto& localstate = state["species"][name];

  // Logarithms used to calculate perpendicular velocity
  // V_perp = -Dnn * ( Grad_perp(Nn)/Nn + Grad_perp(Tn)/Tn )
  //
  // Grad(Pn) / Pn = Grad(Tn)/Tn + Grad(Nn)/Nn
  //               = Grad(logTn + logNn)
  // Field3D logNn = log(Nn);
  // Field3D logTn = log(Tn);

  logPnlim = log(Pnlim);
  logPnlim.applyBoundary();

  ///////////////////////////////////////////////////////
  // Calculate cross-field diffusion from collision frequency
  //
  //

  Field3D mfp_pseudo_nu =
    sqrt(floor(Tn, 1e-5) / AA) / maximum_mfp; // Pseudo-collisionality due to vessel size MFP limit

  if (localstate.isSet("collision_frequency")) {

    // Collisionality
    // Braginskii mode: plasma - self collisions and ei, neutrals - CX, IZ
    if (collision_names.empty()) {     /// Calculate only once - at the beginning

      if (diffusion_collisions_mode == "afn") {
        for (const auto& collision : localstate["collision_frequencies"].getChildren()) {

          std::string collision_name = collision.second.name();

          if (/// Charge exchange
              (collisionSpeciesMatch(    
                collision_name, name, "+", "cx", "partial")) or
              /// Ionisation
              (collisionSpeciesMatch(    
                collision_name, name, "+", "iz", "partial")) or
              /// Neutral-neutral collisions
              (collisionSpeciesMatch(    
                collision_name, name, name, "coll", "exact"))) {
                  collision_names.push_back(collision_name);
                }
        }
      // Legacy mode: all collisions and CX are included
      } else if (diffusion_collisions_mode == "legacy") {
        for (const auto& collision : localstate["collision_frequencies"].getChildren()) {

          std::string collision_name = collision.second.name();

          if (/// Charge exchange
              (collisionSpeciesMatch(    
                collision_name, name, "", "cx", "partial")) or
              /// Any collision (en, in, ee, ii, nn)
              (collisionSpeciesMatch(    
                collision_name, name, "", "coll", "partial"))) {
                  collision_names.push_back(collision_name);
                }
        }
        
      } else {
        throw BoutException("\ndiffusion_collisions_mode for {:s} must be either legacy or braginskii", name);
      }

      /// Write chosen collisions to log file
      output_info.write("\t{:s} neutral collisionality mode: '{:s}' using ",
                      name, diffusion_collisions_mode);
      for (const auto& collision : collision_names) {        
        output_info.write("{:s} ", collision);
      }
      output_info.write("\n");
      }

    /// Collect the collisionalities based on list of names
    nu = 0;
    for (const auto& collision_name : collision_names) {
      nu += GET_VALUE(Field3D, localstate["collision_frequencies"][collision_name]);
    }

    // Dnn = Vth^2 / sigma
    // This thermal speed is isotropic, so sqrt(T/m)
    Dnn_unlimited = (Tn / AA) / (nu + mfp_pseudo_nu);
  } else {
    Dnn_unlimited = (Tn / AA) / mfp_pseudo_nu;
  }

  if (legacy_vth_limiter) {
    vth = sqrt(Tn / AA);
  } else {
    vth = 0.25 * sqrt((8 * Tn) / (PI * AA));
  }

  // Legacy flux limiter: limit Dn upstream
  Dnn = 0;
  if (legacy_limiter) {
    Dmax = advection_limit_alpha * vth / (abs(Grad_perp(logPnlim)) + 1. / maximum_mfp);
    BOUT_FOR(i, Dmax.getRegion("RGN_NOBNDRY")) { Dnn[i] = BOUTMIN(Dnn_unlimited[i], Dmax[i]); }
  } else {
    Dmax = Dnn_unlimited;
    Dnn = Dnn_unlimited;
  }

  if (diffusion_limit > 0.0) {
    // Impose an upper limit on the diffusion coefficient
    BOUT_FOR(i, Dnn.getRegion("RGN_NOBNDRY")) {
      Dnn[i] = BOUTMIN(Dnn[i], diffusion_limit);
    }
  }

  mesh->communicate(Dnn);
  Dnn.clearParallelSlices();
  Dnn.applyBoundary();

  // Neutral diffusion parameters have the same boundary condition as Dnn
  DnnNn = Dnn * Nnlim;
  DnnPn = Dnn * Pnlim;
  DnnNVn = Dnn * NVn;

  DnnPn.applyBoundary();
  DnnNn.applyBoundary();
  DnnNVn.applyBoundary();

  // Heat conductivity 
  // Note: This is kappa_n = (5/2) * Pn / (m * nu)
  //       where nu is the collision frequency used in Dnn

  // Implement separate conduction limiter by limiting Kappa
  // Just like the legacy limiters limit D.
  // Note that kappa is always calculated using unlimited D.

  Field3D DnnNn_unlimited = Dnn_unlimited * Nnlim;
  DnnNn_unlimited.applyBoundary("dirichlet");    // TODO: is this correct?
  kappa_n_unlimited = (5. / 2) * DnnNn_unlimited;
  kappa_n_Dnchained = (5. / 2) * DnnNn;   // Include only limited D, not also limited kappa
  
  if (legacy_limiter and legacy_separate_conduction) {
    Field3D cond_vel = Pnlim * sqrt((2*Tnlim) / (PI*AA));  // 1D heat flux of 3D maxwellian (Stangeby)  
    kappa_n_max = conduction_limit_alpha * cond_vel / (abs(Grad_perp(Tn)) + 1. / maximum_mfp);
    BOUT_FOR(i, kappa_n_max.getRegion("RGN_NOBNDRY")) { kappa_n[i] = BOUTMIN(kappa_n_unlimited[i], kappa_n_max[i]); }

  } else {
    kappa_n = (5. / 2) * DnnNn;
    kappa_n_max = 0;
  }
    

  // Viscosity
  // Relationship between heat conduction and viscosity for neutral
  // gas Chapman, Cowling "The Mathematical Theory of Non-Uniform
  // Gases", CUP 1952 Ferziger, Kaper "Mathematical Theory of
  // Transport Processes in Gases", 1972
  // eta_n = (2. / 5) * m_n * kappa_n;
  //
  eta_n = AA * (2. / 5) * kappa_n_Dnchained;

  // These are for debugging only
  gradlogP = abs(Grad(logPnlim));
  gradperplogP = abs(Grad_perp(logPnlim));


  //// Calculate flux limiting factors for perpendicular transport only
  advection_factor = 1;
  conduction_factor = 1;
  viscosity_factor = 1;

  // Advection (of particles, pressure, momentum)
  if (advection_limit_alpha > 0.0) {
    Vector3D v_perp = -Dnn * Grad_perp(logPnlim);     // vector of perp velocity
    Field3D v_abs = sqrt(v_perp * v_perp);            // magintude: |v dot v|
    advection_flux_abs = Nnlim * v_abs;
    advection_limit = Nnlim * vth;          
    advection_factor = pow(1. + pow(advection_flux_abs / (advection_limit_alpha * advection_limit),
                                          flux_limit_gamma),-1./flux_limit_gamma);
  } else {
    advection_factor = 1;
  }

  // Conduction
  if (conduction_limit_alpha > 0.0 and neutral_conduction) {
    Vector3D heat_flux = -kappa_n * Grad_perp(Tn);  
    Field3D heat_flux_abs = sqrt(heat_flux * heat_flux);
    Field3D heat_limit = Pnlim * sqrt((2*Tnlim) / (PI*AA));  // 1D heat flux of 3D maxwellian (Stangeby)          
    conduction_factor = pow(1. + pow(heat_flux_abs / (conduction_limit_alpha * heat_limit),
                                          flux_limit_gamma),-1./flux_limit_gamma);
  } else {
    conduction_factor = 1;
  }

  // Viscosity
  if (viscosity_limit_alpha > 0.0 and neutral_viscosity) {
    Vector3D momentum_flux = -eta_n * Grad_perp(Vn);     
    Field3D momentum_flux_abs = sqrt(momentum_flux * momentum_flux);
    Field3D momentum_limit = Pnlim;    // Can't have more dynamic pressure than there is static pressure
    viscosity_factor = pow(1. + pow(momentum_flux_abs / (viscosity_limit_alpha * momentum_limit),
                                          flux_limit_gamma),-1./flux_limit_gamma);
  } else {
    viscosity_factor = 1;
  }

  // Force conduction/viscosity to use advection limiter
  if (override_limiter) {
    conduction_factor = advection_factor;
    viscosity_factor = advection_factor;
  }

  // Set all flux factors to 1, use limited Dnn instead (see upstream)
  if (legacy_limiter) {
    advection_factor = 1;
    conduction_factor = 1;
    viscosity_factor = 1;
  }

  

  if (sheath_ydown) {
    for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        Dnn(r.ind, mesh->ystart - 1, jz) = -Dnn(r.ind, mesh->ystart, jz);
        DnnNn(r.ind, mesh->ystart - 1, jz) = -DnnNn(r.ind, mesh->ystart, jz);
        DnnPn(r.ind, mesh->ystart - 1, jz) = -DnnPn(r.ind, mesh->ystart, jz);
        DnnNVn(r.ind, mesh->ystart - 1, jz) = -DnnNVn(r.ind, mesh->ystart, jz);
      }
    }
  }

  if (sheath_yup) {
    for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        Dnn(r.ind, mesh->yend + 1, jz) = -Dnn(r.ind, mesh->yend, jz);
        DnnNn(r.ind, mesh->yend + 1, jz) = -DnnNn(r.ind, mesh->yend, jz);
        DnnPn(r.ind, mesh->yend + 1, jz) = -DnnPn(r.ind, mesh->yend, jz);
        DnnNVn(r.ind, mesh->yend + 1, jz) = -DnnNVn(r.ind, mesh->yend, jz);
      }
    }
  }

  // Sound speed appearing in Lax flux for advection terms
  sound_speed = 0;
  if (lax_flux) {
    sound_speed = sqrt(Tn * (5. / 3) / AA);
  }



  /////////////////////////////////////////////////////
  // Neutral density
  TRACE("Neutral density");

  perp_nn_adv_src = Div_a_Grad_perp_upwind_flows(DnnNn * advection_factor, logPnlim,       // Perpendicular advection
                                   particle_flow_xlow,
                                   particle_flow_ylow);                 

  par_nn_adv_src = FV::Div_par_mod<ParLimiter>(Nn, Vn, sound_speed);    // Parallel advection

  ddt(Nn) =
    - par_nn_adv_src
    + perp_nn_adv_src
    ;

  Sn = density_source; // Save for possible output
  if (localstate.isSet("density_source")) {
    Sn += get<Field3D>(localstate["density_source"]);
  }
  ddt(Nn) += Sn; // Always add density_source

  /////////////////////////////////////////////////////
  // Neutral pressure
  TRACE("Neutral pressure");

  ddt(Pn) = 
    - FV::Div_par_mod<ParLimiter>(Pn, Vn, sound_speed)                  // Parallel advection
    - (2. / 3) * Pn * Div_par(Vn)                                       // Parallel compression
    + Div_a_Grad_perp_upwind_flows(
          (5. / 3) * DnnPn * advection_factor, logPnlim,                // Perpendicular advection
          energy_flow_xlow, energy_flow_ylow)  
     ;
  // The factor here is likely 5/2 as we're advecting internal energy and pressure.
  // Doing this still leaves a heat imbalance factor of 0.11 in the cells, but better than 0.33 with 3/2.
  energy_flow_xlow *= 5/2; 
  energy_flow_ylow *= 5/2;

  gradperpT = sqrt(Grad_perp(Tn) * Grad_perp(Tn));

  if (neutral_conduction) {
    ddt(Pn) += 
      (2. / 3) * Div_a_Grad_perp_upwind_flows(kappa_n * conduction_factor, Tn,
                            conduction_flow_xlow, conduction_flow_ylow)                      // Perpendicular conduction
      + FV::Div_par_K_Grad_par(kappa_n * conduction_factor, Tn)                             // Parallel conduction
      ;

    // The factor here is likely 3/2 as this is pure energy flow, but needs checking.
    conduction_flow_xlow *= 3/2;
    conduction_flow_ylow *= 3/2;
  }

  
  
  Sp = pressure_source;
  if (localstate.isSet("energy_source")) {
    Sp += (2. / 3) * get<Field3D>(localstate["energy_source"]);
  }
  ddt(Pn) += Sp;


  if (evolve_momentum) {

    /////////////////////////////////////////////////////
    // Neutral momentum
    TRACE("Neutral momentum");

    ddt(NVn) =
        -AA * FV::Div_par_fvv<ParLimiter>(Nnlim, Vn, sound_speed)       // Parallel advection
        - Grad_par(Pn)                                                  // Pressure gradient
      + Div_a_Grad_perp_upwind_flows(
            DnnNVn * advection_factor, logPnlim,                    // Perpendicular advection
            momentum_flow_xlow,
            momentum_flow_ylow) 
      ;

    if (neutral_viscosity) {
      // NOTE: The following viscosity terms are not (yet) balanced
      //       by a viscous heating term

      Field3D momentum_source = FV::Div_a_Grad_perp(eta_n * viscosity_factor, Vn)  // Perpendicular viscosity
              + FV::Div_par_K_Grad_par(eta_n * viscosity_factor, Vn)               // Parallel viscosity
      ;

      ddt(NVn) += momentum_source; // Viscosity
      ddt(Pn) += -(2. /3) * Vn * momentum_source;                       // Viscous heating

    }

    if (localstate.isSet("momentum_source")) {
      Snv = get<Field3D>(localstate["momentum_source"]);
      ddt(NVn) += Snv;
    }

  } else {
    ddt(NVn) = 0;
    Snv = 0;
  }

  BOUT_FOR(i, Pn.getRegion("RGN_ALL")) {
    if ((Pn[i] < pn_floor * 1e-2) && (ddt(Pn)[i] < 0.0)) {
      ddt(Pn)[i] = 0.0;
    }
    if ((Nn[i] < nn_floor * 1e-2) && (ddt(Nn)[i] < 0.0)) {
      ddt(Nn)[i] = 0.0;
    }
  }

  // Scale time derivatives
  if (state.isSet("scale_timederivs")) {
    Field3D scale_timederivs = get<Field3D>(state["scale_timederivs"]);
    ddt(Nn) *= scale_timederivs;
    ddt(Pn) *= scale_timederivs;
    ddt(NVn) *= scale_timederivs;
  }

#if CHECKLEVEL >= 1
  for (auto& i : Nn.getRegion("RGN_NOBNDRY")) {
    if (!std::isfinite(ddt(Nn)[i])) {
      throw BoutException("ddt(N{}) non-finite at {}\n", name, i);
    }
    if (!std::isfinite(ddt(Pn)[i])) {
      throw BoutException("ddt(P{}) non-finite at {}\n", name, i);
    }
    if (!std::isfinite(ddt(NVn)[i])) {
      throw BoutException("ddt(NV{}) non-finite at {}\n", name, i);
    }
  }
#endif
}

void NeutralMixed::outputVars(Options& state) {
  // Normalisations
  auto Nnorm = get<BoutReal>(state["Nnorm"]);
  auto Tnorm = get<BoutReal>(state["Tnorm"]);
  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
  auto Cs0 = get<BoutReal>(state["Cs0"]);
  auto rho_s0 = get<BoutReal>(state["rho_s0"]);
  const BoutReal Pnorm = SI::qe * Tnorm * Nnorm;

  state[std::string("N") + name].setAttributes({{"time_dimension", "t"},
                                                {"units", "m^-3"},
                                                {"conversion", Nnorm},
                                                {"standard_name", "density"},
                                                {"long_name", name + " number density"},
                                                {"species", name},
                                                {"source", "neutral_mixed"}});

  state[std::string("P") + name].setAttributes({{"time_dimension", "t"},
                                                {"units", "Pa"},
                                                {"conversion", Pnorm},
                                                {"standard_name", "pressure"},
                                                {"long_name", name + " pressure"},
                                                {"species", name},
                                                {"source", "neutral_mixed"}});

  state[std::string("NV") + name].setAttributes(
      {{"time_dimension", "t"},
       {"units", "kg / m^2 / s"},
       {"conversion", SI::Mp * Nnorm * Cs0},
       {"standard_name", "momentum"},
       {"long_name", name + " parallel momentum"},
       {"species", name},
       {"source", "neutral_mixed"}});

  if (output_ddt) {
    set_with_attrs(
        state[std::string("ddt(N") + name + std::string(")")], ddt(Nn),
        {{"time_dimension", "t"},
         {"units", "m^-3 s^-1"},
         {"conversion", Nnorm * Omega_ci},
         {"long_name", std::string("Rate of change of ") + name + " number density"},
         {"source", "neutral_mixed"}});
    set_with_attrs(state[std::string("ddt(P") + name + std::string(")")], ddt(Pn),
                   {{"time_dimension", "t"},
                    {"units", "Pa s^-1"},
                    {"conversion", Pnorm * Omega_ci},
                    {"source", "neutral_mixed"}});
    set_with_attrs(state[std::string("ddt(NV") + name + std::string(")")], ddt(NVn),
                   {{"time_dimension", "t"},
                    {"units", "kg m^-2 s^-2"},
                    {"conversion", SI::Mp * Nnorm * Cs0 * Omega_ci},
                    {"source", "neutral_mixed"}});
  }
  if (diagnose) {
    set_with_attrs(state[std::string("T") + name], Tn,
                   {{"time_dimension", "t"},
                    {"units", "eV"},
                    {"conversion", Tnorm},
                    {"standard_name", "temperature"},
                    {"long_name", name + " temperature"},
                    {"source", "neutral_mixed"}});
    set_with_attrs(state[std::string("Dnn") + name], Dnn,
                   {{"time_dimension", "t"},
                    {"units", "m^2/s"},
                    {"conversion", Cs0 * Cs0 / Omega_ci},
                    {"standard_name", "diffusion coefficient"},
                    {"long_name", name + " diffusion coefficient"},
                    {"source", "neutral_mixed"}});
    set_with_attrs(state[std::string("Dnn_unlim") + name], Dnn_unlimited,
                   {{"time_dimension", "t"},
                    {"units", "m^2/s"},
                    {"conversion", Cs0 * Cs0 / Omega_ci},
                    {"standard_name", "unlimited diffusion coefficient"},
                    {"long_name", name + " unlimited diffusion coefficient"},
                    {"source", "neutral_mixed"}});
    set_with_attrs(state[std::string("Dmax_") + name], Dmax,
                   {{"time_dimension", "t"},
                    {"units", "m^2/s"},
                    {"conversion", Cs0 * Cs0 / Omega_ci},
                    {"standard_name", "max diffusion coefficient"},
                    {"long_name", name + " max diffusion coefficient"},
                    {"source", "neutral_mixed"}});
    set_with_attrs(state[std::string("kappa_") + name], kappa_n,
                   {{"time_dimension", "t"},
                    {"units", "m^2/s"},
                    {"conversion", Cs0 * Cs0 / Omega_ci},
                    {"standard_name", "conduction coefficient"},
                    {"long_name", name + " conduction coefficient"},
                    {"source", "neutral_mixed"}});
    set_with_attrs(state[std::string("DnnNn_") + name], DnnNn,
                   {{"time_dimension", "t"},
                    {"units", "m^2/s"},
                    {"conversion", Cs0 * Cs0 / Omega_ci * Nnorm},
                    {"standard_name", ""},
                    {"long_name", ""},
                    {"source", "neutral_mixed"}});
    set_with_attrs(state[std::string("gradlogP_") + name], gradlogP,
                   {{"time_dimension", "t"},
                    {"units", "m^-1"},
                    {"conversion", 1 / rho_s0},
                    {"standard_name", "inv. P gradient length scale"},
                    {"long_name", name + " inv. P gradient length scale"},
                    {"source", "neutral_mixed"}});
    set_with_attrs(state[std::string("gradperplogP_") + name], gradperplogP,
                   {{"time_dimension", "t"},
                    {"units", "m^-1"},
                    {"conversion", 1 / rho_s0},
                    {"standard_name", "inv. P perp gradient length scale"},
                    {"long_name", name + " inv. P perp gradient length scale"},
                    {"source", "neutral_mixed"}});
    set_with_attrs(state[std::string("gradperpT") + name], gradperpT,
                   {{"time_dimension", "t"},
                    {"units", "m^-1"},
                    {"conversion", 1 / rho_s0},
                    {"standard_name", "inv. T perp gradient length scale"},
                    {"long_name", name + " inv. T perp gradient length scale"},
                    {"source", "neutral_mixed"}});
    set_with_attrs(state[std::string("advection_factor_") + name], advection_factor,
                   {{"time_dimension", "t"},
                    {"units", ""},
                    {"conversion", 1.0},
                    {"standard_name", "flux factor"},
                    {"long_name", name + " particle flux factor"},
                    {"species", name},
                    {"source", "neutral_mixed"}});
    set_with_attrs(state[std::string("vth_") + name], vth,
                   {{"time_dimension", "t"},
                    {"units", "m / s"},
                    {"conversion", Cs0},
                    {"standard_name", "thermal speed"},
                    {"long_name", name + " thermal speed"},
                    {"species", name},
                    {"source", "neutral_mixed"}});
    set_with_attrs(state[std::string("advection_flux_abs_") + name], advection_flux_abs,
                   {{"time_dimension", "t"},
                    {"units", "m^-2 s^-1"},
                    {"conversion", Nnorm * Cs0},
                    {"standard_name", ""},
                    {"long_name", ""},
                    {"species", ""},
                    {"source", "neutral_mixed"}});
    set_with_attrs(state[std::string("advection_limit_") + name], advection_limit,
                   {{"time_dimension", "t"},
                    {"units", "m^-2 s^-1"},
                    {"conversion", Nnorm * Cs0},
                    {"standard_name", ""},
                    {"long_name", ""},
                    {"species", ""},
                    {"source", "neutral_mixed"}});
    set_with_attrs(state[std::string("conduction_factor_") + name], conduction_factor,
                   {{"time_dimension", "t"},
                    {"units", ""},
                    {"conversion", 1.0},
                    {"standard_name", "flux factor"},
                    {"long_name", name + " conduction_factor"},
                    {"species", name},
                    {"source", "neutral_mixed"}});
    set_with_attrs(state[std::string("viscosity_factor_") + name], viscosity_factor,
                   {{"time_dimension", "t"},
                    {"units", ""},
                    {"conversion", 1.0},
                    {"standard_name", "flux factor"},
                    {"long_name", name + " viscosity_factor"},
                    {"species", name},
                    {"source", "neutral_mixed"}});
    set_with_attrs(state[std::string("SN") + name], Sn,
                   {{"time_dimension", "t"},
                    {"units", "m^-3 s^-1"},
                    {"conversion", Nnorm * Omega_ci},
                    {"standard_name", "density source"},
                    {"long_name", name + " number density source"},
                    {"source", "neutral_mixed"}});
    set_with_attrs(state[std::string("SP") + name], Sp,
                   {{"time_dimension", "t"},
                    {"units", "Pa s^-1"},
                    {"conversion", SI::qe * Tnorm * Nnorm * Omega_ci},
                    {"standard_name", "pressure source"},
                    {"long_name", name + " pressure source"},
                    {"source", "neutral_mixed"}});
    set_with_attrs(state[std::string("SNV") + name], Snv,
                   {{"time_dimension", "t"},
                    {"units", "kg m^-2 s^-2"},
                    {"conversion", SI::Mp * Nnorm * Cs0 * Omega_ci},
                    {"standard_name", "momentum source"},
                    {"long_name", name + " momentum source"},
                    {"source", "neutral_mixed"}});
    set_with_attrs(state[std::string("S") + name + std::string("_src")], density_source,
                   {{"time_dimension", "t"},
                    {"units", "m^-3 s^-1"},
                    {"conversion", Nnorm * Omega_ci},
                    {"standard_name", "density source"},
                    {"long_name", name + " number density source"},
                    {"species", name},
                    {"source", "neutral_mixed"}});
    set_with_attrs(state[std::string("P") + name + std::string("_src")], pressure_source,
                   {{"time_dimension", "t"},
                    {"units", "Pa s^-1"},
                    {"conversion", Pnorm * Omega_ci},
                    {"standard_name", "pressure source"},
                    {"long_name", name + " pressure source"},
                    {"species", name},
                    {"source", "neutral_mixed"}});
    set_with_attrs(state[std::string("S") + name + std::string("_perp_adv")], perp_nn_adv_src,
                   {{"time_dimension", "t"},
                    {"units", "m^-3 s^-1"},
                    {"conversion", Nnorm * Omega_ci},
                    {"standard_name", "density source due to perp advection"},
                    {"long_name", name + " number density source due to perp advection"},
                    {"species", name},
                    {"source", "neutral_mixed"}});
    set_with_attrs(state[std::string("S") + name + std::string("_par_adv")], par_nn_adv_src,
                   {{"time_dimension", "t"},
                    {"units", "m^-3 s^-1"},
                    {"conversion", Nnorm * Omega_ci},
                    {"standard_name", "density source due to par advection"},
                    {"long_name", name + " number density source due to par advection"},
                    {"species", name},
                    {"source", "neutral_mixed"}});

    if (particle_flow_xlow.isAllocated()) {
      set_with_attrs(state[std::string("ParticleFlow_") + name + std::string("_xlow")], particle_flow_xlow,
                   {{"time_dimension", "t"},
                    {"units", "s^-1"},
                    {"conversion", rho_s0 * SQ(rho_s0) * Nnorm * Omega_ci},
                    {"standard_name", "particle flow"},
                    {"long_name", name + " particle flow in X. Note: May be incomplete."},
                    {"species", name},
                    {"source", "neutral_mixed"}});
    }
    if (particle_flow_ylow.isAllocated()) {
      set_with_attrs(state[std::string("ParticleFlow_") + name + std::string("_ylow")], particle_flow_ylow,
                   {{"time_dimension", "t"},
                    {"units", "s^-1"},
                    {"conversion", rho_s0 * SQ(rho_s0) * Nnorm * Omega_ci},
                    {"standard_name", "particle flow"},
                    {"long_name", name + " particle flow in Y. Note: May be incomplete."},
                    {"species", name},
                    {"source", "evolve_density"}});
    }
    if (momentum_flow_xlow.isAllocated()) {
      set_with_attrs(state[std::string("MomentumFlow_") + name + std::string("_xlow")], momentum_flow_xlow,
                   {{"time_dimension", "t"},
                    {"units", "N"},
                    {"conversion", rho_s0 * SQ(rho_s0) * SI::Mp * Nnorm * Cs0 * Omega_ci},
                    {"standard_name", "momentum flow"},
                    {"long_name", name + " momentum flow in X. Note: May be incomplete."},
                    {"species", name},
                    {"source", "evolve_momentum"}});
    }
    if (momentum_flow_ylow.isAllocated()) {
      set_with_attrs(state[std::string("MomentumFlow_") + name + std::string("_ylow")], momentum_flow_ylow,
                   {{"time_dimension", "t"},
                    {"units", "N"},
                    {"conversion", rho_s0 * SQ(rho_s0) * SI::Mp * Nnorm * Cs0 * Omega_ci},
                    {"standard_name", "momentum flow"},
                    {"long_name", name + " momentum flow in Y. Note: May be incomplete."},
                    {"species", name},
                    {"source", "evolve_momentum"}});
    }
    if (energy_flow_xlow.isAllocated()) {
      set_with_attrs(state[std::string("EnergyFlow_") + name + std::string("_xlow")], energy_flow_xlow,
                   {{"time_dimension", "t"},
                    {"units", "W"},
                    {"conversion", rho_s0 * SQ(rho_s0) * Pnorm * Omega_ci},
                    {"standard_name", "power"},
                    {"long_name", name + " power through X cell face. Note: May be incomplete."},
                    {"species", name},
                    {"source", "evolve_pressure"}});
    }
    if (energy_flow_ylow.isAllocated()) {
      set_with_attrs(state[std::string("EnergyFlow_") + name + std::string("_ylow")], energy_flow_ylow,
                   {{"time_dimension", "t"},
                    {"units", "W"},
                    {"conversion", rho_s0 * SQ(rho_s0) * Pnorm * Omega_ci},
                    {"standard_name", "power"},
                    {"long_name", name + " power through Y cell face. Note: May be incomplete."},
                    {"species", name},
                    {"source", "evolve_pressure"}});
    }
    if (conduction_flow_xlow.isAllocated()) {
      set_with_attrs(state[std::string("ConductionFlow_") + name + std::string("_xlow")],conduction_flow_xlow,
                   {{"time_dimension", "t"},
                    {"units", "W"},
                    {"conversion", rho_s0 * SQ(rho_s0) * Pnorm * Omega_ci},
                    {"standard_name", "power"},
                    {"long_name", name + " conducted power through X cell face. Note: May be incomplete."},
                    {"species", name},
                    {"source", "evolve_pressure"}});
    }
    if (conduction_flow_ylow.isAllocated()) {
      set_with_attrs(state[std::string("ConductionFlow_") + name + std::string("_ylow")], conduction_flow_ylow,
                   {{"time_dimension", "t"},
                    {"units", "W"},
                    {"conversion", rho_s0 * SQ(rho_s0) * Pnorm * Omega_ci},
                    {"standard_name", "power"},
                    {"long_name", name + " conducted power through Y cell face. Note: May be incomplete."},
                    {"species", name},
                    {"source", "evolve_pressure"}});
    }
  }
}

void NeutralMixed::precon(const Options& state, BoutReal gamma) {
  if (!precondition) {
    return;
  }

  // Neutral gas diffusion
  // Solve (1 - gamma*Dnn*Delp2)^{-1}

  Field3D coef = -gamma * Dnn;

  if (state.isSet("scale_timederivs")) {
    coef *= get<Field3D>(state["scale_timederivs"]);
  }

  inv->setCoefD(coef);

  ddt(Nn) = inv->solve(ddt(Nn));
  if (evolve_momentum) {
    ddt(NVn) = inv->solve(ddt(NVn));
  }
  ddt(Pn) = inv->solve(ddt(Pn));
}
