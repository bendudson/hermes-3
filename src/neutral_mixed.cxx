
#include <bout/constants.hxx>
#include <bout/fv_ops.hxx>
#include <bout/derivs.hxx>
#include <bout/difops.hxx>
#include <bout/output_bout_types.hxx>

#include "../include/div_ops.hxx"
#include "../include/neutral_mixed.hxx"
#include "../include/hermes_build_config.hxx"

using bout::globals::mesh;

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
  solver->add(NVn, std::string("NV") + name);

  sheath_ydown = options["sheath_ydown"]
                     .doc("Enable wall boundary conditions at ydown")
                     .withDefault<bool>(true);

  sheath_yup = options["sheath_yup"]
                   .doc("Enable wall boundary conditions at yup")
                   .withDefault<bool>(true);

  nn_floor = options["nn_floor"]
                 .doc("A minimum density used when dividing NVn by Nn. "
                      "Normalised units.")
                 .withDefault(1e-5);

  precondition = options["precondition"]
                     .doc("Enable preconditioning in neutral model?")
                     .withDefault<bool>(true);

  flux_limit = options["flux_limit"]
    .doc("Use isotropic flux limiters?")
    .withDefault(true);

  particle_flux_limiter = options["particle_flux_limiter"]
    .doc("Enable particle flux limiter?")
    .withDefault(true);

  heat_flux_limiter = options["heat_flux_limiter"]
    .doc("Enable heat flux limiter?")
    .withDefault(true);

  momentum_flux_limiter = options["momentum_flux_limiter"]
    .doc("Enable momentum flux limiter?")
    .withDefault(true);

  flux_limit_alpha = options["flux_limit_alpha"]
    .doc("Scale flux limits")
    .withDefault(1.0);

  flux_limit_gamma = options["flux_limit_gamma"]
    .doc("Higher values increase sharpness of flux limiting")
    .withDefault(2.0);

  neutral_viscosity = options["neutral_viscosity"]
    .doc("Include neutral gas viscosity?")
    .withDefault<bool>(true);

  maximum_mfp = options["maximum_mfp"]
    .doc("Optional maximum mean free path in [m] for diffusive processes. < 0 is off")
    .withDefault(0.1);

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
  density_source = alloptions[std::string("N") + name]["source"]
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
  alloptions[std::string("Dnn") + name]["bndry_all"] = alloptions[std::string("Dnn") + name]["bndry_all"].withDefault("dirichlet");
  alloptions[std::string("T") + name]["bndry_all"] = alloptions[std::string("T") + name]["bndry_all"].withDefault("neumann");
  alloptions[std::string("P") + name]["bndry_all"] = alloptions[std::string("P") + name]["bndry_all"].withDefault("neumann");
  alloptions[std::string("N") + name]["bndry_all"] = alloptions[std::string("N") + name]["bndry_all"].withDefault("neumann");

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

  // Product of Dnn and another parameter has same BC as Dnn - see eqns to see why this is necessary
  DnnNn.setBoundary(std::string("Dnn") + name);
  DnnPn.setBoundary(std::string("Dnn") + name);
  DnnTn.setBoundary(std::string("Dnn") + name);
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

  Pnlim = floor(Nnlim * Tn, 1e-8);
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
  if (localstate.isSet("collision_frequency")) {
    // Dnn = Vth^2 / sigma

    if (maximum_mfp > 0) {   // MFP limit enabled
        Field3D Rnn = sqrt(Tn / AA) / (maximum_mfp / get<BoutReal>(state["units"]["meters"]));
        Dnn = (Tn / AA) / (get<Field3D>(localstate["collision_frequency"]) + Rnn);
      } else {   // MFP limit disabled
        Dnn = (Tn / AA) / (get<Field3D>(localstate["collision_frequency"]));
      }
      
  } else {  // If no collisions, hardcode max MFP to 0.1m
    output_warn.write("No collisions set for the neutrals, limiting mean free path to 0.1m");
    Field3D Rnn = sqrt(Tn / AA) / (0.1 / get<BoutReal>(state["units"]["meters"]));
    Dnn = (Tn / AA) / Rnn;
  }

  mesh->communicate(Dnn);
  Dnn.clearParallelSlices();
  Dnn.applyBoundary();

  // Neutral diffusion parameters have the same boundary condition as Dnn
  DnnPn = Dnn * Pn;
  DnnPn.applyBoundary();
  DnnNn = Dnn * Nn;
  DnnNn.applyBoundary();
  Field3D DnnNVn = Dnn * NVn;
  DnnNVn.applyBoundary();

  if (sheath_ydown) {
    for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        Dnn(r.ind, mesh->ystart - 1, jz) = -Dnn(r.ind, mesh->ystart, jz);
        DnnPn(r.ind, mesh->ystart - 1, jz) = -DnnPn(r.ind, mesh->ystart, jz);
        DnnNn(r.ind, mesh->ystart - 1, jz) = -DnnNn(r.ind, mesh->ystart, jz);
        DnnNVn(r.ind, mesh->ystart - 1, jz) = -DnnNVn(r.ind, mesh->ystart, jz);
      }
    }
  }

  if (sheath_yup) {
    for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        Dnn(r.ind, mesh->yend + 1, jz) = -Dnn(r.ind, mesh->yend, jz);
        DnnPn(r.ind, mesh->yend + 1, jz) = -DnnPn(r.ind, mesh->yend, jz);
        DnnNn(r.ind, mesh->yend + 1, jz) = -DnnNn(r.ind, mesh->yend, jz);
        DnnNVn(r.ind, mesh->yend + 1, jz) = -DnnNVn(r.ind, mesh->yend, jz);
      }
    }
  }

  // Heat conductivity
  // Note: This is kappa_n = (5/2) * Pn / (m * nu)
  //       where nu is the collision frequency used in Dnn
  kappa_n = (5. / 2) * DnnNn;

  // Viscosity
  eta_n = AA * (2. / 5) * kappa_n;

  // Sound speed appearing in Lax flux for advection terms
  Field3D sound_speed = sqrt(Tn * (5. / 3) / AA);

  // Set factors that multiply the fluxes
  particle_flux_factor = 1.0;
  momentum_flux_factor = 1.0;
  heat_flux_factor = 1.0;
  if (flux_limit) {
    // Apply flux limiters
    // Note: Fluxes calculated here are cell centre, rather than cell edge

    // Cross-field velocity
    Vector3D v_perp = -Dnn * Grad_perp(logPnlim);

    // Total velocity, perpendicular + parallel
    Vector3D v_total = v_perp;
    ASSERT2(v_total.covariant == true);
    auto* coord = mesh->getCoordinates();
    v_total.y = Vn * (coord->J * coord->Bxy); // Set parallel component

    Field3D v_sq = v_perp * v_perp + SQ(Vn); // v dot v
    Field3D v_abs = sqrt(v_sq); // |v|

    // Magnitude of the particle flux
    Field3D particle_flux_abs = Nnlim * v_abs;

    // Normalised particle flux limit
    Field3D particle_limit = Nnlim * 0.25 * sqrt(8 * Tnlim / (PI * AA));

    // Particle flux reduction factor
    if (particle_flux_limiter){
      particle_flux_factor = pow(1. + pow(particle_flux_abs / (flux_limit_alpha * particle_limit),
                                          flux_limit_gamma),
                                -1./flux_limit_gamma);
    } else {
      particle_flux_factor = 1.0;
    }

    // Flux of parallel momentum
    // Question: Should flux-limited particle flux be used here, or original flux?
    // Note: momentum flux here doesn't include pressure gradient term
    //
    // The resulting momentum_flux_factor is applied to the momentum advection and viscosity terms,
    // and so also scales the advection of kinetic energy
    Vector3D momentum_flux = NVn * v_total;
    if (neutral_viscosity) {
      momentum_flux -= eta_n * Grad_perp(Vn);
    }
    Field3D momentum_flux_abs = sqrt(momentum_flux * momentum_flux);
    Field3D momentum_limit = Pnlim;

    if (momentum_flux_limiter) {
      momentum_flux_factor = pow(1. + pow(momentum_flux_abs / (flux_limit_alpha * momentum_limit),
                                          flux_limit_gamma),
                                -1./flux_limit_gamma);
    } else {
      momentum_flux_factor = 1.0;
    }
    // Flux of heat
    // Note:
    //  - Limiting the heat flux, not energy flux
    //  - Advection and heat conduction are both vectors, and e.g
    //    could be in opposite directions.
    //  - Flux doesn't include compression term, or kinetic energy transport
    //    that is in the momentum equation
    Vector3D heat_flux = (3./2) * Pn * v_total + Pn * v_perp - kappa_n * Grad(Tn);

    Field3D heat_flux_abs = sqrt(heat_flux * heat_flux);
    Field3D heat_limit = Pnlim * sqrt(2. * Tnlim / (PI * AA));

    if (heat_flux_limiter) {
      heat_flux_factor = pow(1. + pow(heat_flux_abs / (flux_limit_alpha * heat_limit),
                                      flux_limit_gamma),
                              -1./flux_limit_gamma);
    } else {
      heat_flux_factor = 1.0;
    }
    // Communicate guard cells and apply boundary conditions
    // because the flux factors will be differentiated
    mesh->communicate(particle_flux_factor, momentum_flux_factor, heat_flux_factor);
    particle_flux_factor.applyBoundary("neumann");
    momentum_flux_factor.applyBoundary("neumann");
    heat_flux_factor.applyBoundary("neumann");
  }

  /////////////////////////////////////////////////////
  // Neutral density
  TRACE("Neutral density");

  // Note: Parallel and perpendicular flux scaled by limiter
  ddt(Nn) = -FV::Div_par_mod<hermes::Limiter>(Nn, Vn * particle_flux_factor, sound_speed) // Advection
    + FV::Div_a_Grad_perp(DnnNn * particle_flux_factor, logPnlim) // Perpendicular diffusion
    ;

  Sn = density_source; // Save for possible output
  if (localstate.isSet("density_source")) {
    Sn += get<Field3D>(localstate["density_source"]);
  }
  ddt(Nn) += Sn; // Always add density_source

  /////////////////////////////////////////////////////
  // Neutral momentum
  TRACE("Neutral momentum");

  ddt(NVn) =
      -AA * FV::Div_par_fvv<hermes::Limiter>(Nnlim, Vn * momentum_flux_factor, sound_speed) // Momentum flow
      - Grad_par(Pn)                                                 // Pressure gradient. Not included in flux limit.
      + FV::Div_a_Grad_perp(DnnNVn * momentum_flux_factor, logPnlim) // Perpendicular diffusion
      ;

  if (localstate.isSet("momentum_source")) {
    Snv = get<Field3D>(localstate["momentum_source"]);
    ddt(NVn) += Snv;
  }

  /////////////////////////////////////////////////////
  // Neutral pressure
  TRACE("Neutral pressure");

  ddt(Pn) = -FV::Div_par_mod<hermes::Limiter>(Pn, Vn * heat_flux_factor, sound_speed) // Advection
    - (2. / 3) * Pn * Div_par(Vn)                       // Compression
    + FV::Div_a_Grad_perp((5. / 3) * DnnPn * heat_flux_factor, logPnlim)   // Perpendicular advection: q = 5/2 p u_perp
    + (2. / 3) * (FV::Div_a_Grad_perp(kappa_n * heat_flux_factor, Tn)      // Perpendicular Conduction
                + FV::Div_par_K_Grad_par(kappa_n * heat_flux_factor, Tn))  // Parallel conduction
      ;

  Sp = pressure_source;
  if (localstate.isSet("energy_source")) {
    Sp += (2. / 3) * get<Field3D>(localstate["energy_source"]);
  }
  ddt(Pn) += Sp;

  if (neutral_viscosity) {
    Field3D momentum_source = FV::Div_a_Grad_perp(eta_n * momentum_flux_factor, Vn)    // Perpendicular viscosity
              + FV::Div_par_K_Grad_par(eta_n * momentum_flux_factor, Vn) // Parallel viscosity
      ;

    ddt(NVn) += momentum_source; // Viscosity
    ddt(Pn) -= (2. / 3) * Vn * momentum_source; // Viscous heating
  }

  BOUT_FOR(i, Pn.getRegion("RGN_ALL")) {
    if ((Pn[i] < 1e-9) && (ddt(Pn)[i] < 0.0)) {
      ddt(Pn)[i] = 0.0;
    }
    if ((Nn[i] < 1e-7) && (ddt(Nn)[i] < 0.0)) {
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
    set_with_attrs(state[std::string("eta_") + name], eta_n,
                   {{"time_dimension", "t"},
                    {"units", "Pa s"},
                    {"conversion", SQ(rho_s0) * Omega_ci * SI::Mp * Nnorm},
                    {"standard_name", "viscosity"},
                    {"long_name", name + " viscosity"},
                    {"source", "neutral_mixed"}});
    set_with_attrs(state[std::string("kappa_") + name], kappa_n,
                   {{"time_dimension", "t"},
                    {"units", "W / m / eV"},
                    {"conversion", SI::qe * Nnorm * rho_s0 * Cs0},
                    {"standard_name", "heat conductivity"},
                    {"long_name", name + " heat conductivity"},
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

    set_with_attrs(state[std::string("particle_flux_factor_") + name], particle_flux_factor,
                   {{"time_dimension", "t"},
                    {"units", ""},
                    {"conversion", 1.0},
                    {"standard_name", "flux factor"},
                    {"long_name", name + " particle flux factor"},
                    {"species", name},
                    {"source", "neutral_mixed"}});
    set_with_attrs(state[std::string("momentum_flux_factor_") + name], momentum_flux_factor,
                   {{"time_dimension", "t"},
                    {"units", ""},
                    {"conversion", 1.0},
                    {"standard_name", "flux factor"},
                    {"long_name", name + " momentum flux factor"},
                    {"species", name},
                    {"source", "neutral_mixed"}});
    set_with_attrs(state[std::string("heat_flux_factor_") + name], heat_flux_factor,
                   {{"time_dimension", "t"},
                    {"units", ""},
                    {"conversion", 1.0},
                    {"standard_name", "flux factor"},
                    {"long_name", name + " heat flux factor"},
                    {"species", name},
                    {"source", "neutral_mixed"}});
  }
}

void NeutralMixed::precon(const Options& state, BoutReal gamma) {
  if (!precondition) {
    return;
  }

  // Neutral gas diffusion
  // Solve (1 - gamma*Dnn*Delp2)^{-1}

  Field3D coef = -gamma * Dnn * particle_flux_factor;

  if (state.isSet("scale_timederivs")) {
    coef *= get<Field3D>(state["scale_timederivs"]);
  }

  inv->setCoefD(coef);

  ddt(Nn) = inv->solve(ddt(Nn));
  ddt(NVn) = inv->solve(ddt(NVn));
  ddt(Pn) = inv->solve(ddt(Pn));
}
