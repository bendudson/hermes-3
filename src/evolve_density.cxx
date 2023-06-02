
#include <bout/constants.hxx>
#include <bout/fv_ops.hxx>
#include <bout/output_bout_types.hxx>
#include <bout/derivs.hxx>
#include <bout/difops.hxx>
#include <bout/initialprofiles.hxx>

#include "../include/div_ops.hxx"
#include "../include/evolve_density.hxx"
#include "../include/hermes_utils.hxx"
#include "../include/hermes_build_config.hxx"

using bout::globals::mesh;

EvolveDensity::EvolveDensity(std::string name, Options& alloptions, Solver* solver)
    : name(name) {
  AUTO_TRACE();

  auto& options = alloptions[name];

  bndry_flux = options["bndry_flux"]
                   .doc("Allow flows through radial boundaries")
                   .withDefault<bool>(true);

  poloidal_flows =
      options["poloidal_flows"].doc("Include poloidal ExB flow").withDefault<bool>(true);

  density_floor = options["density_floor"].doc("Minimum density floor").withDefault(1e-5);

  low_n_diffuse = options["low_n_diffuse"]
                      .doc("Parallel diffusion at low density")
                      .withDefault<bool>(false);

  low_n_diffuse_perp = options["low_n_diffuse_perp"]
                           .doc("Perpendicular diffusion at low density")
                           .withDefault<bool>(false);

  hyper_z = options["hyper_z"].doc("Hyper-diffusion in Z").withDefault(-1.0);

  evolve_log = options["evolve_log"]
                   .doc("Evolve the logarithm of density?")
                   .withDefault<bool>(false);

  if (evolve_log) {
    // Evolve logarithm of density
    solver->add(logN, std::string("logN") + name);
    // Save the density to the restart file
    // so the simulation can be restarted evolving density
    // get_restart_datafile()->addOnce(N, std::string("N") + name);

    if (!alloptions["hermes"]["restarting"]) {
      // Set logN from N input options
      initial_profile(std::string("N") + name, N);
      logN = log(N);
    } else {
      // Ignore these settings
      Options::root()[std::string("N") + name].setConditionallyUsed();
    }
  } else {
    // Evolve the density in time
    solver->add(N, std::string("N") + name);
  }

  // Charge and mass
  charge = options["charge"].doc("Particle charge. electrons = -1");
  AA = options["AA"].doc("Particle atomic mass. Proton = 1");

  diagnose =
      options["diagnose"].doc("Output additional diagnostics?").withDefault<bool>(false);

  const Options& units = alloptions["units"];
  const BoutReal Nnorm = units["inv_meters_cubed"];
  const BoutReal Omega_ci = 1. / units["seconds"].as<BoutReal>();

  // Try to read the density source from the mesh
  // Units of particles per cubic meter per second
  source = 0.0;
  mesh->get(source, std::string("N") + name + "_src");
  // Allow the user to override the source
  source = alloptions[std::string("N") + name]["source"]
               .doc("Source term in ddt(N" + name + std::string("). Units [m^-3/s]"))
               .withDefault(source)
           / (Nnorm * Omega_ci);

  if (alloptions[std::string("N") + name]["source_only_in_core"]
      .doc("Zero the source outside the closed field-line region?")
      .withDefault<bool>(false)) {
    for (int x = mesh->xstart; x <= mesh->xend; x++) {
      if (!mesh->periodicY(x)) {
        // Not periodic, so not in core
        for (int y = mesh->ystart; y <= mesh->yend; y++) {
          for (int z = mesh->zstart; z <= mesh->zend; z++) {
            source(x, y, z) = 0.0;
          }
        }
      }
    }
  }

  neumann_boundary_average_z = alloptions[std::string("N") + name]["neumann_boundary_average_z"]
    .doc("Apply neumann boundary with Z average?")
    .withDefault<bool>(false);
}

void EvolveDensity::transform(Options& state) {
  AUTO_TRACE();

  if (evolve_log) {
    // Evolving logN, but most calculations use N
    N = exp(logN);
  }

  mesh->communicate(N);

  if (neumann_boundary_average_z) {
    // Take Z (usually toroidal) average and apply as X (radial) boundary condition
    if (mesh->firstX()) {
      for (int j = mesh->ystart; j <= mesh->yend; j++) {
        BoutReal Navg = 0.0; // Average N in Z
        for (int k = 0; k < mesh->LocalNz; k++) {
          Navg += N(mesh->xstart, j, k);
        }
        Navg /= mesh->LocalNz;

        // Apply boundary condition
        for (int k = 0; k < mesh->LocalNz; k++) {
          N(mesh->xstart - 1, j, k) = 2. * Navg - N(mesh->xstart, j, k);
          N(mesh->xstart - 2, j, k) = N(mesh->xstart - 1, j, k);
        }
      }
    }

    if (mesh->lastX()) {
      for (int j = mesh->ystart; j <= mesh->yend; j++) {
        BoutReal Navg = 0.0; // Average N in Z
        for (int k = 0; k < mesh->LocalNz; k++) {
          Navg += N(mesh->xend, j, k);
        }
        Navg /= mesh->LocalNz;

        for (int k = 0; k < mesh->LocalNz; k++) {
          N(mesh->xend + 1, j, k) = 2. * Navg - N(mesh->xend, j, k);
          N(mesh->xend + 2, j, k) = N(mesh->xend + 1, j, k);
        }
      }
    }
  }

  auto& species = state["species"][name];
  set(species["density"], floor(N, 0.0)); // Density in state always >= 0
  set(species["AA"], AA);                 // Atomic mass
  if (charge != 0.0) {                    // Don't set charge for neutral species
    set(species["charge"], charge);
  }

  if (low_n_diffuse) {
    // Calculate a diffusion coefficient which can be used in N, P and NV equations

    auto* coord = mesh->getCoordinates();

    Field3D low_n_coeff =
        SQ(coord->dy) * coord->g_22
        * log(density_floor / clamp(N, 1e-3 * density_floor, density_floor));
    low_n_coeff.applyBoundary("neumann");
    set(species["low_n_coeff"], low_n_coeff);
  }
}

void EvolveDensity::finally(const Options& state) {
  AUTO_TRACE();

  auto& species = state["species"][name];

  // Get density boundary conditions
  // but retain densities which fall below zero
  N.setBoundaryTo(get<Field3D>(species["density"]));

  if (state.isSection("fields") and state["fields"].isSet("phi")) {
    // Electrostatic potential set -> include ExB flow

    Field3D phi = get<Field3D>(state["fields"]["phi"]);

    ddt(N) = -Div_n_bxGrad_f_B_XPPM(N, phi, bndry_flux, poloidal_flows,
                                    true); // ExB drift
  } else {
    ddt(N) = 0.0;
  }

  if (species.isSet("velocity")) {
    // Parallel velocity set
    Field3D V = get<Field3D>(species["velocity"]);

    // Wave speed used for numerical diffusion
    // Note: For simulations where ion density is evolved rather than electron density,
    // the fast electron dynamics still determine the stability.
    Field3D fastest_wave;
    if (state.isSet("fastest_wave")) {
      fastest_wave = get<Field3D>(state["fastest_wave"]);
    } else {
      Field3D T = get<Field3D>(species["temperature"]);
      BoutReal AA = get<BoutReal>(species["AA"]);
      fastest_wave = sqrt(T / AA);
    }

    ddt(N) -= FV::Div_par_mod<hermes::Limiter>(N, V, fastest_wave);
  }

  if (low_n_diffuse) {
    // Diffusion which kicks in at very low density, in order to
    // help prevent negative density regions

    Field3D low_n_coeff = get<Field3D>(species["low_n_coeff"]);
    ddt(N) += FV::Div_par_K_Grad_par(low_n_coeff, N);
  }
  if (low_n_diffuse_perp) {
    ddt(N) += Div_Perp_Lap_FV_Index(density_floor / floor(N, 1e-3 * density_floor), N,
                                    bndry_flux);
  }

  if (hyper_z > 0.) {
    auto* coord = N.getCoordinates();
    ddt(N) -= hyper_z * SQ(SQ(coord->dz)) * D4DZ4(N);
  }

  Sn = source; // Save for possible output
  if (species.isSet("density_source")) {
    Sn += get<Field3D>(species["density_source"]);
  }
  ddt(N) += Sn;

  // Scale time derivatives
  if (state.isSet("scale_timederivs")) {
    ddt(N) *= get<Field3D>(state["scale_timederivs"]);
  }

  if (evolve_log) {
    ddt(logN) = ddt(N) / N;
  }

#if CHECKLEVEL >= 1
  for (auto& i : N.getRegion("RGN_NOBNDRY")) {
    if (!std::isfinite(ddt(N)[i])) {
      throw BoutException("ddt(N{}) non-finite at {}. Sn={}\n", name, i, Sn[i]);
    }
  }
#endif
}

void EvolveDensity::outputVars(Options& state) {
  // Normalisations
  auto Nnorm = get<BoutReal>(state["Nnorm"]);
  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);

  if (evolve_log) {
    // Save density to output files
    state[std::string("N") + name].force(N);
  }
  state[std::string("N") + name].setAttributes({{"time_dimension", "t"},
                                                {"units", "m^-3"},
                                                {"conversion", Nnorm},
                                                {"standard_name", "density"},
                                                {"long_name", name + " number density"},
                                                {"species", name},
                                                {"source", "evolve_density"}});

  if (diagnose) {
    set_with_attrs(
        state[std::string("ddt(N") + name + std::string(")")], ddt(N),
        {{"time_dimension", "t"},
         {"units", "m^-3 s^-1"},
         {"conversion", Nnorm * Omega_ci},
         {"long_name", std::string("Rate of change of ") + name + " number density"},
         {"species", name},
         {"source", "evolve_density"}});

    set_with_attrs(state[std::string("SN") + name], Sn,
                   {{"time_dimension", "t"},
                    {"units", "m^-3 s^-1"},
                    {"conversion", Nnorm * Omega_ci},
                    {"standard_name", "density source"},
                    {"long_name", name + " number density source"},
                    {"species", name},
                    {"source", "evolve_density"}});

    set_with_attrs(state[std::string("S") + name + std::string("_src")], source,
                   {{"time_dimension", "t"},
                    {"units", "m^-3 s^-1"},
                    {"conversion", Nnorm * Omega_ci},
                    {"standard_name", "density source"},
                    {"long_name", name + " number density source"},
                    {"species", name},
                    {"source", "evolve_density"}});
  }
}
