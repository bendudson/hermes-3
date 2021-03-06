
#include <bout/fv_ops.hxx>
#include <derivs.hxx>
#include <difops.hxx>
#include <bout/constants.hxx>

#include "../include/evolve_density.hxx"
#include "../include/div_ops.hxx"

using bout::globals::mesh;

EvolveDensity::EvolveDensity(std::string name, Options &alloptions, Solver *solver) : name(name) {
  AUTO_TRACE();
  
  // Evolve the density in time
  solver->add(N, std::string("N") + name);

  auto& options = alloptions[name];

  bndry_flux = options["bndry_flux"]
                      .doc("Allow flows through radial boundaries")
                      .withDefault<bool>(true);

  poloidal_flows = options["poloidal_flows"]
                       .doc("Include poloidal ExB flow")
                       .withDefault<bool>(true);

  low_n_diffuse = options["low_n_diffuse"]
                      .doc("Parallel diffusion at low density")
                      .withDefault<bool>(false);

  low_n_diffuse_perp = options["low_n_diffuse_perp"]
                           .doc("Perpendicular diffusion at low density")
                           .withDefault<bool>(false);
  
  hyper_z = options["hyper_z"].doc("Hyper-diffusion in Z").withDefault<bool>(false);

  // Charge and mass, default to electron
  charge = options["charge"].doc("Particle charge. electrons = -1").withDefault(-1.0);
  AA = options["AA"].doc("Particle atomic mass. Proton = 1").withDefault(SI::Me / SI::Mp);

  if (options["diagnose"]
          .doc("Output additional diagnostics?")
          .withDefault<bool>(false)) {
    bout::globals::dump.addRepeat(ddt(N), std::string("ddt(N") + name + std::string(")"));
    bout::globals::dump.addRepeat(Sn, std::string("SN") + name);
    Sn = 0.0;
  }

  const Options& units = alloptions["units"];
  const BoutReal Nnorm = units["inv_meters_cubed"];
  const BoutReal Omega_ci = 1. / units["seconds"].as<BoutReal>();

  source = alloptions[std::string("N") + name]["source"]
               .doc("Source term in ddt(N" + name + std::string("). Units [m^-3/s]"))
               .withDefault(Field3D(0.0))
           / (Nnorm * Omega_ci);
}

void EvolveDensity::transform(Options &state) {
  AUTO_TRACE();
  mesh->communicate(N);
  
  auto& species = state["species"][name];
  set(species["density"], N);
  set(species["AA"], AA); // Atomic mass
  if (charge != 0.0) { // Don't set charge for neutral species
    set(species["charge"], charge);
  }
}

void EvolveDensity::finally(const Options &state) {
  AUTO_TRACE();

  // Get the coordinate system
  auto coord = N.getCoordinates();
  
  auto& species = state["species"][name];

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
    
    // Typical wave speed used for numerical diffusion
    Field3D sound_speed;
    if (state.isSet("sound_speed")) {
      Field3D sound_speed = get<Field3D>(state["sound_speed"]);
    } else {
      Field3D T = get<Field3D>(species["temperature"]);
      sound_speed = sqrt(T);
    }

    if (state.isSection("fields") and state["fields"].isSet("phi")) {
      // Parallel wave speed increased to electron sound speed
      // since electrostatic & electromagnetic waves are supported
      ddt(N) -= FV::Div_par(N, V, sqrt(SI::Me / SI::Mp) * sound_speed);
    } else {
      // Parallel wave speed is ion sound speed
      ddt(N) -= FV::Div_par(N, V, sound_speed);
    }
  }

  if (low_n_diffuse) {
    // Diffusion which kicks in at very low density, in order to
    // help prevent negative density regions
    ddt(N) += FV::Div_par_K_Grad_par(SQ(coord->dy) * coord->g_22 * 1e-4 / N, N);
  }
  if (low_n_diffuse_perp) {
    ddt(N) += Div_Perp_Lap_FV_Index(1e-4 / N, N, bndry_flux);
  }

  if (hyper_z > 0.) {
    ddt(N) -= hyper_z * SQ(SQ(coord->dz)) * D4DZ4(N);
  }

  Sn = source; // Save for possible output
  if (species.isSet("density_source")) {
    Sn += get<Field3D>(species["density_source"]);
  }
  ddt(N) += Sn;

#if CHECK > 1
  bout::checkFinite(ddt(N), std::string("ddt N") + name, "RGN_NOBNDRY");
#endif
}
