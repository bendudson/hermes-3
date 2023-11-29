#pragma once
#ifndef TOKAMAK_CORE_H
#define TOKAMAK_CORE_H

#include <bout/constants.hxx>
#include "component.hxx"

/// Allows to set boundary conditions on the core in a tokamak simulation
struct TokamakCore : public Component {

  /// Comments on inputs will go here
  ///  - <name> (e.g. "d+")
  ///    - input        comment

  TokamakCore(std::string name, Options& alloptions, Solver*) : name(name) {
    const auto& units = alloptions["units"];
    const BoutReal Nnorm = get<BoutReal>(units["inv_meters_cubed"]);
    const BoutReal Tnorm = get<BoutReal>(units["eV"]);
    const BoutReal Omega_ci = 1. / units["seconds"].as<BoutReal>();


    Options& options = alloptions[name];

    power = options["core_power"]
                               .doc("Inflow of energy in the core [W]. -1 = off")
                               .withDefault<BoutReal>(-1.0) / (SI::qe * Tnorm * Omega_ci);

    particle_flow = options["core_particle_flow"]
                               .doc("Inflow of particles in the core [s^-1]. -1 = off")
                               .withDefault<BoutReal>(-1.0) / (Omega_ci);

    diagnose = options["diagnose"]
                   .doc("Output additional diagnostics?")
                   .withDefault<bool>(false);
  }

  /// Inputs
  ///  - <name>
  ///    - density
  ///
  /// Outputs
  ///
  ///  - <name>
  ///    - density_source
  ///
  void transform(Options& state) override;

  void outputVars(Options& state) override {
    AUTO_TRACE();

    auto Nnorm = get<BoutReal>(state["Nnorm"]);
    auto Tnorm = get<BoutReal>(state["Tnorm"]);
    auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
    auto rho_s0 = get<BoutReal>(state["rho_s0"]);
    const BoutReal Pnorm = SI::qe * Tnorm * Nnorm; 

    if (diagnose) {

      if (power > 0) {
        set_with_attrs(state[std::string("E") + name + std::string("_core_src")], energy_source,
                    {{"time_dimension", "t"},
                      {"units", "W m^-3"},
                      {"conversion", Nnorm * Tnorm * SI::qe / (rho_s0*rho_s0*rho_s0)},
                      {"standard_name", "energy source"},
                      {"long_name", name + " core energy source"},
                      {"species", name},
                      {"source", "tokamak_core"}});
      }

      if (particle_flow > 0) {
        set_with_attrs(state[std::string("S") + name + std::string("_core_src")], density_source,
                    {{"time_dimension", "t"},
                      {"units", "m^-3 s^-1"},
                      {"conversion", Omega_ci / (rho_s0*rho_s0*rho_s0)},
                      {"standard_name", "density source"},
                      {"long_name", name + " core number density source"},
                      {"species", name},
                      {"source", "tokamak_core"}});
      }

    }
  }


private:
  std::string name; ///< The species name

  BoutReal power;            ///< Core energy inflow in [W]
  Field3D energy_source;     ///< Core energy source in [W/m^3]
  BoutReal particle_flow;    ///< Core particle inflow in [s^-1]
  Field3D density_source;   ///< Core particle source in [m^-3/s]
  BoutReal core_volume_local, core_volume;      ///< Volume of first core ring in [m^3] for one processor, and for whole domain

  Field2D loopmarker; /// Debug only

  bool diagnose; ///< Output diagnostic information?
};

namespace {
RegisterComponent<TokamakCore> register_tokamak_core("tokamak_core");
}

#endif // TOKAMAK_CORE
