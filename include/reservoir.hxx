#pragma once
#ifndef RESERVOIR_H
#define RESERVOIR_H

#include "component.hxx"

/// Simple 1D reservoir model
///
/// This is a simple reservoir model. It's meant for 1D but will work in 2D and 3D as
/// well. Define the reservoir location by using an analytical expression in the input
/// file e.g. reservoir_location = H(y - y_xpt) where y_xpt is the x-point location would
/// create one between the target and the X-point. The reservoir is set to a constant
/// density and your chosen species will equilibrate with it based on the density
/// difference and a provided timescale. The intended use is to provide some form of
/// cross-field transport for neutrals in 1D.
struct Reservoir : public Component {

  /// # Inputs
  /// - <name>
  ///   - reservoir_density     BoutReal, Constant density in m^-3
  ///   - reservoir_location    Field3D, >0 indicates reservoir location
  ///   - reservoir_timescale   BoutReal, timescale in seconds
  ///   - diagnose              bool. Save sources/sinks
  Reservoir(std::string name, Options& alloptions, Solver*);

  /// # Inputs
  /// - species
  ///   - <name>
  ///     - density
  ///     - pressure
  ///     - momentum
  ///
  /// # Outputs
  /// - species
  ///   - <name>
  ///     - density_source
  ///     - energy_source
  ///     - momentum_source
  void transform(Options& state) override;
  void outputVars(Options& state) override;

private:
  std::string name;           ///< Short name of the species e.g. h+
  Field3D reservoir_location; ///< Indicates reservoir if >0
  BoutReal reservoir_density, reservoir_timescale;
  bool diagnose;

  /// Cell indices where reservoir_location > 0
  Region<Ind3D> reservoir_region;

  Field3D density_source, energy_source, momentum_source;
};

namespace {
RegisterComponent<Reservoir> registercomponentreservoir("reservoir");
}

#endif // RESERVOIR_H
