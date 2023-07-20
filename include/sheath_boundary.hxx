#pragma once
#ifndef SHEATH_BOUNDARY_H
#define SHEATH_BOUNDARY_H

#include "component.hxx"

/// Boundary condition at the wall in Y
///
/// This is a collective component, because it couples all charged species
///
/// These are based on
/// "Boundary conditions for the multi-ion magnetized plasma-wall transition"
///  by D.Tskhakaya, S.Kuhn. JNM 337-339 (2005), 405-409
///
/// Notes:
///   - The approximation used here is for ions having similar
///     gyro-orbit sizes
///   - No boundary condition is applied to neutral species
///   - Boundary conditions are applied to field-aligned fields
///     using to/fromFieldAligned
///
struct SheathBoundary : public Component {
  /// # Input options
  /// - <name>  e.g. "sheath_boundary"
  ///   - lower_y                  Boundary on lower y?
  ///   - upper_y                  Boundary on upper y?
  ///   - wall_potential           Voltage of the wall [Volts]
  ///   - floor_potential          Apply floor to sheath potential?
  ///   - secondary_electron_coef  Effective secondary electron emission coefficient
  ///   - sin_alpha                Sine of the angle between magnetic field line and wall surface (0 to 1)
  ///   - always_set_phi           Always set phi field? Default is to only modify if already set
  SheathBoundary(std::string name, Options &options, Solver *);

  ///
  /// # Inputs
  /// - species
  ///   - e
  ///     - density
  ///     - temperature
  ///     - pressure    Optional
  ///     - velocity    Optional
  ///     - mass        Optional
  ///     - adiabatic   Optional. Ratio of specific heats, default 5/3.
  ///   - <ions>  if charge is set (i.e. not neutrals)
  ///     - charge
  ///     - mass
  ///     - density
  ///     - temperature
  ///     - pressure     Optional
  ///     - velocity     Optional. Default 0
  ///     - momentum     Optional. Default mass * density * velocity
  ///     - adiabatic    Optional. Ratio of specific heats, default 5/3.
  /// - fields
  ///   - phi    Optional. If not set, calculated at boundary (see note below)
  ///
  /// # Outputs
  /// - species
  ///   - e
  ///     - density      Sets boundary
  ///     - temperature  Sets boundary
  ///     - velocity     Sets boundary
  ///     - energy_source
  ///   - <ions>
  ///     - density      Sets boundary
  ///     - temperature  Sets boundary
  ///     - velocity     Sets boundary
  ///     - momentum     Sets boundary
  ///     - energy_source
  /// - fields
  ///   - phi   Sets boundary
  ///
  /// If the field phi is set, then this is used in the boundary condition.
  /// If not set, phi at the boundary is calculated and stored in the state.
  /// Note that phi in the domain will not be set, so will be invalid data.
  ///
  ///
  void transform(Options &state) override;
private:
  BoutReal Ge; // Secondary electron emission coefficient
  BoutReal sin_alpha; // sin of angle between magnetic field and wall.
  
  bool lower_y; // Boundary on lower y?
  bool upper_y; // Boundary on upper y?

  bool always_set_phi; ///< Set phi field?

  Field3D wall_potential; ///< Voltage at the wall. Normalised units.

  bool floor_potential; ///< Apply floor to sheath potential?
};

namespace {
RegisterComponent<SheathBoundary>
    registercomponentsheathboundary("sheath_boundary");
}

#endif // SHEATH_BOUNDARY_H
