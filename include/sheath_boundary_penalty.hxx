#pragma once
#ifndef SHEATH_BOUNDARY_PENALTY_H
#define SHEATH_BOUNDARY_PENALTY_H

#include "component.hxx"
#include <bout/region.hxx>

/// Penalty method for imposing sheath boundary conditions over an
/// arbitrary shaped wall that is immersed inside the structured mesh.
///
/// Notes:
///  - This should be applied AFTER sheath boundaries have been applied
///    An error should be raised if this is done in the wrong order.
///
struct SheathBoundaryPenalty : public Component {
  /// # Input options
  /// - <name> e.g. "sheath_boundary_penalty"
  ///   - penalty_timescale    Timescale in seconds
  ///
  /// # Reads from the mesh
  /// - penalty_mask[x,y,z]    A 3D field defining the shape of the boundary
  ///                          Equal 0 inside the plasma, 1 in the wall
  SheathBoundaryPenalty(std::string name, Options& options, Solver*);

  /// # Inputs
  ///
  /// # Outputs
  ///
  ///
  void transform(Options& state) override;

  void outputVars(Options& state) override;

private:
  /// Mask function that is 0 in the plasma, 1 in the wall
  Field3D penalty_mask;

  /// Cell indices where penalty_mask > 0
  Region<Ind3D> penalty_region;

  BoutReal gamma_e, gamma_i; // Sheath heat transmission factors

  /// Timescale of penalisation [normalised]
  BoutReal penalty_timescale;

  /// Diagnostics for output
  Options diagnostics;
  bool diagnose; ///< Save diagnostics?
};

namespace {
RegisterComponent<SheathBoundaryPenalty>
    registercomponentsheathboundarypenalty("sheath_boundary_penalty");
}

#endif // SHEATH_BOUNDARY_PENALTY_H
