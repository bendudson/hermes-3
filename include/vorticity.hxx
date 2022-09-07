#pragma once
#ifndef VORTICITY_H
#define VORTICITY_H

#include <vector2d.hxx>

#include "component.hxx"

class LaplaceXY;
class Laplacian;

/// Evolve electron density in time
/// 
struct Vorticity : public Component {
  /// Options
  ///
  /// - <name>
  ///   - diamagnetic
  ///   - diamagnetic_polarisation
  ///   - average_atomic_mass
  ///   - bndry_flux
  ///   - poloidal_flows
  ///   - split_n0
  ///   - laplacian
  ///     Options for the Laplacian phi solver
  /// 
  Vorticity(std::string name, Options &options, Solver *solver);

  /// Optional inputs
  ///
  /// - species
  ///   - pressure and charge => Calculates diamagnetic terms [if diamagnetic=true]
  ///   - pressure, charge and mass => Calculates polarisation current terms [if diamagnetic_polarisation=true]
  /// 
  /// Sets in the state
  /// - species
  ///   - [if has pressure and charge]
  ///     - energy_source
  /// - fields
  ///   - vorticity
  ///   - phi         Electrostatic potential
  ///   - DivJdia     Divergence of diamagnetic current [if diamagnetic=true]
  ///
  /// Note: Diamagnetic current calculated here, but could be moved
  ///       to a component with the diamagnetic drift advection terms
  void transform(Options &state) override;

  /// Optional inputs
  /// - fields
  ///   - DivJextra    Divergence of current, including parallel current
  ///                  Not including diamagnetic or polarisation currents
  /// 
  void finally(const Options &state) override;
  
private:
  Field3D Vort; // Evolving vorticity
  
  Field3D phi; // Electrostatic potential
  std::unique_ptr<Laplacian> phiSolver; // Laplacian solver in X-Z

  bool exb_advection; // Include nonlinear ExB advection?
  bool diamagnetic; // Include diamagnetic current?
  bool diamagnetic_polarisation; // Include diamagnetic drift in polarisation current
  bool boussinesq; // Use 'Boussinesq' approximation?
  BoutReal average_atomic_mass; // Weighted average atomic mass, for polarisaion current (Boussinesq approximation)
  bool poloidal_flows;   ///< Include poloidal ExB flow?
  bool bndry_flux;  ///< Allow flows through radial boundaries?

  bool collisional_friction; ///< Damping of vorticity due to collisional friction

  bool sheath_boundary; ///< Set outer boundary to j=0?

  bool vort_dissipation; ///< Parallel dissipation of vorticity
  bool phi_dissipation;  ///< Parallel dissipation of potential

  bool phi_boundary_relax; ///< Relax boundary to zero-gradient
  BoutReal phi_boundary_timescale; ///< Relaxation timescale [normalised]
  BoutReal phi_boundary_last_update; ///< Time when last updated

  bool split_n0; // Split phi into n=0 and n!=0 components
  LaplaceXY* laplacexy; // Laplacian solver in X-Y (n=0)

  Field2D Bsq; // SQ(coord->Bxy)
  Vector2D Curlb_B; // Curvature vector Curl(b/B)
  BoutReal hyper_z; ///< Hyper-viscosity in Z

  // Intermediate variables to carry between transform() and finally()
  Field3D N_sum, Pi_sum;

  // Diagnostic outputs
  Field3D DivJdia, DivJcol; // Divergence of diamagnetic and collisional current
};

namespace {
RegisterComponent<Vorticity> registercomponentvorticity("vorticity");
}

#endif // VORTICITY_H
