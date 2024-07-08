
#pragma once
#ifndef NEUTRAL_MIXED_H
#define NEUTRAL_MIXED_H

#include <memory>
#include <string>

#include <bout/invert_laplace.hxx>

#include "component.hxx"

/// Evolve density, parallel momentum and pressure
/// for a neutral gas species with cross-field diffusion
struct NeutralMixed : public Component {
  ///
  /// @param name     The name of the species e.g. "h"
  /// @param options  Top-level options. Settings will be taken from options[name]
  /// @param solver   Time-integration solver to be used
  NeutralMixed(const std::string& name, Options& options, Solver *solver);
  
  /// Modify the given simulation state
  void transform(Options &state) override;
  
  /// Use the final simulation state to update internal state
  /// (e.g. time derivatives)
  void finally(const Options &state) override;

  /// Add extra fields for output, or set attributes e.g docstrings
  void outputVars(Options &state) override;

  /// Preconditioner
  void precon(const Options &state, BoutReal gamma) override;
private:
  std::string name;  ///< Species name
  
  Field3D Nn, Pn, NVn; // Density, pressure and parallel momentum
  Field3D Vn; ///< Neutral parallel velocity
  Field3D Tn; ///< Neutral temperature
  Field3D Nnlim, Pnlim, logPnlim, Vnlim, Tnlim; // Limited in regions of low density

  BoutReal AA; ///< Atomic mass (proton = 1)

  Field3D Dnn, Dnn_unlimited; ///< Diffusion coefficient
  Field3D DnnNn, DnnPn, DnnNVn;
  Field3D Dmax;
  Field3D gradlogP, gradperplogP, gradperpT;
  Field3D nu; ///< Collisionality to use for diffusion
  std::vector<std::string> collision_names; ///< Collisions used for collisionality
  std::string diffusion_collisions_mode;  ///< Collision selection, either afn or legacy

  bool sheath_ydown, sheath_yup;

  BoutReal nn_floor; ///< Minimum Nn used when dividing NVn by Nn to get Vn.
  BoutReal pn_floor; ///< Minimum Pn used when dividing Pn by Nn to get Tn.

  BoutReal flux_limit; ///< Diffusive flux limit
  BoutReal diffusion_limit;    ///< Maximum diffusion coefficient
  BoutReal maximum_mfp;   ///< Reduce diffusion using physical MFP limit
  Field3D vth;   ///< Thermal speed to use in flux limiter
  bool legacy_vth_limiter;
  bool asymptotic_limiter, asymptotic_limiter_advection;
  bool constant_conduction;
  bool override_limiter;   // Force conduction and viscosity limiters to use the advection limiter?
  bool legacy_limiter; 
  bool legacy_separate_conduction;
  bool first_RHS {true};
  
  // Limit perpendicular advection fluxes to fraction of thermal speed.
  BoutReal advection_limit_alpha, conduction_limit_alpha, viscosity_limit_alpha;  // Limiters
  BoutReal flux_limit_gamma;   ///< Sharpness of perpendicular flux limitation. Usually 1 or 2
  Field3D advection_factor, conduction_factor, viscosity_factor;   // Multiplier to achieve limit
  Field3D advection_flux_abs, advection_limit;

  bool neutral_viscosity; ///< include viscosity?
  bool neutral_conduction; ///< Include heat conduction?
  bool evolve_momentum; ///< Evolve parallel momentum?
  
  Field3D kappa_n, eta_n; ///< Neutral conduction and viscosity
  Field3D kappa_n_unlimited, kappa_n_Dnchained, kappa_n_max;

  bool precondition {true}; ///< Enable preconditioner?
  bool lax_flux; ///< Use Lax flux for advection terms
  
  std::unique_ptr<Laplacian> inv; ///< Laplacian inversion used for preconditioning

  Field3D density_source, pressure_source; ///< External input source
  Field3D Sn, Sp, Snv; ///< Particle, pressure and momentum source
  Field3D sound_speed; ///< Sound speed for use with Lax flux
  Field3D perp_nn_adv_src; ///< Source due to perpendicular advection operator
  Field3D par_nn_adv_src; ///< Source due to parallel advection operator

  bool output_ddt; ///< Save time derivatives?
  bool diagnose; ///< Save additional diagnostics?

  // Flow diagnostics
  Field3D particle_flow_xlow, particle_flow_ylow;
  Field3D momentum_flow_xlow, momentum_flow_ylow;
  Field3D energy_flow_xlow, energy_flow_ylow;
  Field3D conduction_flow_xlow, conduction_flow_ylow;
};

namespace {
RegisterComponent<NeutralMixed> registersolverneutralmixed("neutral_mixed");
}

#endif // NEUTRAL_MIXED_H
