/// Ion viscosity model

#include <bout/constants.hxx>
#include <bout/fv_ops.hxx>
#include <bout/difops.hxx>
#include <bout/output_bout_types.hxx>
#include <bout/mesh.hxx>

#include "../include/ion_viscosity.hxx"
#include "../include/div_ops.hxx"

using bout::globals::mesh;

IonViscosity::IonViscosity(std::string name, Options& alloptions, Solver*) {
  auto& options = alloptions[name];

  eta_limit_alpha = options["eta_limit_alpha"]
    .doc("Viscosity flux limiter coefficient. <0 = turned off")
    .withDefault(-1.0);

  diagnose = options["diagnose"]
      .doc("Output additional diagnostics?")
      .withDefault<bool>(false);

  perpendicular = options["perpendicular"]
    .doc("Include perpendicular flow? (Requires phi)")
    .withDefault<bool>(false);
  
  if (perpendicular) {
    // Read curvature vector
    try {
      Curlb_B.covariant = false; // Contravariant
      mesh->get(Curlb_B, "bxcv");
    } catch (BoutException& e) {
      // May be 2D, reading as 3D
      Vector2D curv2d;
      curv2d.covariant = false;
      mesh->get(curv2d, "bxcv");
      Curlb_B = curv2d;
    }

    if (Options::root()["mesh"]["paralleltransform"]["type"].as<std::string>()
        == "shifted") {
      Field2D I;
      mesh->get(I, "sinty");
      Curlb_B.z += I * Curlb_B.x;
    }

    // Normalisations
    const Options& units = alloptions["units"];
    const BoutReal Bnorm = units["Tesla"];
    const BoutReal Lnorm = units["meters"];

    auto coord = mesh->getCoordinates();

    Curlb_B.x /= Bnorm;
    Curlb_B.y *= SQ(Lnorm);
    Curlb_B.z *= SQ(Lnorm);
    
    Curlb_B *= 2. / coord->Bxy;
  }
}

void IonViscosity::transform(Options &state) {
  AUTO_TRACE();

  Options& allspecies = state["species"];

  auto coord = mesh->getCoordinates();
  const Field3D Bxy = coord->Bxy;
  const Field3D sqrtB = sqrt(Bxy);
  const Field3D Grad_par_logB = Grad_par(log(Bxy));

  // Loop through all species
  for (auto& kv : allspecies.getChildren()) {
    if (kv.first == "e") {
      continue; // Skip electrons -> only ions
    }
    const auto& species_name = kv.first;

    Options& species = allspecies[species_name];

    if (!(isSetFinal(species["pressure"], "ion_viscosity") and
          isSetFinal(species["velocity"], "ion_viscosity") and
          isSetFinal(species["charge"], "ion_viscosity"))) {
      // Species doesn't have a pressure, velocity and charge => Skip
      continue;
    }

    if (std::fabs(get<BoutReal>(species["charge"])) < 1e-3) {
      // No charge
      continue;
    }

    const Field3D tau = 1. / get<Field3D>(species["collision_frequency"]);
    const Field3D P = get<Field3D>(species["pressure"]);
    const Field3D V = get<Field3D>(species["velocity"]);

    // Parallel ion viscosity (4/3 * 0.96 coefficient)
    Field3D eta = 1.28 * P * tau;

    if (eta_limit_alpha > 0.) {
      // SOLPS-style flux limiter
      // Values of alpha ~ 0.5 typically

      const Field3D q_cl = eta * Grad_par(V);   // Collisional value
      const Field3D q_fl = eta_limit_alpha * P; // Flux limit

      eta = eta / (1. + abs(q_cl / q_fl));

      eta.getMesh()->communicate(eta);
      eta.applyBoundary("neumann");
    }
    
    // This term is the parallel flow part of
    // -(2/3) B^(3/2) Grad_par(Pi_ci / B^(3/2))
    const Field3D div_Pi_cipar = sqrtB * FV::Div_par_K_Grad_par(eta / Bxy, sqrtB * V);

    add(species["momentum_source"], div_Pi_cipar);
    subtract(species["energy_source"], V * div_Pi_cipar); // Internal energy

    if (!perpendicular) {
      continue; // Skip perpendicular flow parts below
    }

    // Need electrostatic potential for ExB flow
    const Field3D phi = get<Field3D>(state["fields"]["phi"]);

    // Density and temperature needed for diamagnetic flows
    const Field3D N = get<Field3D>(species["density"]);
    const Field3D T = get<Field3D>(species["temperature"]);

    // Parallel ion stress tensor component
    Field3D Pi_cipar = -0.96 * P * tau *
      (2. * Grad_par(V) + V * Grad_par_logB);
    // Could also be written as:
    // Pi_cipar = -0.96*Pi*tau*2.*Grad_par(sqrt(Bxy)*Vi)/sqrt(Bxy);

    // Perpendicular ion stress tensor
    // 0.96 P tau kappa * (V_E + V_di + 1.61 b x Grad(T)/B )
    // Note: Heat flux terms are neglected for now
    Field3D Pi_ciperp = -0.5 * 0.96 * P * tau *
      (Curlb_B * Grad(phi + 1.61 * T) - Curlb_B * Grad(P) / N);

    // Limit size of stress tensor components
    // If the off-diagonal components of the pressure tensor are large compared
    // to the scalar pressure, then the model is probably breaking down.
    // This can happen in very low density regions
    BOUT_FOR(i, Pi_ciperp.getRegion("RGN_NOBNDRY")) {
      if (fabs(Pi_ciperp[i]) > P[i]) {
        Pi_ciperp[i] = SIGN(Pi_ciperp[i]) * P[i];
      }
      if (fabs(Pi_cipar[i]) > P[i]) {
        Pi_cipar[i] = SIGN(Pi_cipar[i]) * P[i];
      }
    }

    // Apply perpendicular boundary conditions
    Pi_ciperp.applyBoundary("neumann");
    Pi_cipar.applyBoundary("neumann");

    // Apply parallel boundary conditions
    Pi_ciperp = toFieldAligned(Pi_ciperp);
    Pi_cipar = toFieldAligned(Pi_cipar);
    {
      int jy = 1;
      for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
	for (int jz = 0; jz < mesh->LocalNz; jz++) {
	  Pi_ciperp(r.ind, jy, jz) = Pi_ciperp(r.ind, jy + 1, jz);
	  Pi_cipar(r.ind, jy, jz) = Pi_cipar(r.ind, jy + 1, jz);
	}
      }
      jy = mesh->yend + 1;
      for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
	for (int jz = 0; jz < mesh->LocalNz; jz++) {
	  Pi_ciperp(r.ind, jy, jz) = Pi_ciperp(r.ind, jy - 1, jz);
	  Pi_cipar(r.ind, jy, jz) = Pi_cipar(r.ind, jy - 1, jz);
	}
      }
    }
    Pi_ciperp = fromFieldAligned(Pi_ciperp);
    Pi_cipar = fromFieldAligned(Pi_cipar);

    mesh->communicate(Pi_ciperp, Pi_cipar);

    const Field3D div_Pi_ciperp = - (2. / 3) * Grad_par(Pi_ciperp) + Pi_ciperp * Grad_par_logB;
    //const Field3D div_Pi_ciperp = - (2. / 3) * B32 * Grad_par(Pi_ciperp / B32);

    add(species["momentum_source"], div_Pi_ciperp);
    subtract(species["energy_source"], V * div_Pi_ciperp);

    // Total scalar ion viscous pressure
    Field3D Pi_ci = Pi_cipar + Pi_ciperp;

#if CHECKLEVEL >= 1
    for (auto& i : N.getRegion("RGN_NOBNDRY")) {
      if (!std::isfinite(Pi_cipar[i])) {
        throw BoutException("{} Pi_cipar non-finite at {}.\n", species_name, i);
      }
      if (!std::isfinite(Pi_ciperp[i])) {
        throw BoutException("{} Pi_ciperp non-finite at {}.\n", species_name, i);
      }
    }
#endif

    // Divergence of current in vorticity equation
    Field3D DivJ = Div(0.5 * Pi_ci * Curlb_B) -
      Div_n_bxGrad_f_B_XPPM(1. / 3, Pi_ci, false);
    add(state["fields"]["DivJextra"], DivJ);

    // Transfer of energy between ion internal energy and ExB flow
    subtract(species["energy_source"], 0.5 * Pi_ci * Curlb_B * Grad(phi + P)
             - (1 / 3) * bracket(Pi_ci, phi + P, BRACKET_ARAKAWA));

    if (diagnose) {
      // Find the diagnostics struct for this species
      auto search = diagnostics.find(species_name);
      if (search == diagnostics.end()) {
        // First time, create diagnostic
        diagnostics.emplace(species_name, Diagnostics {Pi_ciperp, Pi_cipar});
      } else {
        // Update diagnostic values
        auto& d = search->second;
        d.Pi_ciperp = Pi_ciperp;
        d.Pi_cipar = Pi_cipar;
      }
    }
  }
}

void IonViscosity::outputVars(Options &state) {
  AUTO_TRACE();

  if (diagnose) {
    // Normalisations
    auto Nnorm = get<BoutReal>(state["Nnorm"]);
    auto Tnorm = get<BoutReal>(state["Tnorm"]);
    BoutReal Pnorm = SI::qe * Tnorm * Nnorm; // Pressure normalisation

    for (const auto& it : diagnostics) {
      const std::string& species_name = it.first;
      const auto& d = it.second;

      set_with_attrs(state[std::string("P") + species_name + std::string("_ciperp")], d.Pi_ciperp,
                     {{"time_dimension", "t"},
                      {"units", "Pa"},
                      {"conversion", Pnorm},
                      {"standard_name", "Viscous pressure"},
                      {"long_name", species_name + " perpendicular collisional viscous pressure"},
                      {"species", species_name},
                      {"source", "ion_viscosity"}});

      set_with_attrs(state[std::string("P") + species_name + std::string("_cipar")], d.Pi_cipar,
                     {{"time_dimension", "t"},
                      {"units", "Pa"},
                      {"conversion", Pnorm},
                      {"standard_name", "Viscous pressure"},
                      {"long_name", species_name + " parallel collisional viscous pressure"},
                      {"species", species_name},
                      {"source", "ion_viscosity"}});
    }
  }
}
