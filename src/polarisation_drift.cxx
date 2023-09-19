#include "../include/polarisation_drift.hxx"

#include <bout/constants.hxx>
#include <bout/fv_ops.hxx>
#include <bout/output_bout_types.hxx>
#include <bout/invert_laplace.hxx>
#include <bout/difops.hxx>

using bout::globals::mesh;

PolarisationDrift::PolarisationDrift(std::string name,
                                     Options &alloptions,
                                     Solver *UNUSED(solver)) {
  AUTO_TRACE();

  // Get options for this component
  auto& options = alloptions[name];

  // Cache the B^2 value
  auto coord = mesh->getCoordinates();
  Bsq = SQ(coord->Bxy);

  phiSolver = Laplacian::create(&options["laplacian"]);

  // For zero polarisation drift fluxes through the boundary the
  // radial electric field at the boundary should be constant
  // (e.g. zero), so zero-gradient radial BC on phi_G.

  phiSolver->setInnerBoundaryFlags(0);
  phiSolver->setOuterBoundaryFlags(0);

  boussinesq = options["boussinesq"]
    .doc("Assume a uniform mass density in calculating the polarisation drift")
    .withDefault<bool>(true);

  if (boussinesq) {
    average_atomic_mass =
      options["average_atomic_mass"]
      .doc("Weighted average atomic mass, for polarisaion current "
           "(Boussinesq approximation)")
      .withDefault<BoutReal>(2.0); // Deuterium
  } else {
    average_atomic_mass = 1.0;
    // Use a density floor to prevent divide-by-zero errors
    density_floor = options["density_floor"].doc("Minimum density floor").withDefault(1e-5);
  }

  advection = options["advection"]
    .doc("Include advection by polarisation drift, using potential flow approximation?")
    .withDefault<bool>(true);

  diamagnetic_polarisation =
      options["diamagnetic_polarisation"]
          .doc("Include diamagnetic drift in polarisation current?")
          .withDefault<bool>(true);

  diagnose = options["diagnose"]
    .doc("Output additional diagnostics?")
    .withDefault<bool>(false);
}

void PolarisationDrift::transform(Options &state) {
  AUTO_TRACE();

  // Iterate through all subsections
  Options& allspecies = state["species"];

  // Calculate divergence of all currents except the polarisation current

  if (IS_SET(state["fields"]["DivJdia"])) {
    DivJ = get<Field3D>(state["fields"]["DivJdia"]);
  } else {
    DivJ = 0.0;
  }

  // Parallel current due to species parallel flow
  for (auto& kv : allspecies.getChildren()) {
    const Options& species = kv.second;

    if (!species.isSet("charge") or !species.isSet("momentum")) {
      continue; // Not charged, or no parallel flow
    }
    const BoutReal Z = get<BoutReal>(species["charge"]);
    if (fabs(Z) < 1e-5) {
      continue; // Not charged
    }

    const Field3D NV = GET_VALUE(Field3D, species["momentum"]);
    const BoutReal A = get<BoutReal>(species["AA"]);

    // Note: Using NV rather than N*V so that the cell boundary flux is correct
    DivJ += Div_par((Z / A) * NV);
  }

  if (diamagnetic_polarisation) {
    // Compression of ion diamagnetic contribution to polarisation velocity
    // Note: For now this ONLY includes the parallel and diamagnetic current terms
    //       Other terms e.g. ion viscous current, are in their separate components
    if (!boussinesq) {
      throw BoutException("diamagnetic_polarisation not implemented for non-Boussinesq");
    }

    // Calculate energy exchange term nonlinear in pressure
    // (3 / 2) ddt(Pi) += (Pi / n0) * Div((Pe + Pi) * Curlb_B + Jpar);
    for (auto& kv : allspecies.getChildren()) {
      Options& species = allspecies[kv.first]; // Note: need non-const

      if (!(IS_SET_NOBOUNDARY(species["pressure"]) and species.isSet("charge")
            and species.isSet("AA"))) {
        // No pressure, charge or mass -> no polarisation current due to
        // diamagnetic flow
        continue;
      }

      const auto charge = get<BoutReal>(species["charge"]);
      if (fabs(charge) < 1e-5) {
        // No charge
        continue;
      }

      const auto P = GET_NOBOUNDARY(Field3D, species["pressure"]);
      const auto AA = get<BoutReal>(species["AA"]);

      add(species["energy_source"],
          P * (AA / average_atomic_mass / charge) * DivJ);
    }
  }

  if (!advection) {
    return;
  }
  // Calculate advection terms using a potential-flow approximation

  // Calculate the total mass density of species
  // which contribute to polarisation current
  Field3D mass_density;
  if (boussinesq) {
    mass_density = average_atomic_mass;
  } else {
    mass_density = 0.0;
    for (auto& kv : allspecies.getChildren()) {
      const Options& species = kv.second;

      if (!(species.isSet("charge") and species.isSet("AA"))) {
        continue; // No charge or mass -> no current
      }
      if (fabs(get<BoutReal>(species["charge"])) < 1e-5) {
        continue; // Zero charge
      }

      const BoutReal A = get<BoutReal>(species["AA"]);
      const Field3D N = GET_NOBOUNDARY(Field3D, species["density"]);
      mass_density += A * N;
    }

    // Apply a floor to prevent divide-by-zero errors
    mass_density = floor(mass_density, density_floor);
  }
  
  if (IS_SET(state["fields"]["DivJextra"])) {
    DivJ += get<Field3D>(state["fields"]["DivJextra"]);
  }
  if (IS_SET(state["fields"]["DivJcol"])) {
    DivJ += get<Field3D>(state["fields"]["DivJcol"]);
  }

  // Solve for time derivative of potential
  // Using Div(mass_density / B^2 Grad_perp(dphi/dt)) = DivJ

  phiSolver->setCoefC(mass_density / Bsq);

  // Calculate time derivative of generalised potential
  // The assumption is that the polarisation drift can be parameterised
  // as a potential flow
  phi_pol = phiSolver->solve(DivJ * Bsq / mass_density);

  // Ensure that potential is set in communication guard cells
  mesh->communicate(phi_pol);

  // Zero flux
  phi_pol.applyBoundary("neumann");

  // Polarisation drift is given by
  //
  // v_p = - (m_i / (Z_i * B^2)) * Grad(phi_pol)
  //
  for (auto& kv : allspecies.getChildren()) {
    Options& species = allspecies[kv.first]; // Note: need non-const

    if (!(species.isSet("charge") and species.isSet("AA"))) {
      continue; // No charge or mass -> no current
    }
    const BoutReal Z = get<BoutReal>(species["charge"]);
    if (fabs(Z) < 1e-5) {
      continue; // Not charged
    }
    const BoutReal A = get<BoutReal>(species["AA"]);

    // Shared coefficient in polarisation velocity
    Field3D coef = (A / Z) / Bsq;

    if (IS_SET(species["density"])) {
      auto N = GET_VALUE(Field3D, species["density"]);
      add(species["density_source"], FV::Div_a_Grad_perp(N * coef, phi_pol));
    }

    if (IS_SET(species["pressure"])) {
      auto P = GET_VALUE(Field3D, species["pressure"]);
      add(species["energy_source"], (5. / 2) * FV::Div_a_Grad_perp(P * coef, phi_pol));
    }

    if (IS_SET(species["momentum"])) {
      auto NV = GET_VALUE(Field3D, species["momentum"]);
      add(species["momentum_source"], FV::Div_a_Grad_perp(NV * coef, phi_pol));
    }
  }
}

void PolarisationDrift::outputVars(Options &state) {
  AUTO_TRACE();
  // Normalisations
  auto Nnorm = get<BoutReal>(state["Nnorm"]);
  auto Tnorm = get<BoutReal>(state["Tnorm"]);
  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);

  if (diagnose) {
    set_with_attrs(state["DivJpol"], -DivJ,
                   {{"time_dimension", "t"},
                    {"units", "A m^-3"},
                    {"conversion", SI::qe * Nnorm * Omega_ci},
                    {"long_name", "Divergence of polarisation current"},
                    {"source", "polarisation_drift"}});

    set_with_attrs(state["phi_pol"], phi_pol,
                   {{"time_dimension", "t"},
                    {"units", "V / s"},
                    {"conversion", Tnorm * Omega_ci},
                    {"standard_name", "flow potential"},
                    {"long_name", "polarisation flow potential"},
                    {"source", "polarisation_drift"}});
  }
}
