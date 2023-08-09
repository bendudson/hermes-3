
#include "../include/full_velocity.hxx"
#include "../include/div_ops.hxx"

#include "bout/mesh.hxx"
#include "bout/solver.hxx"
#if not BOUT_USE_METRIC_3D 
using bout::globals::mesh;

NeutralFullVelocity::NeutralFullVelocity(const std::string& name, Options& options, Solver *solver) : name(name) {

  // This is used in both transform and finally functions
  coord = mesh->getCoordinates();
  
  dy2D = DC(coord->dy);
  dx2D = DC(coord->dx);
  J2D = DC(coord->J);
  g_22_2D = DC(coord->g_22);
  Bxy2D = DC(coord->Bxy); // redundant name...
  
  AA = options["AA"].doc("Atomic mass number").withDefault(2.0);

  gamma_ratio =
      options["gamma_ratio"].doc("Ratio of specific heats").withDefault(5. / 3);

  neutral_viscosity = options["viscosity"]
                          .doc("Normalised kinematic viscosity")
                          .withDefault(1e-2);

  neutral_bulk =
      options["bulk"].doc("Normalised bulk viscosity").withDefault(1e-2);

  neutral_conduction =
      options["conduction"].doc("Normalised heat conduction").withDefault(1e-2);

  outflow_ydown = options["outflow_ydown"]
                      .doc("Allow outflowing neutrals?")
                      .withDefault<bool>(false);

  neutral_gamma = options["neutral_gamma"]
                      .doc("Surface heat transmission coefficient")
                      .withDefault(5. / 4);

  // Get the normalisations
  auto& units = Options::root()["units"];
  BoutReal Lnorm = units["meters"];
  BoutReal Bnorm = units["Tesla"];
  Tnorm = units["eV"];
  
  // Evolve 2D density, pressure, and velocity
  solver->add(Nn2D, "Nn");
  solver->add(Pn2D, "Pn");
  solver->add(Vn2D, "Vn");

  DivV2D.setBoundary("Pn"); // Same boundary condition as Pn

  // Load necessary metrics for non-orth calculation
  Field2D etaxy, cosbeta;
  if (mesh->get(etaxy, "etaxy")) {
    etaxy = 0.0;
  }
  cosbeta = sqrt(1. - SQ(etaxy));

  // Calculate transformation to Cartesian coordinates
  Field2D Rxy, Zxy, hthe, Bpxy;

  if (mesh->get(Rxy, "Rxy")) {
    throw BoutException("Fluid neutrals model requires Rxy");
  }
  if (mesh->get(Zxy, "Zxy")) {
    throw BoutException("Fluid neutrals model requires Zxy");
  }
  if (mesh->get(hthe, "hthe")) {
    throw BoutException("Fluid neutrals model requires hthe");
  }
  if (mesh->get(Bpxy, "Bpxy")) {
    throw BoutException("Fluid neutrals model requires Bpxy");
  }

  // Normalise
  Rxy /= Lnorm;
  Zxy /= Lnorm;
  hthe /= Lnorm;
  Bpxy /= Bnorm;

  // Axisymmetric neutrals simplifies things considerably...

  Urx.allocate();
  Ury.allocate();
  Uzx.allocate();
  Uzy.allocate();

  Txr.allocate();
  Txz.allocate();
  Tyr.allocate();
  Tyz.allocate();
  
  for (int i = 0; i < mesh->LocalNx; i++)
    for (int j = mesh->ystart; j <= mesh->yend; j++) {
      // Central differencing of coordinates
      BoutReal dRdtheta, dZdtheta;
      if (j == mesh->ystart) {
        dRdtheta = (Rxy(i, j + 1) - Rxy(i, j)) / (dy2D(i, j));
        dZdtheta = (Zxy(i, j + 1) - Zxy(i, j)) / (dy2D(i, j));
      } else if (j == mesh->yend) {
        dRdtheta = (Rxy(i, j) - Rxy(i, j - 1)) / (dy2D(i, j));
        dZdtheta = (Zxy(i, j) - Zxy(i, j - 1)) / (dy2D(i, j));
      } else {
        dRdtheta = (Rxy(i, j + 1) - Rxy(i, j - 1)) / (2. * dy2D(i, j));
        dZdtheta = (Zxy(i, j + 1) - Zxy(i, j - 1)) / (2. * dy2D(i, j));
      }

      // Match to hthe, 1/|Grad y|
      BoutReal h = sqrt(SQ(dRdtheta) + SQ(dZdtheta));
      BoutReal grady = 1.0 / hthe(i, j);
      dRdtheta = dRdtheta / grady / h;
      dZdtheta = dZdtheta / grady / h;

      BoutReal dRdpsi, dZdpsi;
      if (i == 0) {
        // One-sided differences
        dRdpsi = (Rxy(i + 1, j) - Rxy(i, j)) / (dx2D(i, j));
        dZdpsi = (Zxy(i + 1, j) - Zxy(i, j)) / (dx2D(i, j));
      } else if (i == (mesh->LocalNx - 1)) {
        // One-sided differences
        dRdpsi = (Rxy(i, j) - Rxy(i - 1, j)) / (dx2D(i, j));
        dZdpsi = (Zxy(i, j) - Zxy(i - 1, j)) / (dx2D(i, j));
      } else {
        dRdpsi = (Rxy(i + 1, j) - Rxy(i - 1, j)) / (2. * dx2D(i, j));
        dZdpsi = (Zxy(i + 1, j) - Zxy(i - 1, j)) / (2. * dx2D(i, j));
      }

      // Match to Bp, |Grad psi|. NOTE: this only works if
      // X and Y are orthogonal.
      BoutReal dldpsi =
          sqrt(SQ(dRdpsi) + SQ(dZdpsi)) * cosbeta(i, j); // ~ 1/(R*Bp)
      dRdpsi /= dldpsi * Bpxy(i, j) * Rxy(i, j);
      dZdpsi /= dldpsi * Bpxy(i, j) * Rxy(i, j);

      Urx(i, j) = dRdpsi;
      Ury(i, j) = dRdtheta;
      Uzx(i, j) = dZdpsi;
      Uzy(i, j) = dZdtheta;

      // Poloidal (R,Z) transformation Jacobian
      BoutReal J = dRdpsi * dZdtheta - dZdpsi * dRdtheta;

      Txr(i, j) = dZdtheta / J;
      Txz(i, j) = -dRdtheta / J;
      Tyr(i, j) = -dZdpsi / J;
      Tyz(i, j) = dRdpsi / J;
    }

  Urx.applyBoundary("neumann");
  Ury.applyBoundary("neumann");
  Uzx.applyBoundary("neumann");
  Uzy.applyBoundary("neumann");

  Txr.applyBoundary("neumann");
  Txz.applyBoundary("neumann");
  Tyr.applyBoundary("neumann");
  Tyz.applyBoundary("neumann");
}

/// Modify the given simulation state
void NeutralFullVelocity::transform(Options &state) {
  mesh->communicate(Nn2D, Vn2D, Pn2D);
  
  // Navier-Stokes for axisymmetric neutral gas profiles
  // Nn2D, Pn2D and Tn2D are unfloored
  Tn2D = Pn2D / Nn2D;

  // Floored fields, used for rate coefficients
  Field2D Nn = floor(Nn2D, 1e-8);
  Field2D Tn = floor(Tn2D, 0.01 / Tnorm);

  //////////////////////////////////////////////////////
  // 2D (X-Y) full velocity model
  //
  // Evolves density Nn2D, velocity vector Vn2D and pressure Pn2D
  //

  if (outflow_ydown) {
    // Outflowing boundaries at ydown. If flow direction is
    // into domain then zero value is set. If flow is out of domain
    // then Neumann conditions are set

    for (RangeIterator idwn = mesh->iterateBndryLowerY(); !idwn.isDone();
         idwn.next()) {

      if (Vn2D.y(idwn.ind, mesh->ystart, 0) < 0.0) {
        // Flowing out of domain
        Vn2D.y(idwn.ind, mesh->ystart - 1) = Vn2D.y(idwn.ind, mesh->ystart);
      } else {
        // Flowing into domain
        Vn2D.y(idwn.ind, mesh->ystart - 1) = -Vn2D.y(idwn.ind, mesh->ystart);
      }
      // Neumann boundary condition on X and Z components
      Vn2D.x(idwn.ind, mesh->ystart - 1) = Vn2D.x(idwn.ind, mesh->ystart);
      Vn2D.z(idwn.ind, mesh->ystart - 1) = Vn2D.z(idwn.ind, mesh->ystart);

      // Neumann conditions on density and pressure
      Nn2D(idwn.ind, mesh->ystart - 1) = Nn2D(idwn.ind, mesh->ystart);
      Pn2D(idwn.ind, mesh->ystart - 1) = Pn2D(idwn.ind, mesh->ystart);
    }
  }

  // Exchange of parallel momentum. This could be done
  // in a couple of ways, but here we use the fact that
  // Vn2D is covariant and b = e_y / (JB) to write:
  //
  // V_{||n} = b dot V_n = Vn2D.y / (JB)
  Field2D Vnpar = Vn2D.y / (J2D * Bxy2D);

  // Set values in the state
  auto& localstate = state["species"][name];
  set(localstate["density"], Nn);
  set(localstate["pressure"], Pn2D);
  set(localstate["momentum"], Vnpar * Nn * AA);
  set(localstate["velocity"], Vnpar); // Parallel velocity
  set(localstate["temperature"], Tn);
}

/// Use the final simulation state to update internal state
/// (e.g. time derivatives)
void NeutralFullVelocity::finally(const Options &state) {
  
  // Density
  ddt(Nn2D) = -Div(Vn2D, Nn2D);

  Field2D Nn2D_floor = floor(Nn2D, 1e-2);
  
  // Velocity
  ddt(Vn2D) = -Grad(Pn2D) / (AA * Nn2D_floor);

  //////////////////////////////////////////////////////
  // Momentum advection

  // Convert to cylindrical coordinates for velocity
  // advection term. This is to avoid Christoffel symbol
  // terms in curvilinear geometry
  Field2D vr = Txr * Vn2D.x + Tyr * Vn2D.y; // Grad R component
  Field2D vz = Txz * Vn2D.x + Tyz * Vn2D.y; // Grad Z component

  // Advect as scalars (no Christoffel symbols needed)
  ddt(vr) = -V_dot_Grad(Vn2D, vr);
  ddt(vz) = -V_dot_Grad(Vn2D, vz);

  // Convert back to field-aligned coordinates
  ddt(Vn2D).x += Urx * ddt(vr) + Uzx * ddt(vz);
  ddt(Vn2D).y += Ury * ddt(vr) + Uzy * ddt(vz);

  //////////////////////////////////////////////////////
  // Viscosity
  // This includes dynamic ( neutral_viscosity)
  // and bulk/volume viscosity ( neutral_bulk )

  ddt(vr) = Laplace_FV(neutral_viscosity, vr);
  ddt(vz) = Laplace_FV(neutral_viscosity, vz);

  ddt(Vn2D).x += Urx * ddt(vr) + Uzx * ddt(vz);
  ddt(Vn2D).y += Ury * ddt(vr) + Uzy * ddt(vz);

  DivV2D = Div(Vn2D);
  DivV2D.applyBoundary(0.0);
  mesh->communicate(DivV2D);

  // ddt(Vn2D) += Grad( (neutral_viscosity/3. + neutral_bulk) * DivV2D ) /
  // Nn2D_floor;

  //////////////////////////////////////////////////////
  // Pressure
  ddt(Pn2D) = -Div(Vn2D, Pn2D) -
              (gamma_ratio - 1.) * Pn2D * DivV2D * floor(Nn2D, 0) / Nn2D_floor +
              Laplace_FV(neutral_conduction, Pn2D / Nn2D);

  ///////////////////////////////////////////////////////////////////
  // Boundary condition on fluxes

  for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {

    // Loss of thermal energy to the target.
    // This depends on the reflection coefficient
    // and is controlled by the option neutral_gamma
    //         q = neutral_gamma * n * T * cs

    // Density at the target
    BoutReal Nnout =
        0.5 * (Nn2D(r.ind, mesh->ystart) + Nn2D(r.ind, mesh->ystart - 1));
    if (Nnout < 0.0)
      Nnout = 0.0;
    // Temperature at the target
    BoutReal Tnout =
        0.5 * (Tn2D(r.ind, mesh->ystart) + Tn2D(r.ind, mesh->ystart - 1));
    if (Tnout < 0.0)
      Tnout = 0.0;

    // gamma * n * T * cs
    BoutReal q = neutral_gamma * Nnout * Tnout * sqrt(Tnout);
    // Multiply by cell area to get power
    BoutReal heatflux =
        q * (J2D(r.ind, mesh->ystart) + J2D(r.ind, mesh->ystart - 1)) /
        (sqrt(g_22_2D(r.ind, mesh->ystart)) +
         sqrt(g_22_2D(r.ind, mesh->ystart - 1)));

    // Divide by volume of cell, and multiply by 2/3 to get pressure
    ddt(Pn2D)(r.ind, mesh->ystart) -=
        (2. / 3) * heatflux /
        (dy2D(r.ind, mesh->ystart) * J2D(r.ind, mesh->ystart));
  }

  for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {

    // Loss of thermal energy to the target.
    // This depends on the reflection coefficient
    // and is controlled by the option neutral_gamma
    //         q = neutral_gamma * n * T * cs

    // Density at the target
    BoutReal Nnout =
        0.5 * (Nn2D(r.ind, mesh->yend) + Nn2D(r.ind, mesh->yend + 1));
    if (Nnout < 0.0)
      Nnout = 0.0;
    // Temperature at the target
    BoutReal Tnout =
        0.5 * (Tn2D(r.ind, mesh->yend) + Tn2D(r.ind, mesh->yend + 1));
    if (Tnout < 0.0)
      Tnout = 0.0;

    // gamma * n * T * cs
    BoutReal q = neutral_gamma * Nnout * Tnout * sqrt(Tnout);
    // Multiply by cell area to get power
    BoutReal heatflux =
        q * (J2D(r.ind, mesh->yend) + J2D(r.ind, mesh->yend + 1)) /
        (sqrt(g_22_2D(r.ind, mesh->yend)) +
         sqrt(g_22_2D(r.ind, mesh->yend + 1)));

    // Divide by volume of cell, and multiply by 2/3 to get pressure
    ddt(Pn2D)(r.ind, mesh->yend) -=
        (2. / 3) * heatflux /
        (dy2D(r.ind, mesh->yend) * J2D(r.ind, mesh->yend));
  }

  /////////////////////////////////////////////////////
  // Atomic processes
  
  auto& localstate = state["species"][name];

  // Particles
  if (localstate.isSet("density_source")) {
    ddt(Nn2D) += get<Field2D>(localstate["density_source"]);
  }

  // Momentum. Note need to turn back into covariant form
  if (localstate.isSet("momentum_source")) {
    ddt(Vn2D).y += get<Field2D>(localstate["momentum_source"])
      * (J2D * Bxy2D) / (AA * Nn2D_floor);
  }

  // Energy
  if (localstate.isSet("energy_source")) {
    ddt(Pn2D) += (2. / 3) * get<Field2D>(localstate["energy_source"]);
  }

  // Density evolution
  for (auto &i : Nn2D.getRegion("RGN_ALL")) {
    if ((Nn2D[i] < 1e-8) && (ddt(Nn2D)[i] < 0.0)) {
      ddt(Nn2D)[i] = 0.0;
    }
  }
}

/// Add extra fields for output, or set attributes e.g docstrings
void NeutralFullVelocity::outputVars(Options &state) {
  // Normalisations
  auto Nnorm = get<BoutReal>(state["Nnorm"]);
  auto Tnorm = get<BoutReal>(state["Tnorm"]);
  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
  auto Cs0 = get<BoutReal>(state["Cs0"]);

  set_with_attrs(state["DivV2D"], DivV2D, {
      {"time_dimension", "t"},
      {"units", "s^-1"},
      {"conversion", Omega_ci}
    });

  set_with_attrs(state["Urx"], Urx, {});
  set_with_attrs(state["Ury"], Ury, {});
  set_with_attrs(state["Uzx"], Uzx, {});
  set_with_attrs(state["Uzy"], Uzy, {});
  set_with_attrs(state["Txr"], Txr, {});
  set_with_attrs(state["Txz"], Txz, {});
  set_with_attrs(state["Tyr"], Tyr, {});
  set_with_attrs(state["Tyz"], Tyz, {});
}

#endif
