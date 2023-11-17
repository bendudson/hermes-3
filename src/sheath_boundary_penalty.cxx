
#include "../include/sheath_boundary_penalty.hxx"

#include <bout/globals.hxx>
#include <bout/mesh.hxx>
using bout::globals::mesh;

SheathBoundaryPenalty::SheathBoundaryPenalty(std::string name, Options& alloptions,
                                             Solver*) {
  AUTO_TRACE();

  Options& options = alloptions[name];

  gamma_e = options["gamma_e"]
                .doc("Electron sheath heat transmission coefficient")
                .withDefault(3.5);

  gamma_i =
      options["gamma_i"].doc("Ion sheath heat transmission coefficient").withDefault(3.5);

  std::string mask_name = options["mask_name"]
                              .doc("Name of the mesh variable containing penalty mask")
                              .withDefault<std::string>("penalty_mask");

  if (mesh->get(penalty_mask, mask_name) != 0) {
    throw BoutException("Could not read penalty mask variable '{}'", mask_name);
  }

  // Find every cell that has penalty_mask > 0
  // so we can efficiently iterate over them later
  Region<Ind3D>::RegionIndices indices;
  BOUT_FOR_SERIAL(i, penalty_mask.getRegion("RGN_NOBNDRY")) {
    if (penalty_mask[i] > 1e-5) {
      // Add this cell to the iteration
      indices.push_back(i);
    }
  }
  penalty_region = Region<Ind3D>(indices);
}

void SheathBoundaryPenalty::transform(Options& state) {
  AUTO_TRACE();

  Options& allspecies = state["species"];

  // Electrons
  auto& electrons = allspecies["e"];
  auto Ne = getNoBoundary<Field3D>(electrons["density"]);
  auto Te = getNoBoundary<Field3D>(electrons["temperature"]);
  const BoutReal Me = get<BoutReal>(electrons["AA"]);

  {
    Field3D Pe = electrons.isSet("pressure")
                     ? getNoBoundary<Field3D>(electrons["pressure"])
                     : Ne * Te;

    Field3D Ve = electrons.isSet("velocity")
                     ? getNoBoundary<Field3D>(electrons["velocity"])
                     : zeroFrom(Ne);

    Field3D NVe = electrons.isSet("momentum")
                      ? getNoBoundary<Field3D>(electrons["momentum"])
                      : Me * Ne * Ve;

    Field3D density_source = electrons.isSet("density_source")
                                 ? getNonFinal<Field3D>(electrons["density_source"])
                                 : zeroFrom(Ne);

    Field3D momentum_source = electrons.isSet("momentum_source")
                                  ? getNonFinal<Field3D>(electrons["momentum_source"])
                                  : zeroFrom(Ne);

    Field3D energy_source = electrons.isSet("energy_source")
                                ? getNonFinal<Field3D>(electrons["energy_source"])
                                : zeroFrom(Ne);

    BOUT_FOR(i, penalty_region) {
      BoutReal mask = penalty_mask[i]; // 1 in boundary
      density_source[i] = (1 - mask) * density_source[i] - mask * Ne[i];
      momentum_source[i] = (1 - mask) * momentum_source[i] - mask * NVe[i];
      energy_source[i] = (1 - mask) * energy_source[i] - mask * gamma_e * Pe[i];
    }

    set(electrons["density_source"], density_source);
    set(electrons["momentum_source"], momentum_source);
    set(electrons["energy_source"], energy_source);
  }

  for (auto& kv : allspecies.getChildren()) {
    if (kv.first == "e") {
      continue; // Skip electrons
    }
    Options& species = allspecies[kv.first]; // Note: Need non-const

    // Characteristics of this species
    const BoutReal Mi = get<BoutReal>(species["AA"]);
    const BoutReal Zi = get<BoutReal>(species["charge"]);

    if (fabs(Zi) < 1e-5) {
      // Skip neutrals
      continue;
    }

    Field3D Ni = floor(getNoBoundary<Field3D>(species["density"]), 0.0);
    Field3D Ti = getNoBoundary<Field3D>(species["temperature"]);

    Field3D Pi =
        species.isSet("pressure") ? getNoBoundary<Field3D>(species["pressure"]) : Ni * Ti;

    Field3D Vi = species.isSet("velocity") ? getNoBoundary<Field3D>(species["velocity"])
                                           : zeroFrom(Ni);

    Field3D NVi = species.isSet("momentum") ? getNoBoundary<Field3D>(species["momentum"])
                                            : Mi * Ni * Vi;

    // Get the particle source to modify
    Field3D density_source = species.isSet("density_source")
                                 ? getNonFinal<Field3D>(species["density_source"])
                                 : zeroFrom(Ni);

    Field3D momentum_source = species.isSet("momentum_source")
                                  ? getNonFinal<Field3D>(species["momentum_source"])
                                  : zeroFrom(Ni);

    Field3D energy_source = species.isSet("energy_source")
                                ? getNonFinal<Field3D>(species["energy_source"])
                                : zeroFrom(Ni);

    BOUT_FOR(i, penalty_region) {
      BoutReal mask = penalty_mask[i]; // 1 in boundary
      density_source[i] = (1 - mask) * density_source[i] - mask * Ni[i];
      momentum_source[i] = (1 - mask) * momentum_source[i] - mask * NVi[i];

      //
      BoutReal u = momentum_source[i] / (Mi * density_source[i]);

      energy_source[i] = (1 - mask) * energy_source[i] - mask * gamma_i * Pi[i];
    }

    set(species["density_source"], density_source);
    set(species["momentum_source"], momentum_source);
    set(species["energy_source"], energy_source);
  }
}
