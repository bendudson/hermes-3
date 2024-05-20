#include "../include/tokamak_core_boundary.hxx"

#include <bout/mesh.hxx>
#include <bout/smoothing.hxx>
using bout::globals::mesh;

void TokamakCoreBoundary::transform(Options& state) {
  AUTO_TRACE();

  if (!mesh->firstX() || !mesh->periodicY(mesh->xstart)) {
    return; // This processor is not on the core boundary
  }

  // Operates on every species
  Options& allspecies = state["species"];
  for (auto& kv : allspecies.getChildren()) {
    Options& species = allspecies[kv.first]; // Note: need non-const

    if (IS_SET_NOBOUNDARY(species["density"])) {
      const Field3D N = GET_NOBOUNDARY(Field3D, species["density"]);

      // Flux-surface average in core region
      const Field2D N_avg = averageY(DC(N));

      Field3D density_source = (species.isSet("density_source"))
        ? getNonFinal<Field3D>(species["density_source"])
        : zeroFrom(N);

      // Add source/sink that damps towards average value
      const int x = mesh->xstart;
      for (int y = mesh->ystart; y <= mesh->yend; y++) {
        for (int z = mesh->zstart; z <= mesh->zend; z++) {
          density_source(x, y, z) += damping_rate * (N_avg(x, y) - N(x, y, z));
        }
      }
      set(species["density_source"], density_source);
    }

    if (IS_SET_NOBOUNDARY(species["pressure"])) {
      const Field3D P = GET_NOBOUNDARY(Field3D, species["pressure"]);

      // Flux-surface average in core region
      const Field2D P_avg = averageY(DC(P));

      Field3D energy_source = (species.isSet("energy_source"))
        ? getNonFinal<Field3D>(species["energy_source"])
        : zeroFrom(P);

      // Add source/sink that damps towards average value
      const int x = mesh->xstart;
      for (int y = mesh->ystart; y <= mesh->yend; y++) {
        for (int z = mesh->zstart; z <= mesh->zend; z++) {
          energy_source(x, y, z) += damping_rate * (3./2) * (P_avg(x, y) - P(x, y, z));
        }
      }
      set(species["energy_source"], energy_source);
    }
  }
}
