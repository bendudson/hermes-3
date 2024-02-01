
#include "../include/reservoir.hxx"
#include <bout/coordinates.hxx>
#include <bout/mesh.hxx>
#include <bout/constants.hxx>

using bout::globals::mesh;

Reservoir::Reservoir(std::string name, Options& alloptions, Solver*) {
  AUTO_TRACE();

  const Options& units = alloptions["units"];
  const BoutReal Tnorm = units["eV"];
  const BoutReal Nnorm = units["inv_meters_cubed"];

  // Get the options for this species
  Options& options = alloptions[name];

  reservoir_density = options["reservoir_density"]
                  .doc("Set the density of the reservoir in [m^-3]")
                  .withDefault<BoutReal>(1e19)
                  /Nnorm;    

  reservoir_location = options["reservoir_location"]
                  .doc("Indicates reservoir location if >0")
                  .withDefault(0.0);  

  diagnose =
      options["diagnose"].doc("Save additional diagnostics?").withDefault<bool>(false);
          
}


void Reservoir::transform(Options& state) {
  AUTO_TRACE();

  // Get metric tensor components
  Coordinates* coord = mesh->getCoordinates();
  const Field2D& J = coord->J;
  const Field2D& dy = coord->dy;
  const Field2D& dx = coord->dx;
  const Field2D& dz = coord->dz;

  const Field2D dv = dx * dy * dz * J;



}

void Reservoir::outputVars(Options& state) {
  
  AUTO_TRACE();
  // Normalisations
  // auto Nnorm = get<BoutReal>(state["Nnorm"]);
  // auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
  // auto Tnorm = get<BoutReal>(state["Tnorm"]);
  // BoutReal Pnorm = SI::qe * Tnorm * Nnorm; // Pressure normalisation

  if (diagnose) {

      set_with_attrs(state[{std::string("reservoir_location_") + name}], reservoir_location,
                              {{"standard_name", name + "reservoir location"},
                              {"long_name", name + "reservoir location"},
                              {"source", "reservoir"}}
      );
  };


}



