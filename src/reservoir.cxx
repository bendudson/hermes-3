
#include "../include/reservoir.hxx"
#include "../include/hermes_utils.hxx" // For indexAt
#include <bout/constants.hxx>
#include <bout/coordinates.hxx>
#include <bout/mesh.hxx>
#include <iostream> // For outputting to log file

using bout::globals::mesh;

Reservoir::Reservoir(std::string name, Options& alloptions, Solver*) : name(name) {
  AUTO_TRACE();

  const Options& units = alloptions["units"];
  const BoutReal Nnorm = units["inv_meters_cubed"];
  const BoutReal Omega_ci = 1. / units["seconds"].as<BoutReal>();

  // Get the options for this species
  Options& options = alloptions[name];

  reservoir_density = options["reservoir_density"]
                          .doc("Set the density of the reservoir in [m^-3]. Default 1e19")
                          .withDefault<BoutReal>(1e19)
                      / Nnorm;

  reservoir_location = options["reservoir_location"]
                           .doc("Indicates reservoir location if >0")
                           .withDefault<Field3D>(0.0);

  reservoir_timescale =
      options["reservoir_timescale"]
          .doc("Set the timescale of the reservoir in [s]. Default 1e-6")
          .withDefault<BoutReal>(1e-6)
      * Omega_ci;

  diagnose =
      options["diagnose"].doc("Save additional diagnostics?").withDefault<bool>(false);

  reservoir_sink_only = 
    options["reservoir_sink_only"].doc("Set reservoir to only take particles away?").withDefault<bool>(true);


  xpoint_position =
      options["xpoint_position"]
          .doc("Parallel position of X-point [m]")
          .withDefault<BoutReal>(0);

  baffle_position =
      options["baffle_position"]
          .doc("Parallel position of border between upstream-SOL and divertor-SOL reservoirs [m]")
          .withDefault<BoutReal>(xpoint_position);


  // Get lpar
  const int MYPE = BoutComm::rank();   // Current rank
  const int NPES = BoutComm::size();    // Number of procs
  const int NYPE = NPES / mesh->NXPE;    // Number of procs in Y
  Coordinates *coord = mesh->getCoordinates();

  lpar = 0;
  BoutReal offset = 0;   // Offset to ensure ylow domain boundary starts at 0
  auto dy = coord->dy;

  lpar(0,0,0) = 0.5 * dy(0,0,0);
  for (int id = 0; id <= NYPE-1; ++id) {   // Iterate through each proc
    for (int j = 0; j <= mesh->LocalNy; ++j) {
          if (j > 0) {
            lpar(0, j, 0) = lpar(0, j-1, 0) + 0.5 * dy(0, j-1, 0) + 0.5 * dy(0, j, 0);
          }
    }
    mesh->communicate(lpar);  // Communicate across guard cells so other procs keep counting
    
    if (MYPE == 0) {
      offset = (lpar(0,1,0) + lpar(0,2,0)) / 2;  // Offset lives on first proc
    }
  }

  MPI_Bcast(&offset, 1, MPI_DOUBLE, 0, BoutComm::get());  // Ensure all procs get offset
  lpar -= offset;

  // Find every cell that has reservoir_location > 0
  // so we can efficiently iterate over them later
  Region<Ind3D>::RegionIndices indices;
  BOUT_FOR_SERIAL(i, reservoir_location.getRegion("RGN_NOBNDRY")) {
    if (reservoir_location[i] > 0.0) {
      // Add this cell to the iteration
      indices.push_back(i);
    }
  }
  reservoir_region = Region<Ind3D>(indices);
}

void Reservoir::transform(Options& state) {
  AUTO_TRACE();

  // We are operating on only one species
  auto& species = state["species"][name];

  // These are the sources we are computing which we will add to state sources at the end
  density_source = 0;
  energy_source = 0;
  momentum_source = 0;

  // Get conditions. Boundary conditions do not need to be set
  // because we don't use the state in the boundary cells.
  Field3D N = GET_NOBOUNDARY(Field3D, species["density"]);
  Field3D P = GET_NOBOUNDARY(Field3D, species["pressure"]);
  Field3D NV = GET_NOBOUNDARY(Field3D, species["momentum"]);

  // Iterate over the cells where reservoir_location > 0
  BOUT_FOR(i, reservoir_region) {
    // Particle transfer rate proportional to difference in density over timescale
    BoutReal Nrate = (reservoir_density - N[i]) / reservoir_timescale;

    if (reservoir_sink_only && Nrate > 0) {
      // If we only want to remove particles, set the rate to 0
      Nrate = 0;
    };


    // Pressure and momentum flows proportional to density
    // When flow is reversed and these become sources, the new particles have
    // the same pressure and momentum as the local particles.
    BoutReal Prate = P[i] / N[i] * Nrate;
    BoutReal NVrate = NV[i] / N[i] * Nrate;

    density_source[i] += Nrate;
    energy_source[i] += (3. / 2) * Prate;
    momentum_source[i] += NVrate;
  }

  add(species["density_source"], density_source);
  add(species["energy_source"], energy_source);
  add(species["momentum_source"], momentum_source);
}

void Reservoir::outputVars(Options& state) {
  AUTO_TRACE();

  // Normalisations
  auto Nnorm = get<BoutReal>(state["Nnorm"]);
  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
  auto Tnorm = get<BoutReal>(state["Tnorm"]);
  auto Cs0 = get<BoutReal>(state["Cs0"]);
  BoutReal Pnorm = SI::qe * Tnorm * Nnorm; // Pressure normalisation

  // Always save reservoir location (time independent)
  set_with_attrs(state[{std::string("reservoir_location_") + name}], reservoir_location,
                 {{"standard_name", name + std::string("reservoir location")},
                  {"long_name", name + std::string("reservoir location")},
                  {"source", "reservoir"}});
  if (diagnose) {
    set_with_attrs(state[{std::string("S") + name + std::string("_rsv")}], density_source,
                   {{"time_dimension", "t"},
                    {"units", "m^-3 s^-1"},
                    {"conversion", Nnorm * Omega_ci},
                    {"standard_name", "particle transfer"},
                    {"long_name", name + std::string(" density transfer from reservoir")},
                    {"source", "reservoir"}});

    set_with_attrs(state[{std::string("E") + name + std::string("_rsv")}], energy_source,
                   {{"time_dimension", "t"},
                    {"units", "W m^-3"},
                    {"conversion", Pnorm * Omega_ci},
                    {"standard_name", "energy transfer"},
                    {"long_name", name + std::string(" energy transfer from reservoir")},
                    {"source", "reservoir"}});

    set_with_attrs(
        state[{std::string("F") + name + std::string("_rsv")}], momentum_source,
        {{"time_dimension", "t"},
         {"units", "kg m^-2 s^-2"},
         {"conversion", SI::Mp * Nnorm * Cs0 * Omega_ci},
         {"standard_name", "momentum transfer"},
         {"long_name", name + std::string(" momentum transfer from reservoir")},
         {"source", "reservoir"}});

    
    set_with_attrs(    // Doesn't need unnormalising for some reason
          state[std::string("lpar")], lpar,
          {{"units", "m"},
          {"standard_name", "Parallel distance from midplane"},
          {"long_name", "Parallel distance from midplane"},
          {"source", "reservoir"}});


  }
}
