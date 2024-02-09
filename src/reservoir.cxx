
#include "../include/reservoir.hxx"
#include <bout/coordinates.hxx>
#include <bout/mesh.hxx>
#include <bout/constants.hxx>
#include "../include/hermes_utils.hxx"  // For indexAt
#include <iostream>   // For outputting to log file

using bout::globals::mesh;

// This is a simple reservoir model. It's meant for 1D but will work in 2D and 3D as well.
// Define the reservoir location by using an analytical expression in the input file
// e.g. reservoir_location = H(y - y_xpt) where y_xpt is the x-point location would create
// one between the target and the X-point. 
// The reservoir is set to a constant density and your chosen species will equilibrate with it
// based on the density difference and a provided timescale. 
// The intended use is to provide some form of cross-field transport for neutrals in 1D.

Reservoir::Reservoir(std::string name, Options& alloptions, Solver*) 
    :name(name){
  AUTO_TRACE();

  const Options& units = alloptions["units"];
  const BoutReal Tnorm = units["eV"];
  const BoutReal Nnorm = units["inv_meters_cubed"];
  const BoutReal Omega_ci = 1. / units["seconds"].as<BoutReal>();

  

  // Get the options for this species
  Options& options = alloptions[name];

  reservoir_density = options["reservoir_density"]
                  .doc("Set the density of the reservoir in [m^-3]. Default 1e19")
                  .withDefault<BoutReal>(1e19)
                  /Nnorm;    

  reservoir_location = options["reservoir_location"]
                  .doc("Indicates reservoir location if >0")
                  .withDefault<Field3D>(0.0);  

  reservoir_timescale = options["reservoir_timescale"]
                  .doc("Set the timescale of the reservoir in [s]. Default 1e-6")
                  .withDefault<BoutReal>(1e-6)
                  *Omega_ci;

  diagnose =
      options["diagnose"].doc("Save additional diagnostics?").withDefault<bool>(false);
          
}


void Reservoir::transform(Options& state) {
  AUTO_TRACE();

  // We are operating on only one species
  auto& species = state["species"][name];


  // Get metric tensor components
  Coordinates* coord = mesh->getCoordinates();
  const Field2D& J = coord->J;
  const Field2D& dy = coord->dy;
  const Field2D& dx = coord->dx;
  const Field2D& dz = coord->dz;

  const Field2D dv = dx * dy * dz * J;

  // Get state sources - we will then add to them and pass them back to the state so they get added to the RHS
  // state_density_source = species.isSet("density_source")
  //                               ? getNonFinal<Field3D>(species["density_source"])
  //                               : 0.0;
  // state_energy_source = species.isSet("energy_source")
  //                             ? getNonFinal<Field3D>(species["energy_source"])
  //                             : 0.0;
  // state_momentum_source = species.isSet("momentum_source")
  //                             ? getNonFinal<Field3D>(species["momentum_source"])
  //                             : 0.0;

  // These are the sources we are computing which we will add to state sources at the end
  density_source = 0;
  energy_source = 0;
  momentum_source = 0;

  // Get conditions
  Field3D N = get<Field3D>(species["density"]);
  Field3D P = get<Field3D>(species["pressure"]);
  Field3D NV = get<Field3D>(species["momentum"]);


  for(int ix=0; ix < mesh->LocalNx ; ix++){
      for(int iy=0; iy < mesh->LocalNy ; iy++){
          for(int iz=0; iz < mesh->LocalNz; iz++){

            auto i = indexAt(N, ix, iy, iz);

            // Particle transfer rate proportional to difference in density over timescale
            BoutReal Nrate = (reservoir_density - N[i]) / reservoir_timescale;

            // Pressure and momentum flows proportional to density
            // When flow is reversed and these become sources, the new particles have 
            // the same pressure and momentum as the local particles.
            BoutReal Prate = P[i] / N[i] * Nrate;
            BoutReal NVrate = NV[i] / N[i] * Nrate;



            if (reservoir_location[i] > 0) {

              density_source[i] += Nrate;
              energy_source[i] += (3. / 2) * Prate;
              momentum_source[i] += NVrate;
            }

          }
      }
  }

  // Put the updated sources back into the state
  // set<Field3D>(species["density_source"], state_density_source + density_source);
  // set<Field3D>(species["energy_source"], state_energy_source + energy_source);
  // set<Field3D>(species["momentum_source"], state_momentum_source + momentum_source);

  // Field3D neutral_density1 = getNoBoundary<Field3D>(species["density"]);

  // for(int ix=0; ix < mesh->LocalNx ; ix++){
  //     for(int iy=0; iy < mesh->LocalNy ; iy++){
  //         for(int iz=0; iz < mesh->LocalNz; iz++){

  //           // output << "("" << ix << "Y:" << iy << "Z:" << iz << "T:" << Tn(ix, iy, iz) << "  ";
  //           std::string string_count = std::string("(") + std::to_string(ix) + std::string(",") + std::to_string(iy)+ std::string(",") + std::to_string(iz) + std::string(")");
  //           output << string_count + std::string(": ") + std::to_string(neutral_density1(ix,iy,iz)) + std::string("; ");
  //           output << "\n";
  //         }
  //     }
  //   output << "\n";
  //   }

  add(species["density_source"], density_source);
  add(species["energy_source"], energy_source);
  add(species["momentum_source"], momentum_source);

  // for(int ix=0; ix < mesh->LocalNx ; ix++){
  //     for(int iy=0; iy < mesh->LocalNy ; iy++){
  //         for(int iz=0; iz < mesh->LocalNz; iz++){

  //           // output << "("" << ix << "Y:" << iy << "Z:" << iz << "T:" << Tn(ix, iy, iz) << "  ";
  //           std::string string_count = std::string("(") + std::to_string(ix) + std::string(",") + std::to_string(iy)+ std::string(",") + std::to_string(iz) + std::string(")");
  //           output << string_count + std::string(": ") + std::to_string(density_source(ix,iy,iz)) + std::string("; ");
  //           output << "\n";
  //         }
  //     }
  //   output << "\n";
  //   }

  // Field3D neutral_density2 = getNoBoundary<Field3D>(species["density"]);

  // for(int ix=0; ix < mesh->LocalNx ; ix++){
  //     for(int iy=0; iy < mesh->LocalNy ; iy++){
  //         for(int iz=0; iz < mesh->LocalNz; iz++){

  //           // output << "("" << ix << "Y:" << iy << "Z:" << iz << "T:" << Tn(ix, iy, iz) << "  ";
  //           std::string string_count = std::string("(") + std::to_string(ix) + std::string(",") + std::to_string(iy)+ std::string(",") + std::to_string(iz) + std::string(")");
  //           output << string_count + std::string(": ") + std::to_string(neutral_density2(ix,iy,iz)) + std::string("; ");
  //           output << "\n";
  //         }
  //     }
  //   output << "\n";
  //   }

}

void Reservoir::outputVars(Options& state) {
  
  AUTO_TRACE();
  // Normalisations
  auto Nnorm = get<BoutReal>(state["Nnorm"]);
  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
  auto Tnorm = get<BoutReal>(state["Tnorm"]);
  auto Cs0 = get<BoutReal>(state["Cs0"]);
  BoutReal Pnorm = SI::qe * Tnorm * Nnorm; // Pressure normalisation

  if (diagnose) {

      set_with_attrs(state[{std::string("reservoir_location_") + name}], reservoir_location,
                              {{"standard_name", name + std::string("reservoir location")},
                              {"long_name", name + std::string("reservoir location")},
                              {"source", "reservoir"}});

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

      set_with_attrs(state[{std::string("F") + name + std::string("_rsv")}], momentum_source,
                              {{"time_dimension", "t"},
                              {"units", "kg m^-2 s^-2"},
                              {"conversion", SI::Mp * Nnorm * Cs0 * Omega_ci},
                              {"standard_name", "momentum transfer"},
                              {"long_name", name + std::string(" momentum transfer from reservoir")},
                              {"source", "reservoir"}});

  };


}



