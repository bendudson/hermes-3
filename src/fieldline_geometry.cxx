#include "../include/fieldline_geometry.hxx"
#include <bout/constants.hxx>
#include <bout/coordinates.hxx>
#include <bout/mesh.hxx>
#include <bout/field_factory.hxx>
#include <iostream> // For outputting to log file

using bout::globals::mesh;

FieldlineGeometry::FieldlineGeometry(std::string, Options& options, Solver*) {

    Options& geo_options = options["fieldline_geometry"];

    const auto& units = options["units"];
    Lnorm = get<BoutReal>(units["meters"]);

    // lpar = (geo_options["parallel_length"]
    //     .doc("Parallel length as a function of grid index [m]")
    //     .withDefault(Field3D(0.0))
    // ) / Lnorm;

    const int MYPE = BoutComm::rank();   // Current rank
    const int NPES = BoutComm::size();    // Number of procs
    const int NYPE = NPES / mesh->NXPE;    // Number of procs in Y
    Coordinates *coord = mesh->getCoordinates();
    lpar = 0;
    auto dy = coord->dy;
    lpar(0,0,0) = 0.5 * dy(0,0,0);
    for (int id = 0; id <= NYPE-1; ++id) {   // Iterate through each proc
        for (int j = 0; j <= mesh->LocalNy; ++j) {
            if (j > 0) {
                lpar(0, j, 0) = lpar(0, j-1, 0) + 0.5 * dy(0, j-1, 0) + 0.5 * dy(0, j, 0);
            }
        }
        mesh->communicate(lpar);  // When proc is done, communicate to other ones
    }

    auto str = geo_options["B_poloidal"]
      .doc("Function for B_poloidal(lpar) [T].")
      .as<std::string>();
    B_poloidal_function = FieldFactory::get()->parse(str, &geo_options);

    for (int j = mesh->ystart; j <= mesh->yend; ++j) {
        B_poloidal(0, j, 0) = B_poloidal_function ->generate(bout::generator::Context().set("lpar", lpar(0, j, 0)));
    }

    diagnose = geo_options["diagnose"]
                    .doc("Output additional diagnostics?")
                    .withDefault<bool>(false);

}

void FieldlineGeometry::transform(Options& state) {

}

void FieldlineGeometry::outputVars(Options& state) {
    AUTO_TRACE();
    if (diagnose) {

        set_with_attrs(
            state[std::string("fieldline_geometry_parallel_length")], lpar,
            {{"units", "m"},
            // {"conversion", Lnorm},
            {"long_name", "Parallel length"},
            {"source", "fieldline_geometry"}});

}}