#pragma once
#ifndef fieldline_geometry_H
#define fieldline_geometry_H

#include "component.hxx"
#include <bout/constants.hxx>

struct FieldlineGeometry : public Component {

  FieldlineGeometry(std::string name, Options& alloptions, Solver*) : name(name) {

    Options& options = alloptions[name];

    const auto& units = alloptions["units"];
    Nnorm = get<BoutReal>(units["inv_meters_cubed"]);
    Omega_ci = 1.0 / get<BoutReal>(units["seconds"]);

    residence_time = (options["residence_time"]
        .doc("Pumping time-constant. Units [s]")
        .as<BoutReal>()
    ) * Omega_ci;

    sink_shape = (options["sink_shape"]
        .doc("Shape of pumping sink.")
        .withDefault(Field3D(0.0))
    );

    diagnose = options["diagnose"]
                  .doc("Output additional diagnostics?")
                  .withDefault<bool>(false);

    };

    void transform(Options& state) override {

        Field3D species_density = getNoBoundary<Field3D>(state["species"][name]["density"]);

        pumping_sink = (sink_shape * species_density) * (-1.0 / residence_time);

        add(state["species"][name]["density_source"], pumping_sink);

    };

    void outputVars(Options& state) override {
    AUTO_TRACE();
    if (diagnose) {

      set_with_attrs(
          state[std::string("fieldline_geometry_src_shape_" + name)], sink_shape,
          {{"long_name", "simple pump source shape"},
           {"source", "fieldline_geometry"}});
    
      set_with_attrs(
          state[std::string("fieldline_geometry_sink_" + name)], pumping_sink,
          {{"time_dimension", "t"},
           {"units", "m^-3 / s"},
           {"conversion", Nnorm * Omega_ci},
           {"long_name", "simple pump source shape"},
           {"source", "fieldline_geometry"}});

    }}

  private:
    std::string name; ///< The species name
    Field3D sink_shape;
    Field3D pumping_sink;
    BoutReal Nnorm;
    BoutReal Omega_ci;
    BoutReal residence_time;
    bool diagnose;

};

namespace {
  RegisterComponent<FieldlineGeometry> register_fieldline_geometry("fieldline_geometry");
}

#endif // fieldline_geometry_H
