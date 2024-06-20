#pragma once
#ifndef fieldline_geometry_H
#define fieldline_geometry_H

#include "component.hxx"
#include <bout/constants.hxx>

struct FieldlineGeometry : public Component {

  FieldlineGeometry(std::string, Options& options, Solver*) {

    Options& geo_options = options["fieldline_geometry"];

    const auto& units = options["units"];
    Lnorm = get<BoutReal>(units["meters"]);

    parallel_length = (geo_options["parallel_length"]
      .doc("Parallel length as a function of grid index [m]")
      .withDefault(Field3D(0.0))
    ) / Lnorm;

    diagnose = geo_options["diagnose"]
                  .doc("Output additional diagnostics?")
                  .withDefault<bool>(false);

    };

    void transform(Options& state) override {

    };

    void outputVars(Options& state) override {
    AUTO_TRACE();
    if (diagnose) {

      set_with_attrs(
          state[std::string("fieldline_geometry_parallel_length")], parallel_length,
          {{"units", "m"},
           {"conversion", Lnorm},
           {"long_name", "Parallel length"},
           {"source", "fieldline_geometry"}});

      // set_with_attrs(
      //     state[std::string("fieldline_geometry_src_shape_" + name)], sink_shape,
      //     {{"long_name", "simple pump source shape"},
      //      {"source", "fieldline_geometry"}});
    
      // set_with_attrs(
      //     state[std::string("fieldline_geometry_sink_" + name)], pumping_sink,
      //     {{"time_dimension", "t"},
      //      {"units", "m^-3 / s"},
      //      {"conversion", Nnorm * Omega_ci},
      //      {"long_name", "simple pump source shape"},
      //      {"source", "fieldline_geometry"}});

    }}

  private:
    BoutReal Lnorm;
    bool diagnose;

    Field3D parallel_length;

};

namespace {
  RegisterComponent<FieldlineGeometry> register_fieldline_geometry("fieldline_geometry");
}

#endif // fieldline_geometry_H
