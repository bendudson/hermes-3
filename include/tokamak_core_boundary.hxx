#pragma once
#ifndef TOKAMAK_CORE_BOUNDARY_H
#define TOKAMAK_CORE_BOUNDARY_H

#include "component.hxx"

/// Add damping of poloidal variations on core boundary
struct TokamakCoreBoundary : public Component {
  /// Options
  ///
  /// - <name>
  ///   - damping_timescale [seconds]
  ///
  TokamakCoreBoundary(std::string name, Options& alloptions, Solver* UNUSED(solver)) {
    AUTO_TRACE();

    Options& options = alloptions[name];

    damping_rate = get<BoutReal>(alloptions["units"]["seconds"]) /
      options["damping_timescale"]
      .doc("Timescale for core boundary relaxation [seconds]")
      .withDefault(1e-6);

    output.write("\tNormalised core boundary damping rate: {}", damping_rate);
  }

  /// Sets density and energy sources for all species
  /// if they have density and pressure set
  void transform(Options& state) override;

private:
  BoutReal damping_rate;
};

namespace {
RegisterComponent<TokamakCoreBoundary> registercomponenttokamakcoreboundary("tokamak_core_boundary");
}

#endif // TOKAMAK_CORE_BOUNDARY_H
