#pragma once
#ifndef RESERVOIR_H
#define RESERVOIR_H

#include "component.hxx"

/// Simple 1D reservoir model
struct Reservoir : public Component {
  
  Reservoir(std::string name, Options &alloptions, Solver *);

  void transform(Options &state) override;
  void outputVars(Options &state) override;

private:

  std::string name;                ///< Short name of the species e.g. h+
  Field3D reservoir_location;      ///< Indicates reservoir if >0
  BoutReal reservoir_density, reservoir_timescale;
  bool diagnose;

  Field3D state_density_source, state_energy_source, state_momentum_source;
  Field3D density_source, energy_source, momentum_source;
  Field3D N, P, NV;

};

namespace {
RegisterComponent<Reservoir> registercomponentreservoir("reservoir");
}

#endif // RESERVOIR_H
