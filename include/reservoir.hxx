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
  BoutReal reservoir_density;
  bool diagnose;

};

namespace {
RegisterComponent<Reservoir> registercomponentreservoir("reservoir");
}

#endif // RESERVOIR_H
