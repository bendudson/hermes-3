#pragma once
#ifndef fieldline_geometry_H
#define fieldline_geometry_H

#include "component.hxx"
#include <bout/constants.hxx>

struct FieldlineGeometry : public Component {

  FieldlineGeometry(std::string, Options& options, Solver*);

  void transform(Options& state) override;
  void outputVars(Options& state) override;

  private:
    BoutReal Lnorm;
    bool diagnose;

    Field3D lpar{0.0};
    FieldGeneratorPtr B_poloidal_function;
    Field3D B_poloidal{0.0};

};

namespace {
  RegisterComponent<FieldlineGeometry> register_fieldline_geometry("fieldline_geometry");
}

#endif // fieldline_geometry_H
