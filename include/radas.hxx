#pragma once
#ifndef RADAS_H
#define RADAS_H

#include "component.hxx"

#include <options_netcdf.hxx>

/// Use a cooling curve as produced by RADAS
/// https://github.com/cfs-energy/radas
///
/// NetCDF files contain a 3D array
/// double noncoronal_electron_emission_prefactor(
///    dim_electron_density,
///    dim_electron_temperature,
///    dim_residence_time)
///
struct Radas : public Component {
  ///
  /// # Inputs
  ///
  ///  - <name>
  ///    - dataset: string e.g. "dataset.nc"
  ///      A NetCDF file containing the RADAS dataset
  ///
  Radas(std::string name, Options& alloptions, Solver* UNUSED(solver)) : name(name) {
    auto& options = alloptions[name];

    // Read the cooling curve from file
    auto filename = options["dataset"]
                        .doc("NetCDF file containing the RADAS dataset.")
                        .as<std::string>();
    auto data = bout::OptionsNetCDF(filename).read();

    // Extract the cooling curve as a 3D array
    Lz = data["noncoronal_electron_emission_prefactor"].as<Tensor<BoutReal>>();
    // Axes of the 3D tensor are 1D arrays
    ne_axis = data["electron_density"].as<Array<BoutReal>>();
    output.write("\tNe: {} ({} -> {})\n", ne_axis.size(), ne_axis[0],
                 ne_axis[ne_axis.size() - 1]);

    te_axis = data["electron_temperature"].as<Array<BoutReal>>();
    output.write("\tTe: {} ({} -> {})\n", te_axis.size(), te_axis[0],
                 ne_axis[te_axis.size() - 1]);

    for (int i = 0; i < te_axis.size(); ++i) {
      output.write("Te {} : {}\n", i, te_axis[i]);
    }

    tau_axis = data["residence_time"].as<Array<BoutReal>>();
    output.write("\ttau: {} ({} -> {})\n", tau_axis.size(), tau_axis[0],
                 ne_axis[tau_axis.size() - 1]);

    auto shape = Lz.shape();
    if (std::get<0>(shape) != ne_axis.size()) {
      throw BoutException("noncoronal_electron_emission_prefactor Ne axis error");
    }
    if (std::get<1>(shape) != te_axis.size()) {
      throw BoutException("noncoronal_electron_emission_prefactor Te axis error");
    }
    if (std::get<2>(shape) != tau_axis.size()) {
      throw BoutException("noncoronal_electron_emission_prefactor tau axis error");
    }
  }

  void transform(Options& state) override {}

  void outputVars(Options& state) override {}

private:
  std::string name;

  /// [ne, te, tau]   [W m^3]
  Tensor<BoutReal> Lz;

  Array<BoutReal> ne_axis;  ///< Density axis values [m^-3]
  Array<BoutReal> te_axis;  ///< Temperature axis values [eV]
  Array<BoutReal> tau_axis; ///< Residence time values [s]
};

namespace {
RegisterComponent<Radas> registercomponentradas("radas");
}

#endif // RADAS_H
