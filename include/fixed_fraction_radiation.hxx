#pragma once
#ifndef FIXED_FRACTION_RADIATION_H
#define FIXED_FRACTION_RADIATION_H

#include <bout/constants.hxx>
#include "component.hxx"

namespace {
  /// Carbon in coronal equilibrium
  /// From I.H.Hutchinson Nucl. Fusion 34 (10) 1337 - 1348 (1994)
  ///
  /// Input: Electron temperature [eV]
  /// Output: Radiation cooling curve [Wm^3]
  ///         (multiply by Ne * Ni to get W/m^3)
  struct HutchinsonCarbon {
    BoutReal curve(BoutReal Te) {
      if (Te < 0.0) {
        return 0.0;
      }
      return 2e-31*pow(Te/10., 3) / (1. + pow(Te/10., 4.5));
    }
  };

  /// Below fits are from ADAS data by Mike Kryjak 03/06/2023 and supersede above deprecated fits
  /// Data generated using radas https://github.com/cfs-energy/radas
  /// Chose N = 1E20m-3 and tau = 0.5ms based on Moulton, 2021 (DOI: 10.1088/1741-4326/abe4b2)
  /// Those values are applicable in the ITER scrape-off layer but may not be valid in other conditions.
  /// The fits are 10 coefficient polynomials fitted in log-log space like the AMJUEL database in EIRENE

  /// Argon
  struct Argon_adas{
    BoutReal curve(BoutReal Te) {
      BoutReal logT = log(Te);
      BoutReal log_out = 0;

      if (Te >= 1.5 and Te <= 1500) {
        log_out = log_out
        -8.45410692e+01 * pow(logT, 0)
        +1.57727040e+01 * pow(logT, 1)
        -1.54264860e+01 * pow(logT, 2)
        +1.49409902e+01 * pow(logT, 3)
        -1.04815113e+01 * pow(logT, 4)
        +5.00924595e+00 * pow(logT, 5)
        -1.60029106e+00 * pow(logT, 6)
        +3.29455609e-01 * pow(logT, 7)
        -4.14036827e-02 * pow(logT, 8)
        +2.87063206e-03 * pow(logT, 9)
        -8.38888002e-05 * pow(logT, 10);
        return exp(log_out);

    } else if (Te < 1.5) {
        return 1.95353412e-35;   
    } else {
        return 1.22649600e-32;
    }
    }
  };

  /// Neon
  struct Neon_adas{
    BoutReal curve(BoutReal Te) {
      BoutReal logT = log(Te);
      BoutReal log_out = 0;

      if (Te >= 2 and Te <= 1000) {
        log_out = log_out
        -8.21475117e+01 * pow(logT, 0)
        +1.28929854e+01 * pow(logT, 1)
        -4.74266289e+01 * pow(logT, 2)
        +7.45222324e+01 * pow(logT, 3)
        -5.75710722e+01 * pow(logT, 4)
        +2.57375965e+01 * pow(logT, 5)
        -7.12758563e+00 * pow(logT, 6)
        +1.24287546e+00 * pow(logT, 7)
        -1.32943407e-01 * pow(logT, 8)
        +7.97368445e-03 * pow(logT, 9)
        -2.05487897e-04 * pow(logT, 10);
        return exp(log_out);

    } else if (Te < 2) {
        return 6.35304113e-36;   
    } else {
        return 1.17894628e-32;
    }
    }
  };


  /// Nitrogen
  struct Nitrogen_adas{
    BoutReal curve(BoutReal Te) {
      BoutReal logT = log(Te);
      BoutReal log_out = 0;

      if (Te >= 2 and Te <= 500) {
        log_out = log_out
        -5.01649969e+01 * pow(logT, 0)
        -1.35749724e+02 * pow(logT, 1)
        +2.73509608e+02 * pow(logT, 2)
        -2.92109992e+02 * pow(logT, 3)
        +1.90120639e+02 * pow(logT, 4)
        -7.95164871e+01 * pow(logT, 5)
        +2.17762218e+01 * pow(logT, 6)
        -3.88334992e+00 * pow(logT, 7)
        +4.34730098e-01 * pow(logT, 8)
        -2.77683605e-02 * pow(logT, 9)
        +7.72720422e-04 * pow(logT, 10);
        return exp(log_out);

    } else if (Te < 2) {
        return 4.34835380e-34;   
    } else {
        return 8.11096182e-33;
    }
    }
  };


  /// Carbon
  struct Carbon_adas{
    BoutReal curve(BoutReal Te) {
      BoutReal logT = log(Te);
      BoutReal log_out = 0;

      if (Te >= 1 and Te <= 500) {
        log_out = log_out
        -7.87837896e+01 * pow(logT, 0)
        +1.55326376e+00 * pow(logT, 1)
        +1.65898194e+01 * pow(logT, 2)
        -3.23804546e+01 * pow(logT, 3)
        +3.12784663e+01 * pow(logT, 4)
        -1.74826039e+01 * pow(logT, 5)
        +5.91393245e+00 * pow(logT, 6)
        -1.22974105e+00 * pow(logT, 7)
        +1.54004499e-01 * pow(logT, 8)
        -1.06797106e-02 * pow(logT, 9)
        +3.15657594e-04 * pow(logT, 10);
        return exp(log_out);

    } else if (Te < 1) {
        return 6.00623928e-35;   
    } else {
        return 4.53057707e-33;
    }
    }
  };

  /// Argon simplified 1
  /// Based on the ADAS curve above but simplified as a linear interpolation
  /// between the LHS minimum, peak, bottom of RHS slope and final value at Te = 3000eV.
  /// RHS shoulder is preserved.
  /// Helpful for studying impact of cooling curve nonlinearity 
  struct Argon_simplified1{
    BoutReal curve(BoutReal Te) {
      
     if (Te < 0.52) { 
        return 0; 
    
     } else if (Te >= 0.52 and Te < 2.50) {
        return (1.953534e-35*(2.50 - Te) + 6.053680e-34*(Te - 0.52)) / (2.50 - 0.52);

     } else if (Te >= 2.50 and Te < 19.72) {
        return (6.053680e-34*(19.72 - Te) + 2.175237e-31*(Te - 2.50)) / (19.72 - 2.50);

     } else if (Te >= 19.72 and Te < 59.98) {
        return (2.175237e-31*(59.98 - Te) + 4.190918e-32*(Te - 19.72)) / (59.98 - 19.72);

     } else if (Te >= 59.98 and Te < 3000.00) {
        return (4.190918e-32*(3000.00 - Te) + 1.226496e-32*(Te - 59.98)) / (3000.00 - 59.98);

     } else {
        return 1.226496e-32;
     }
    }
  };

  /// Argon simplified 2
  /// Based on the ADAS curve above but simplified as a linear interpolation
  /// between the LHS minimum, peak, and the bottom of RHS slope.
  /// RHS shoulder is eliminated, radiation becomes 0 at 60eV.
  /// Helpful for studying impact of cooling curve nonlinearity 
  struct Argon_simplified2{
    BoutReal curve(BoutReal Te) {
      
     if (Te < 0.52) { 
        return 0; 
    
     } else if (Te >= 0.52 and Te < 2.50) {
        return (1.953534e-35*(2.50 - Te) + 6.053680e-34*(Te - 0.52)) / (2.50 - 0.52);

     } else if (Te >= 2.50 and Te < 19.72) {
        return (6.053680e-34*(19.72 - Te) + 2.175237e-31*(Te - 2.50)) / (19.72 - 2.50);

     } else if (Te >= 19.72 and Te < 59.98) {
        return (2.175237e-31*(59.98 - Te) + 0.000000e+00*(Te - 19.72)) / (59.98 - 19.72);

     } else {
        return 0.0;
     }
    }
  };

  /// Argon simplified 3
  /// Based on the ADAS curve above but simplified as a linear interpolation
  /// between the LHS minimum, peak, and the bottom of RHS slope.
  /// RHS shoulder is eliminated, radiation becomes 0 at 60eV.
  /// LHS / RHS asymmetry is eliminated.
  /// Helpful for studying impact of cooling curve nonlinearity 
  struct Argon_simplified3{
    BoutReal curve(BoutReal Te) {
      
     if (Te < 0.52) { 
        return 0; 
    
     } else if (Te >= 0.52 and Te < 2.50) {
        return (1.953534e-35*(2.50 - Te) + 6.053680e-34*(Te - 0.52)) / (2.50 - 0.52);

     } else if (Te >= 2.50 and Te < 19.72) {
        return (6.053680e-34*(19.72 - Te) + 2.175237e-31*(Te - 2.50)) / (19.72 - 2.50);

     } else if (Te >= 19.72 and Te < 38.02) {
        return (2.175237e-31*(38.02 - Te) + 0.000000e+00*(Te - 19.72)) / (38.02 - 19.72);

     } else {
        return 0.0;
     };
    }
  };

  /// Krypton
  ///
  /// Radas version d50d6c3b (Oct 27, 2023)
  /// Using N = 1E20m-3 and tau = 0.5ms
  struct Krypton_adas {
    BoutReal curve(BoutReal Te) {
      if (Te >= 2 and Te <= 1500) {
        BoutReal logT = log(Te);
        BoutReal log_out =
        -2.97405917e+01 * pow(logT, 0)
        -2.25443986e+02 * pow(logT, 1)
        +4.08713640e+02 * pow(logT, 2)
        -3.86540549e+02 * pow(logT, 3)
        +2.19566710e+02 * pow(logT, 4)
        -7.93264990e+01 * pow(logT, 5)
        +1.86185949e+01 * pow(logT, 6)
        -2.82487288e+00 * pow(logT, 7)
        +2.67070863e-01 * pow(logT, 8)
        -1.43001273e-02 * pow(logT, 9)
        +3.31179737e-04 * pow(logT, 10);
        return exp(log_out);

      } else if (Te < 2) {
        return 8.01651285e-35;
      } else {
        return 6.17035971e-32;
      }
    }
  };

  /// Xenon
  ///
  /// Radas version d50d6c3b (Oct 27, 2023)
  /// Using N = 1E20m-3 and tau = 0.5ms
  ///
  /// Note: Requires more than 10 coefficients to capture multiple
  /// radiation peaks.
  struct Xenon_adas {
    BoutReal curve(BoutReal Te) {
      if (Te >= 2 and Te <= 1300) {
        BoutReal logT = log(Te);
        BoutReal log_out =
        +3.80572137e+02 * pow(logT, 0)
        -3.05839745e+03 * pow(logT, 1)
        +9.01104594e+03 * pow(logT, 2)
        -1.55327244e+04 * pow(logT, 3)
        +1.75819719e+04 * pow(logT, 4)
        -1.38754866e+04 * pow(logT, 5)
        +7.90608220e+03 * pow(logT, 6)
        -3.32041900e+03 * pow(logT, 7)
        +1.03912363e+03 * pow(logT, 8)
        -2.42980719e+02 * pow(logT, 9)
        +4.22211209e+01 * pow(logT, 10)
        -5.36813849e+00 * pow(logT, 11)
        +4.84652106e-01 * pow(logT, 12)
        -2.94023979e-02 * pow(logT, 13)
        +1.07416308e-03 * pow(logT, 14)
        -1.78510623e-05 * pow(logT, 15);
        return exp(log_out);

      } else if (Te < 2) {
        return 1.87187135e-33;
      } else {
        return 1.29785519e-31;
      }
    }
  };

  /// Tungsten
  ///
  /// Radas version d50d6c3b (Oct 27, 2023)
  /// Using N = 1E20m-3 and tau = 0.5ms
  struct Tungsten_adas {
    BoutReal curve(BoutReal Te) {
      if (Te >= 1.25 and Te <= 1500) {
        BoutReal logT = log(Te);
        BoutReal log_out =
          -7.24602210e+01 * pow(logT, 0)
          -2.17524363e+01 * pow(logT, 1)
          +1.90745408e+02 * pow(logT, 2)
          -7.57571067e+02 * pow(logT, 3)
          +1.84119395e+03 * pow(logT, 4)
          -2.99842204e+03 * pow(logT, 5)
          +3.40395125e+03 * pow(logT, 6)
          -2.76328977e+03 * pow(logT, 7)
          +1.63368844e+03 * pow(logT, 8)
          -7.11076320e+02 * pow(logT, 9)
          +2.28027010e+02 * pow(logT, 10)
          -5.30145974e+01 * pow(logT, 11)
          +8.46066686e+00 * pow(logT, 12)
          -7.54960450e-01 * pow(logT, 13)
          -1.61054010e-02 * pow(logT, 14)
          +1.65520152e-02 * pow(logT, 15)
          -2.70054697e-03 * pow(logT, 16)
          +2.50286873e-04 * pow(logT, 17)
          -1.43319310e-05 * pow(logT, 18)
          +4.75258630e-07 * pow(logT, 19)
          -7.03012454e-09 * pow(logT, 20);
        return exp(log_out);

      } else if (Te < 1.25) {
        return 2.09814651e-32;
      } else {
        return 1.85078767e-31;
      }
    }
  };
}

/// Set ion densities from electron densities
///
template <typename CoolingCurve>
struct FixedFractionRadiation : public Component {
  /// Inputs
  /// - <name>
  ///   - fraction
  FixedFractionRadiation(std::string name, Options &alloptions, Solver *UNUSED(solver)) : name(name) {
    auto& options = alloptions[name];

    fraction = options["fraction"]
      .doc("Impurity ion density as fraction of electron density")
      .withDefault(0.0);

    diagnose = options["diagnose"]
      .doc("Output radiation diagnostic?")
      .withDefault<bool>(false);

    radiation_multiplier = options["R_multiplier"]
      .doc("Scale the radiation rate by this factor")
      .withDefault<BoutReal>(1.0);

    // Get the units
    auto& units = alloptions["units"];
    Tnorm = get<BoutReal>(units["eV"]);
    Nnorm = get<BoutReal>(units["inv_meters_cubed"]);
    FreqNorm = 1. / get<BoutReal>(units["seconds"]);
  }

  /// Required inputs
  ///
  /// - species
  ///   - e
  ///     - density
  ///     - temperature
  ///
  /// Sets the electron energy loss
  ///
  /// - species
  ///   - e
  ///     - energy_source
  ///
  void transform(Options &state) override {
    auto& electrons = state["species"]["e"];
    // Don't need boundary cells
    const Field3D Ne = GET_NOBOUNDARY(Field3D, electrons["density"]);
    const Field3D Te = GET_NOBOUNDARY(Field3D, electrons["temperature"]);
    
    radiation = cellAverage(
                            [&](BoutReal ne, BoutReal te) {
                              if (ne < 0.0 or te < 0.0) {
                                return 0.0;
                              }
                              // Fixed fraction ions
                              const BoutReal ni = fraction * ne;
                              // cooling in Wm^3 so normalise.
                              // Note factor of qe due to Watts rather than eV
                              return ne * ni * cooling.curve(te * Tnorm) * radiation_multiplier * 
                              Nnorm / (SI::qe * Tnorm * FreqNorm);
                            },
                            Ne.getRegion("RGN_NOBNDRY"))(Ne, Te);

    // Remove radiation from the electron energy source
    subtract(electrons["energy_source"], radiation);
  }

  void outputVars(Options& state) override {
    AUTO_TRACE();

    if (diagnose) {
      set_with_attrs(state[std::string("R") + name], -radiation,
                     {{"time_dimension", "t"},
                      {"units", "W / m^3"},
                      {"conversion", SI::qe * Tnorm * Nnorm * FreqNorm},
                      {"long_name", std::string("Radiation cooling ") + name},
                      {"source", "fixed_fraction_radiation"}});
    }
  }
 private:
  std::string name;

  CoolingCurve cooling; ///< The cooling curve L(T) -> Wm^3
  BoutReal fraction; ///< Fixed fraction

  bool diagnose; ///< Output radiation diagnostic?
  BoutReal radiation_multiplier; ///< Scale the radiation rate by this factor
  Field3D radiation; ///< For output diagnostic

  // Normalisations
  BoutReal Tnorm, Nnorm, FreqNorm;
};

namespace {
  RegisterComponent<FixedFractionRadiation<HutchinsonCarbon>>
    registercomponentfixedfractionhutchinsonbcarbon("fixed_fraction_hutchinson_carbon");
    
  RegisterComponent<FixedFractionRadiation<Carbon_adas>>
    registercomponentfixedfractioncarbon("fixed_fraction_carbon");  

  RegisterComponent<FixedFractionRadiation<Nitrogen_adas>>
    registercomponentfixedfractionnitrogen("fixed_fraction_nitrogen");

  RegisterComponent<FixedFractionRadiation<Neon_adas>>
    registercomponentfixedfractionneon("fixed_fraction_neon");

  RegisterComponent<FixedFractionRadiation<Argon_adas>>
    registercomponentfixedfractionargon("fixed_fraction_argon");

  RegisterComponent<FixedFractionRadiation<Argon_simplified1>>
    registercomponentfixedfractionargonsimplified1("fixed_fraction_argon_simplified1");

  RegisterComponent<FixedFractionRadiation<Argon_simplified2>>
    registercomponentfixedfractionargonsimplified2("fixed_fraction_argon_simplified2");

  RegisterComponent<FixedFractionRadiation<Argon_simplified3>>
    registercomponentfixedfractionargonsimplified3("fixed_fraction_argon_simplified3");

  RegisterComponent<FixedFractionRadiation<Krypton_adas>>
    registercomponentfixedfractionkrypton("fixed_fraction_krypton");

  RegisterComponent<FixedFractionRadiation<Xenon_adas>>
    registercomponentfixedfractionxenon("fixed_fraction_xenon");

  RegisterComponent<FixedFractionRadiation<Tungsten_adas>>
    registercomponentfixedfractiontungsten("fixed_fraction_tungsten");
}

#endif // FIXED_FRACTION_IONS_H
