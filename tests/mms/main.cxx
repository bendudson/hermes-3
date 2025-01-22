#include "bout/field_factory.hxx"

//#include "hermes-2.hxx"
#include "bout/version.hxx"
#include "div_ops.hxx"
#include "bout/fv_ops.hxx"
#include "bout/difops.hxx"
#include "../include/div_ops.hxx"

// Class used to store function of two arguments
// in operator list below
class nameandfunction2 {
public:
  std::string name;
  std::function<Field3D(const Field3D&, const Field3D&)> func;
};

// List of tested operators
const auto differential_operators = {
    nameandfunction2{"FV::Div_a_Grad_perp(a, f)",
                 [](const Field3D &a, const Field3D &f) {
                   return FV::Div_a_Grad_perp(a, f);
                 }},
    nameandfunction2{"Div_a_Grad_perp_nonorthog(a, f)",
                 [](const Field3D &a, const Field3D &f) {
                   return Div_a_Grad_perp_nonorthog(a, f);
                 }},
};

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);
  
  int n_operators = Options::root()["mesh"]["n_operators"].withDefault(1);
  
  Mesh* mesh = Mesh::create(&Options::root()["mesh"]);
  mesh->load();

  Options dump;

  // Note that sticking all the inputs in the [mesh] section is a bit of a hack, but makes
  // for a more compact input file than using the 'standard' variable-initialisation
  // routines, which would require a separate section for each variable.

  // Load and save coordinate variables
  Field3D xx{mesh}, yy{mesh}, zz{mesh};
  mesh->get(xx, "x_input", 0.0, false);
  mesh->get(yy, "y_input", 0.0, false);
  mesh->get(zz, "z_input", 0.0, false);
  dump["x_input"] = xx;
  dump["y_input"] = yy;
  dump["z_input"] = zz;

  // Get coefficient from input file
  Field3D a{mesh};
  mesh->get(a, "a", 1.0, false);
  mesh->communicate(a);
  dump["a"] = a;

  // Get test variable from input file
  Field3D f{mesh};
  mesh->get(f, "f", 0.0, false);
  mesh->communicate(f);
  dump["f"] = f;

  for (int i = 0; i < n_operators; i++){
      std::string inputname = "differential_operator_name_"+std::to_string(i);
      std::string expectedname = "expected_result_"+std::to_string(i);
      std::string outname = "result_"+std::to_string(i);
      std::string differential_operator_name = Options::root()["mesh"][inputname].withDefault("FV::Div_a_Grad_perp(a, f)");
      // the for loop and if statement below should be replaced
      // by a neater indexing syntax below if possible
      for (const auto& difop: differential_operators) {
          if (difop.name.compare(differential_operator_name) == 0){
              // Get result of applying the named differential operator
              Field3D result = difop.func(a, f);
              dump[outname] = result;
              dump[outname].setAttributes({
                    {"operator", difop.name},
                });
              // Get expected result from input file
              Field3D expected_result{mesh};
              mesh->get(expected_result, expectedname, 0.0, false);
              dump[expectedname] = expected_result;
          }
      }
  }
  //Field3D result = FV::Div_a_Grad_perp(a, f);
  //dump["result"] = result;

  //Field3D result_nonorthog = Div_a_Grad_perp_nonorthog(a, f);
  //dump["result_nonorthog"] = result_nonorthog;

  mesh->outputVars(dump);

  std::string outname = fmt::format(
      "{}/BOUT.{}.nc",
      Options::root()["datadir"].withDefault<std::string>("data"), BoutComm::rank());
  
  bout::OptionsIO::create(outname)->write(dump);

  BoutFinalise();

  return 0;
}
