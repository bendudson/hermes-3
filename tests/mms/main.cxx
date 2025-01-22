#include "bout/field_factory.hxx"

//#include "hermes-2.hxx"
#include "bout/version.hxx"
#include "div_ops.hxx"
#include "bout/fv_ops.hxx"
#include "bout/difops.hxx"
#include "../include/div_ops.hxx"

// Class used to store function in operator list below
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
  
  std::string differential_operator_name = Options::root()["mesh"]["differential_operator_name"].withDefault("FV::Div_a_Grad_perp(a, f)");
  
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

  // Get expected result from input file
  Field3D expected_result{mesh};
  mesh->get(expected_result, "expected_result", 0.0, false);
  dump["expected_result"] = expected_result;

  for (const auto& difop: differential_operators) {
	    if ( difop.name.compare(differential_operator_name) == 0){
          Field3D result = difop.func(a, f);
          dump["result"] = result;
          dump["result"].setAttributes({
                {"operator", difop.name},
            });
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
