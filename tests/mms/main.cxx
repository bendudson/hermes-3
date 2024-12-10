#include "bout/field_factory.hxx"

//#include "hermes-2.hxx"
#include "bout/version.hxx"
#include "div_ops.hxx"
#include "bout/fv_ops.hxx"
#include "bout/difops.hxx"

std::vector<std::string> getAll(std::string str) {
  std::vector<std::string> out{};
  int i = 0;
  while (true) {
    std::string section = fmt::format(str, i++);

    output.write("Trying '{}'\n", section);
    //auto sec = Options::root()[section];
    //if (sec.getChildren().empty()){
    if (not Options::root().isSection(section)) {
      return out;
    }
    out.push_back(section);
  }
}

class nameandfunc2 {
public:
  std::string name;
  std::function<Field3D(const Field3D&, const Field3D&)> func;
};

class nameandfunc1 {
public:
  std::string name;
  std::function<Field3D(const Field3D&)> func;
};

const auto functions = {
    nameandfunc2{"r", [](const Field3D &R,
			 const Field3D &Z) { return sqrt(Z*Z + (R-5)*(R-5)); }},
    nameandfunc2{"sin(theta)", [](const Field3D &R,
				  const Field3D &Z) { return Z / sqrt(Z*Z + (R-5)*(R-5)); }},
    nameandfunc2{"R", [](const Field3D &R, const Field3D &Z) { return R; }},
    nameandfunc2{"R²",
                 [](const Field3D &R, const Field3D &Z) { return R * R; }},
    nameandfunc2{"sin(R)",
                 [](const Field3D &R, const Field3D &Z) { return sin(R); }},
    nameandfunc2{"sin(Z)",
                 [](const Field3D &R, const Field3D &Z) { return sin(Z); }},
    nameandfunc2{
        "sin(Z)*sin(R)",
        [](const Field3D &R, const Field3D &Z) { return sin(Z) * sin(R); }},
    nameandfunc2{"sin(10*R)", [](const Field3D &R,
                                 const Field3D &Z) { return sin(10 * R); }},
    nameandfunc2{"sin(100*R)", [](const Field3D &R,
                                  const Field3D &Z) { return sin(100 * R); }},
    nameandfunc2{"sin(1000*R)", [](const Field3D &R,
                                   const Field3D &Z) { return sin(1000 * R); }},
};

const auto difops = {
    nameandfunc2{"FV::Div_a_Grad_perp(1, f)",
                 [](const Field3D &a, const Field3D &f) {
                   return FV::Div_a_Grad_perp(a, f);
                 }},
    nameandfunc2{"Delp2(f)",
                 [](const Field3D &a, const Field3D &f) { return Delp2(f); }},
    nameandfunc2{"Laplace(f)",
                 [](const Field3D &a, const Field3D &f) { return Laplace(f); }},
    // nameandfunc2{"newDelp2(f)", [] (const Field3D& a, const Field3D& f) {
    // return newDelp2.apply(f); }}, nameandfunc2("bracket(a, f)", [] (const
    // Field3D& a, const Field3D& f) { return bracket(a, f); }},
};

#include <list>
#include <tuple>

const std::list<std::tuple<nameandfunc2, nameandfunc2>> functions2 = {
    // const auto functions2 = {
    {{"R", [](const Field3D &R, const Field3D &Z) { return R; }},
     {"Z", [](const Field3D &R, const Field3D &Z) { return Z; }}},
    {{"R", [](const Field3D &R, const Field3D &Z) { return R; }},
     {"R", [](const Field3D &R, const Field3D &Z) { return R; }}},
    {{"Z", [](const Field3D &R, const Field3D &Z) { return Z; }},
     {"Z", [](const Field3D &R, const Field3D &Z) { return Z; }}},
    {{"Z", [](const Field3D &R, const Field3D &Z) { return Z; }},
     {"R", [](const Field3D &R, const Field3D &Z) { return R; }}},
    {{"sin(R)", [](const Field3D &R, const Field3D &Z) { return sin(R); }},
     {"sin(Z)", [](const Field3D &R, const Field3D &Z) { return sin(Z); }}},
    {{"sin(R*10)",
      [](const Field3D &R, const Field3D &Z) { return sin(R * 10); }},
     {"sin(Z*10)",
      [](const Field3D &R, const Field3D &Z) { return sin(Z * 10); }}},
    {{"sin(R*100)",
      [](const Field3D &R, const Field3D &Z) { return sin(R * 100); }},
     {"sin(Z*100)",
      [](const Field3D &R, const Field3D &Z) { return sin(Z * 100); }}},
    // {{"sin(R*1000)", [] (const Field3D& R, const Field3D& Z) {return
    // sin(R*1000); }}, {"sin(Z*1000)", [] (const Field3D& R, const Field3D& Z)
    // {return sin(Z*1000); }}},
};

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);

  auto meshes = getAll("mesh_{}");
  auto fields = getAll("field_{}");
  
  output.write("Found {:d} meshes and {:d} fields", meshes.size(), fields.size());
  for (const auto& meshname: meshes) {
    Mesh* mesh = Mesh::create(&Options::root()[meshname]);
    mesh->load();

    // FCI::dagp dagp(*mesh);
    // FCI::dagp_fv dagp(*mesh);

    Options dump;
    Field3D R{mesh}, Z{mesh};
    mesh->get(R, "Rxy", 0.0, false);
    mesh->get(Z, "Zxy", 0.0, false);

    // dump["R_3d"] = R;
    // dump["Z_3d"] = Z;


    auto coord = mesh->getCoordinates();
    //coord->g23, coord->g_23, coord->dy, coord->dz, coord->Bxy, coord->J)

    
    Field3D a{1.0, mesh};
    
    mesh->communicate(a, coord->g12, coord->g_12, coord->g23, coord->g_23, coord->dy, coord->dz, coord->Bxy, coord->J);
    a.applyParallelBoundary("parallel_neumann_o2");
    coord->g23.applyParallelBoundary("parallel_neumann_o2");
    coord->g_23.applyParallelBoundary("parallel_neumann_o2");
    coord->g12.applyParallelBoundary("parallel_neumann_o2");
    coord->g_12.applyParallelBoundary("parallel_neumann_o2");
    coord->dy.applyParallelBoundary("parallel_neumann_o2");
    coord->dz.applyParallelBoundary("parallel_neumann_o2");
    coord->Bxy.applyParallelBoundary("parallel_neumann_o2");
    coord->J.applyParallelBoundary("parallel_neumann_o2");
    
    int i = 0;
    for (const auto& func: functions) {
      auto f = func.func(R, Z);
      mesh->communicate(f);
      f.applyParallelBoundary("parallel_neumann_o2");
      for (const auto& dif: difops) {
	auto outname = fmt::format("out_{}", i++);
	dump[outname] = dif.func(a, f);
        dump[outname].setAttributes({
            {"operator", dif.name},
            {"function", func.name},
            {"f", func.name},
            {"inp", func.name},
        });
      }
     // ////
     // {
     //   const auto outname = fmt::format("out_{}", i++);
     //   dump[outname] = dagp(a, f);
     //   const auto opname = "FCI::dagp_fv(1, f)";
     //   dump[outname].setAttributes({
     //       {"operator", opname},
     //       {"f", func.name},
     //       {"inp", func.name},
     //   });
     // }
     // ////
    }

    for (const auto &func : functions2) {
      auto a = std::get<0>(func).func(R, Z);
      auto f = std::get<1>(func).func(R, Z);
      mesh->communicate(a, f);
      a.applyParallelBoundary("parallel_neumann_o2");
      f.applyParallelBoundary("parallel_neumann_o2");
      {
        const auto outname = fmt::format("out_{}", i++);
        dump[outname] =
            sqrt(coord->g_22) / coord->J * bracket(a, f, BRACKET_ARAKAWA);
        const auto opname = "bracket(a, f)";
        dump[outname].setAttributes({
            {"operator", opname},
            {"a", std::get<0>(func).name},
            {"f", std::get<1>(func).name},
            {"inp", fmt::format("{}, {}", std::get<0>(func).name,
                                std::get<1>(func).name)},
        });
      }
      {
        const auto outname = fmt::format("out_{}", i++);
        dump[outname] =
            sqrt(coord->g_22) / coord->J * bracket(a, f, BRACKET_ARAKAWA_OLD);
        const auto opname = "bracket(a, f, OLD)";
        dump[outname].setAttributes({
            {"operator", opname},
            {"a", std::get<0>(func).name},
            {"f", std::get<1>(func).name},
            {"inp", fmt::format("{}, {}", std::get<0>(func).name,
                                std::get<1>(func).name)},
        });
      }
     // {
     //   const auto outname = fmt::format("out_{}", i++);
     //   dump[outname] = FCI::Div_a_Grad_perp(a, f);
     //   const auto opname = "FCI::Div_a_Grad_perp(a, f)";
     //   dump[outname].setAttributes({
     //       {"operator", opname},
     //       {"a", std::get<0>(func).name},
     //       {"f", std::get<1>(func).name},
     //       {"inp", fmt::format("{}, {}", std::get<0>(func).name,
     //                           std::get<1>(func).name)},
     //   });
     // }
     // {
     //   const auto outname = fmt::format("out_{}", i++);
     //   dump[outname] = dagp(a, f);
     //   const auto opname = "FCI::dagp_fv(f)";
     //   dump[outname].setAttributes({
     //       {"operator", opname},
     //       {"a", std::get<0>(func).name},
     //       {"f", std::get<1>(func).name},
     //       {"inp", fmt::format("{}, {}", std::get<0>(func).name,
     //                           std::get<1>(func).name)},
     //   });
     // }
    }
    if (mesh) {
      mesh->outputVars(dump);
      dump["BOUT_VERSION"].force(bout::version::as_double);
    }

      std::string outname = fmt::format(
          "{}/BOUT.{}.{}.nc",
          Options::root()["datadir"].withDefault<std::string>("data"), meshname, BoutComm::rank());
      
      bout::OptionsIO::create(outname)->write(dump);
      
  };
  
  BoutFinalise()    ;
}
//   std::vector<Field3D> fields;
//   fields.resize(static_cast<int>(BoundaryParType::SIZE));
//   Options dump;
//   for (int i=0; i< fields.size(); i++){
//     fields[i] = Field3D{0.0};
//     mesh->communicate(fields[i]);
//     for (const auto &bndry_par : mesh->getBoundariesPar(static_cast<BoundaryParType>(i))) {
//       output.write("{:s} region\n", toString(static_cast<BoundaryParType>(i)));
//       for (bndry_par->first(); !bndry_par->isDone(); bndry_par->next()) {
//         fields[i][bndry_par->ind()] += 1;
//         output.write("{:s} increment\n", toString(static_cast<BoundaryParType>(i)));
//       }
//     }
//     output.write("{:s} done\n", toString(static_cast<BoundaryParType>(i)));
//     dump[fmt::format("field_{:s}", toString(static_cast<BoundaryParType>(i)))] = fields[i];
//   }

//   bout::writeDefaultOutputFile(dump);

//   BoutFinalise();
// }




