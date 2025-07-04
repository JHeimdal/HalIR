//#include <cstring>
#include <vector>
#include <string>
#include <boost/json/src.hpp>

extern "C" {
#include "HalIR/halir.h"
}

namespace json = boost::json;

int main(int argc, char **argv)
{
  enum TEST {
    TEMP_F_TO_K,
  };

  using namespace std;

  vector<string> args(argv + 1, argv + argc);
  TEST test = (TEST)stoi(args[0]);

  json::value inp1 = {
    { "input", {
      {"project", {
        {"pname", "CO_Test2"},
        {"rootDir", "/home/jimmy/Programs/HalIR/tests/TestCalc/"},
        {"hapi_db", ""},
        {"pcomments", "CO Test on some different concentrations"},
        {"pfiles", {nullptr}}}},
      {"sampleEnv", {
        {"temp", 303.15},
        {"tempU", "K"},
        {"press", 1.0},
        {"pressU", "atm"},
        {"pathL", 300.0},
        {"pathLU", "m"},
        {"ROI", {2000.0, 2245.0}},
        {"res", 0.1},
        {"apod", "Boxcar"},
        {"fov", 0.1},
        {"ftype", "Transmission"},
        {"bgfile", ""}}},
      {"composition", { {
        {"molec", "CO"},
        {"isotop", "Natural"},
        {"vmr", 5.921539600296e-05},
        {"prmfile", ""} }, {
        {"molec", "H2O"},
        {"isotop", "Natural"},
        {"vmr", 1.92e-03},
        {"prmfile", ""}}}}
      }
    }
  };

  halir_workspace *work_p1;
  std::string ref1_str = json::serialize(inp1);
  halir_workspace *work_ref = halir_parseJSONinput(ref1_str.c_str());
  string err_msg;

  switch (test) {
    case TEMP_F_TO_K:
      inp1.at_pointer("/input/sampleEnv/temp") = 86;
      inp1.at_pointer("/input/sampleEnv/tempU") = "F";
      err_msg = "Error in temperature conversion F->K\n";
      break;
  }
  return 0;
}
