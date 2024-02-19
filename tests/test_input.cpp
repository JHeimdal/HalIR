#include <boost/json/serialize.hpp>
#include <boost/json/value.hpp>
#include <cstring>
#include <iostream>
#include <boost/json/src.hpp>
#include <limits>
#include <string>

extern "C" {
#include "HalIR/halir.h"
}

namespace json = boost::json;

bool compare_doubles(double *A, double *B) {
  double EPSILON = std::numeric_limits<double>::epsilon();
  if (std::abs(*A-*B) <= EPSILON * std::max(1.0, *A))
    return true;
  else
    return false;
}

bool compare_halir_workspace(halir_workspace *ref, halir_workspace *test) {
  if (strncmp(ref->pname, test->pname, PATH_LEN) != 0)
    return false;
  else if (strncmp(ref->rootDir, test->rootDir, PATH_LEN) != 0)
    return false;
  else if (strncmp(ref->pcomments, test->pcomments, COMMENT_LEN) != 0)
    return false;
  else if (strncmp(ref->bgfile, test->bgfile, PATH_LEN) != 0)
    return false;
  else if (ref->tempU != test->tempU)
    return false;
  else if (ref->pressU != test->pressU)
    return false;
  else if (ref->pathLU != test->pathLU)
    return false;
  else if (ref->apod != test->apod)
    return false;
  else if (ref->ftype != test->ftype)
    return false;
  //char **pfiles;
  //size_t pfiles_length;
  else if (!compare_doubles(&ref->temp, &test->temp))
    return false;
  else if (!compare_doubles(&ref->press, &test->press))
    return false;
  else if (!compare_doubles(&ref->pathL, &test->pathL))
    return false;
  else if (!compare_doubles(&ref->res, &test->res))
    return false;
  else if (!compare_doubles(&ref->fov, &test->fov))
    return false;
  else if (!compare_doubles(&ref->ROI[0], &test->ROI[0]))
    return false;
  else if (!compare_doubles(&ref->ROI[1], &test->ROI[1]))
    return false;
  //// Composition of system
  //size_t composition_length;
  //halir_compound **composition;
  return true;
}

int main()
{
  json::value ref1 = {
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

  inp1.at_pointer("/input/sampleEnv/temp") = 86;
  inp1.at_pointer("/input/sampleEnv/tempU") = "F";
  std::string inp2_str = json::serialize(inp1);
  work_p1 = halir_parseJSONinput(inp2_str.c_str());
  if (!compare_halir_workspace(work_ref, work_p1)) {
    std::cout << "Error in temperature conversion F->K\n";
    return 1;
  }
  inp1.at_pointer("/input/sampleEnv/temp") = 30;
  inp1.at_pointer("/input/sampleEnv/tempU") = "C";
  inp2_str = json::serialize(inp1);
  work_p1 = halir_parseJSONinput(inp2_str.c_str());
  if (!compare_halir_workspace(work_ref, work_p1)) {
    std::cout << "Error in temperature conversion C->K\n";
    return 1;
  }
  //halir_print_workspace(work_ref);
  return 0;
}
