#include <boost/json/serialize.hpp>
#include <boost/json/value.hpp>
#include <cstddef>
#include <cstring>
#include <ios>
#include <iostream>
#include <boost/json/src.hpp>
#include <limits>
#include <string>

extern "C" {
#include "HalIR/halir.h"
}

namespace json = boost::json;

bool compare_halir_num(halir_num *A, halir_num *B) {
  //double EPSILON = std::numeric_limits<double>::epsilon();
  halir_num EPSILON = 1e-9;
  if (std::abs(*A-*B) <= EPSILON * std::max(1.0f, *A))
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
  else if (!compare_halir_num(&ref->temp, &test->temp))
    return false;
  else if (!compare_halir_num(&ref->press, &test->press))
    return false;
  else if (!compare_halir_num(&ref->pathL, &test->pathL))
    return false;
  else if (!compare_halir_num(&ref->res, &test->res))
    return false;
  else if (!compare_halir_num(&ref->fov, &test->fov))
    return false;
  else if (!compare_halir_num(&ref->ROI[0], &test->ROI[0]))
    return false;
  else if (!compare_halir_num(&ref->ROI[1], &test->ROI[1]))
    return false;
  //// Composition of system
  if (ref->composition_length == test->composition_length) {
    // Safe to do comparison.
    for (size_t comp=0; comp < ref->composition_length; comp++) {
      if (!compare_halir_num(&ref->composition[comp]->conc, &test->composition[comp]->conc))
        return false;
      else if (ref->composition[comp]->concU != test->composition[comp]->concU)
        return false;
      else if (!compare_halir_num(&ref->composition[comp]->vmr, &test->composition[comp]->vmr))
        return false;
    }
  } else {
    return false;
  }
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
  inp1.at_pointer("/input/sampleEnv/press") = 1.0132500;
  inp1.at_pointer("/input/sampleEnv/pressU") = "bar";
  inp2_str = json::serialize(inp1);
  work_p1 = halir_parseJSONinput(inp2_str.c_str());
  if (!compare_halir_workspace(work_ref, work_p1)) {
    std::cout << "Error in temperature conversion bar->atm\n";
    return 1;
  }
  inp1.at_pointer("/input/sampleEnv/press") = 1.0132500e3;
  inp1.at_pointer("/input/sampleEnv/pressU") = "mbar";
  inp2_str = json::serialize(inp1);
  work_p1 = halir_parseJSONinput(inp2_str.c_str());
  if (!compare_halir_workspace(work_ref, work_p1)) {
    std::cout << "Error in temperature conversion mbar->atm\n";
    return 1;
  }
  inp1.at_pointer("/input/sampleEnv/press") = 1.0132500e5;
  inp1.at_pointer("/input/sampleEnv/pressU") = "pa";
  inp2_str = json::serialize(inp1);
  work_p1 = halir_parseJSONinput(inp2_str.c_str());
  if (!compare_halir_workspace(work_ref, work_p1)) {
    std::cout << "Error in temperature conversion pa->atm\n";
    return 1;
  }
  inp1.at_pointer("/input/sampleEnv/press") = 1.0132500e3;
  inp1.at_pointer("/input/sampleEnv/pressU") = "hpa";
  inp2_str = json::serialize(inp1);
  work_p1 = halir_parseJSONinput(inp2_str.c_str());
  if (!compare_halir_workspace(work_ref, work_p1)) {
    std::cout << "Error in temperature conversion hpa->atm\n";
    return 1;
  }
  inp1.at_pointer("/input/sampleEnv/press") = 760.00000;
  inp1.at_pointer("/input/sampleEnv/pressU") = "mmhg";
  inp2_str = json::serialize(inp1);
  work_p1 = halir_parseJSONinput(inp2_str.c_str());
  if (!compare_halir_workspace(work_ref, work_p1)) {
    std::cout.precision(13);
    std::cout << std::scientific;
    std::cout << "Error in temperature conversion mmHg->atm\n";
    std::cout << "ref: " << work_ref->press << " inp: " << work_p1->press << std::endl;
    return 1;
  }
  //halir_print_workspace(work_ref);
  return 0;
}
