#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

extern "C" {
#include "HalIR/halir.h"
#include "HalIR/spc.h"
}

bool compare_halir_num(halir_num *A, halir_num *B) {
  //double EPSILON = std::numeric_limits<double>::epsilon();
  halir_num EPSILON = 1e-9;
  if (std::abs(*A-*B) <= EPSILON * std::max(1.0, *A))
    return true;
  else
    return false;
}

bool compare_halir_workspace(halir_simulation_setup *ref, halir_simulation_setup *test) {
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

bool compare_halir_spectra_values(const halir_spectra *ref, const halir_spectra *test) {
  if ((ref == nullptr) || (test == nullptr))
    return false;
  if (ref->ndatapnts != test->ndatapnts)
    return false;

  for (size_t i = 0; i < ref->ndatapnts; i++) {
    if (!compare_halir_num(&ref->wavenum[i], &test->wavenum[i]))
      return false;
    if (!compare_halir_num(&ref->data[i], &test->data[i]))
      return false;
  }

  return true;
}

bool compare_halir_result_values(const halir_result *ref, const halir_result *test) {
  if ((ref == nullptr) || (test == nullptr))
    return false;
  if (ref->nspectra != test->nspectra)
    return false;

  for (size_t i = 0; i < ref->nspectra; i++) {
    if (!compare_halir_spectra_values(&ref->spectra[i], &test->spectra[i]))
      return false;
  }

  return true;
}

bool compare_spc_roundtrip_values(const halir_spectra *ref, const halir_spectra *test) {
  if ((ref == nullptr) || (test == nullptr))
    return false;
  if (ref->ndatapnts != test->ndatapnts)
    return false;

  for (size_t i = 0; i < ref->ndatapnts; i++) {
    halir_num expected_x = static_cast<halir_num>(static_cast<float>(ref->wavenum[i]));
    halir_num expected_y = static_cast<halir_num>(static_cast<float>(ref->data[i]));
    if (!compare_halir_num(&expected_x, &test->wavenum[i]))
      return false;
    if (!compare_halir_num(&expected_y, &test->data[i]))
      return false;
  }

  return true;
}

bool compare_spc_roundtrip_result(const halir_result *ref, const halir_result *test) {
  if ((ref == nullptr) || (test == nullptr))
    return false;
  if (ref->nspectra != test->nspectra)
    return false;

  for (size_t i = 0; i < ref->nspectra; i++) {
    if (!compare_spc_roundtrip_values(&ref->spectra[i], &test->spectra[i]))
      return false;
  }

  return true;
}

bool replace_once(std::string &text, const std::string &search_text, const std::string &replacement) {
  const std::size_t pos = text.find(search_text);
  if (pos == std::string::npos)
    return false;

  text.replace(pos, search_text.size(), replacement);
  return true;
}

std::string build_input_json(const std::string &pname,
                             const std::string &temp,
                             const std::string &temp_unit,
                             const std::string &press,
                             const std::string &press_unit,
                             const std::string &roi)
{
  // Keep test fixture generation dependency-free so CI does not need Boost.
  std::string json_text = R"({
    "input": {
      "project": {
        "pname": "__PNAME__",
        "rootDir": "TestCalc",
        "hapi_db": "",
        "pcomments": "CO Test on some different concentrations",
        "pfiles": [null]
      },
      "sampleEnv": {
        "temp": __TEMP__,
        "tempU": "__TEMP_UNIT__",
        "press": __PRESS__,
        "pressU": "__PRESS_UNIT__",
        "pathL": 300.0,
        "pathLU": "m",
        "ROI": __ROI__,
        "res": 0.1,
        "apod": "Boxcar",
        "fov": 0.1,
        "ftype": "Transmission",
        "bgfile": ""
      },
      "composition": [
        {
          "molec": "CO",
          "isotop": "Natural",
          "vmr": 5.921539600296e-05,
          "prmfile": ""
        },
        {
          "molec": "H2O",
          "isotop": "Natural",
          "vmr": 1.92e-03,
          "prmfile": ""
        }
      ]
    }
  })";

  if (!replace_once(json_text, "__PNAME__", pname) ||
      !replace_once(json_text, "__TEMP__", temp) ||
      !replace_once(json_text, "__TEMP_UNIT__", temp_unit) ||
      !replace_once(json_text, "__PRESS__", press) ||
      !replace_once(json_text, "__PRESS_UNIT__", press_unit) ||
      !replace_once(json_text, "__ROI__", roi)) {
    std::cerr << "Failed to populate JSON test fixture template" << std::endl;
    return {};
  }

  return json_text;
}

int main(int argc, char **argv)
{
  enum TEST {
    TEMP_F_TO_K,
    TEMP_C_TO_K,
    PRESS_BAR_TO_ATM,
    PRESS_MBAR_TO_ATM,
    PRESS_PA_TO_ATM,
    PRESS_HPA_TO_ATM,
    PRESS_MMHG_TO_ATM,
    ROI_TOO_MANY_VALUES,
    ROI_NON_NUMERIC,
    PNAME_TOO_LONG,
    WORKSPACE_CREATE_FREE,
    PROJECT_SETTER_API,
    PROJECT_FILE_APPEND_API,
    PROJECT_API_INVALID_ARGS,
    SAMPLE_ENV_SETTER_API,
    SAMPLE_ENV_UNIT_NORMALIZATION,
    SAMPLE_ENV_INVALID_ARGS,
    COMPOSITION_ADD_API,
    COMPOSITION_SET_VMR_API,
    COMPOSITION_SET_CONCENTRATION_API,
    COMPOSITION_API_INVALID_ARGS,
    WORKSPACE_VALIDATION_SUCCESS,
    WORKSPACE_VALIDATION_INVALID_ARGS,
    RESULT_CREATE_FREE_API,
    RESULT_API_INVALID_ARGS,
    SPECTRA_SPC_ROUNDTRIP,
    RESULT_SPC_ROUNDTRIP,
    RESULT_SPC_DIRECTORY_READ,
  };

  using namespace std;

  vector<string> args(argv + 1, argv + argc);
  TEST test = (TEST)stoi(args[0]);

  if (test == WORKSPACE_CREATE_FREE) {
    halir_simulation_setup *work = halir_simulation_setup_create();
    if (work == nullptr) {
      std::cout << "halir_simulation_setup_create returned NULL" << std::endl;
      return 1;
    }
    halir_simulation_setup_free(work);
    return 0;
  }

  if (test == PROJECT_SETTER_API) {
    halir_simulation_setup *work = halir_simulation_setup_create();
    if (work == nullptr) {
      std::cout << "halir_simulation_setup_create returned NULL" << std::endl;
      return 1;
    }

    int rc = halir_simulation_setup_set_project(work, "MyProject", "/tmp", "A comment");
    if (rc != 0) {
      std::cout << "halir_simulation_setup_set_project failed" << std::endl;
      halir_simulation_setup_free(work);
      return 1;
    }

    if (strncmp(work->pname, "MyProject", PATH_LEN) != 0 ||
        strncmp(work->rootDir, "/tmp", PATH_LEN) != 0 ||
        strncmp(work->pcomments, "A comment", COMMENT_LEN) != 0) {
      std::cout << "halir_simulation_setup_set_project did not persist values" << std::endl;
      halir_simulation_setup_free(work);
      return 1;
    }

    halir_simulation_setup_free(work);
    return 0;
  }

  if (test == PROJECT_FILE_APPEND_API) {
    halir_simulation_setup *work = halir_simulation_setup_create();
    if (work == nullptr) {
      std::cout << "halir_simulation_setup_create returned NULL" << std::endl;
      return 1;
    }

    if (halir_simulation_setup_add_project_file(work, "a.txt") != 0 ||
        halir_simulation_setup_add_project_file(work, "b.txt") != 0) {
      std::cout << "halir_simulation_setup_add_project_file failed" << std::endl;
      halir_simulation_setup_free(work);
      return 1;
    }

    if (work->pfiles_length != 2 ||
        strncmp(work->pfiles[0], "a.txt", PATH_LEN) != 0 ||
        strncmp(work->pfiles[1], "b.txt", PATH_LEN) != 0) {
      std::cout << "halir_simulation_setup_add_project_file did not append correctly" << std::endl;
      halir_simulation_setup_free(work);
      return 1;
    }

    halir_simulation_setup_free(work);
    return 0;
  }

  if (test == PROJECT_API_INVALID_ARGS) {
    halir_simulation_setup *work = halir_simulation_setup_create();
    if (work == nullptr) {
      std::cout << "halir_simulation_setup_create returned NULL" << std::endl;
      return 1;
    }

    if (halir_simulation_setup_set_project(nullptr, "a", "b", "c") == 0 ||
        halir_simulation_setup_set_project(work, nullptr, "b", "c") == 0 ||
        halir_simulation_setup_add_project_file(nullptr, "a.txt") == 0 ||
        halir_simulation_setup_add_project_file(work, nullptr) == 0) {
      std::cout << "project API invalid argument checks failed" << std::endl;
      halir_simulation_setup_free(work);
      return 1;
    }

    halir_simulation_setup_free(work);
    return 0;
  }

  if (test == SAMPLE_ENV_SETTER_API) {
    halir_simulation_setup *work = halir_simulation_setup_create();
    if (work == nullptr) {
      std::cout << "halir_simulation_setup_create returned NULL" << std::endl;
      return 1;
    }

    int rc = halir_simulation_setup_set_sample_env(work,
                                            303.15, K,
                                            1.0, ATM,
                                            300.0, M,
                                            2000.0, 2245.0,
                                            0.1, 0.1,
                                            HALIR_BOXCAR,
                                            HALIR_TRANSMISSION,
                                            "");
    if (rc != 0) {
      std::cout << "halir_simulation_setup_set_sample_env failed" << std::endl;
      halir_simulation_setup_free(work);
      return 1;
    }

    if (work->tempU != K || work->pressU != ATM || work->pathLU != CM ||
        work->apod != HALIR_BOXCAR || work->ftype != HALIR_TRANSMISSION) {
      std::cout << "halir_simulation_setup_set_sample_env did not set canonical metadata" << std::endl;
      halir_simulation_setup_free(work);
      return 1;
    }

    halir_simulation_setup_free(work);
    return 0;
  }

  if (test == SAMPLE_ENV_UNIT_NORMALIZATION) {
    halir_simulation_setup *work = halir_simulation_setup_create();
    int read_err = 0;
    double temp = 30.0;
    double press = 1.01325e3;
    double pathL = 300.0;
    halir_Units temp_u = C;
    halir_Units press_u = MBAR;
    halir_Units path_u = M;
    double expected_temp;
    double expected_press;
    double expected_pathL;

    if (work == nullptr) {
      std::cout << "halir_simulation_setup_create returned NULL" << std::endl;
      return 1;
    }

    expected_temp = halir_Units_to_Hitran(&temp_u, &temp, &read_err);
    expected_press = halir_Units_to_Hitran(&press_u, &press, &read_err);
    expected_pathL = halir_Units_to_Hitran(&path_u, &pathL, &read_err);

    int rc = halir_simulation_setup_set_sample_env(work,
                                            temp, C,
                                            press, MBAR,
                                            pathL, M,
                                            2000.0, 2245.0,
                                            0.1, 0.1,
                                            HALIR_TRIANGLE,
                                            HALIR_ABSORBANCE,
                                            "bg.dat");
    if (rc != 0) {
      std::cout << "halir_simulation_setup_set_sample_env conversion case failed" << std::endl;
      halir_simulation_setup_free(work);
      return 1;
    }

    if (!compare_halir_num(&work->temp, &expected_temp) ||
        !compare_halir_num(&work->press, &expected_press) ||
        !compare_halir_num(&work->pathL, &expected_pathL) ||
        work->tempU != K || work->pressU != ATM || work->pathLU != CM ||
        work->ROI[0] != 2000.0 || work->ROI[1] != 2245.0 ||
        work->apod != HALIR_TRIANGLE || work->ftype != HALIR_ABSORBANCE ||
        strncmp(work->bgfile, "bg.dat", PATH_LEN) != 0) {
      std::cout << "halir_simulation_setup_set_sample_env normalization mismatch" << std::endl;
      halir_simulation_setup_free(work);
      return 1;
    }

    halir_simulation_setup_free(work);
    return 0;
  }

  if (test == SAMPLE_ENV_INVALID_ARGS) {
    halir_simulation_setup *work = halir_simulation_setup_create();
    if (work == nullptr) {
      std::cout << "halir_simulation_setup_create returned NULL" << std::endl;
      return 1;
    }

    if (halir_simulation_setup_set_sample_env(nullptr,
                                       300.0, K,
                                       1.0, ATM,
                                       100.0, CM,
                                       2000.0, 2100.0,
                                       0.1, 0.0,
                                       HALIR_BOXCAR,
                                       HALIR_TRANSMISSION,
                                       "") == 0 ||
        halir_simulation_setup_set_sample_env(work,
                                       300.0, K,
                                       1.0, ATM,
                                       100.0, CM,
                                       2100.0, 2000.0,
                                       0.1, 0.0,
                                       HALIR_BOXCAR,
                                       HALIR_TRANSMISSION,
                                       "") == 0 ||
        halir_simulation_setup_set_sample_env(work,
                                       300.0, K,
                                       1.0, ATM,
                                       100.0, CM,
                                       2000.0, 2100.0,
                                       -1.0, 0.0,
                                       HALIR_BOXCAR,
                                       HALIR_TRANSMISSION,
                                       "") == 0 ||
        halir_simulation_setup_set_sample_env(work,
                                       300.0, K,
                                       1.0, ATM,
                                       100.0, CM,
                                       2000.0, 2100.0,
                                       0.1, 0.0,
                                       HALIR_BOXCAR,
                                       HALIR_TRANSMISSION,
                                       nullptr) == 0) {
      std::cout << "sample env invalid argument checks failed" << std::endl;
      halir_simulation_setup_free(work);
      return 1;
    }

    halir_simulation_setup_free(work);
    return 0;
  }

  if (test == COMPOSITION_ADD_API) {
    halir_simulation_setup *work = halir_simulation_setup_create();
    size_t idx1, idx2;
    if (work == nullptr) {
      std::cout << "halir_simulation_setup_create returned NULL" << std::endl;
      return 1;
    }

    if (halir_simulation_setup_add_composition(work, "CO", "Natural", &idx1) != 0) {
      std::cout << "halir_simulation_setup_add_composition failed for first composition" << std::endl;
      halir_simulation_setup_free(work);
      return 1;
    }

    if (idx1 != 0 || work->composition_length != 1) {
      std::cout << "halir_simulation_setup_add_composition index/length mismatch" << std::endl;
      halir_simulation_setup_free(work);
      return 1;
    }

    if (halir_simulation_setup_add_composition(work, "H2O", "Natural", &idx2) != 0) {
      std::cout << "halir_simulation_setup_add_composition failed for second composition" << std::endl;
      halir_simulation_setup_free(work);
      return 1;
    }

    if (idx2 != 1 || work->composition_length != 2) {
      std::cout << "halir_simulation_setup_add_composition second add failed" << std::endl;
      halir_simulation_setup_free(work);
      return 1;
    }

    if (strncmp(work->composition[0]->molec, "CO", 6) != 0 ||
        strncmp(work->composition[0]->isotop, "Natural", 12) != 0 ||
        strncmp(work->composition[1]->molec, "H2O", 6) != 0 ||
        strncmp(work->composition[1]->isotop, "Natural", 12) != 0) {
      std::cout << "halir_simulation_setup_add_composition molec/isotop mismatch" << std::endl;
      halir_simulation_setup_free(work);
      return 1;
    }

    halir_simulation_setup_free(work);
    return 0;
  }

  if (test == COMPOSITION_SET_VMR_API) {
    halir_simulation_setup *work = halir_simulation_setup_create();
    size_t idx;
    double expected_vmr = 5.921539600296e-05;
    if (work == nullptr) {
      std::cout << "halir_simulation_setup_create returned NULL" << std::endl;
      return 1;
    }

    if (halir_simulation_setup_add_composition(work, "CO", "Natural", &idx) != 0) {
      halir_simulation_setup_free(work);
      return 1;
    }

    if (halir_compound_set_vmr(work, idx, expected_vmr) != 0) {
      std::cout << "halir_compound_set_vmr failed" << std::endl;
      halir_simulation_setup_free(work);
      return 1;
    }

    if (!compare_halir_num(&work->composition[idx]->vmr, &expected_vmr)) {
      std::cout << "halir_compound_set_vmr value mismatch" << std::endl;
      halir_simulation_setup_free(work);
      return 1;
    }

    halir_simulation_setup_free(work);
    return 0;
  }

  if (test == COMPOSITION_SET_CONCENTRATION_API) {
    halir_simulation_setup *work = halir_simulation_setup_create();
    size_t idx;
    double expected_conc;
    int read_err = 0;
    double input_conc = 100.0;
    halir_Units conc_units = PPM;

    if (work == nullptr) {
      std::cout << "halir_simulation_setup_create returned NULL" << std::endl;
      return 1;
    }

    if (halir_simulation_setup_add_composition(work, "CO2", "Natural", &idx) != 0) {
      halir_simulation_setup_free(work);
      return 1;
    }

    expected_conc = halir_Units_to_Hitran(&conc_units, &input_conc, &read_err);

    if (halir_compound_set_concentration(work, idx, input_conc, PPM) != 0) {
      std::cout << "halir_compound_set_concentration failed" << std::endl;
      halir_simulation_setup_free(work);
      return 1;
    }

    if (!compare_halir_num(&work->composition[idx]->conc, &expected_conc) ||
        work->composition[idx]->concU != PPM) {
      std::cout << "halir_compound_set_concentration value/unit mismatch" << std::endl;
      halir_simulation_setup_free(work);
      return 1;
    }

    halir_simulation_setup_free(work);
    return 0;
  }

  if (test == COMPOSITION_API_INVALID_ARGS) {
    halir_simulation_setup *work = halir_simulation_setup_create();
    size_t idx;
    if (work == nullptr) {
      std::cout << "halir_simulation_setup_create returned NULL" << std::endl;
      return 1;
    }

    if (halir_simulation_setup_add_composition(nullptr, "CO", "Natural", &idx) == 0 ||
        halir_simulation_setup_add_composition(work, nullptr, "Natural", &idx) == 0 ||
        halir_simulation_setup_add_composition(work, "CO", nullptr, &idx) == 0 ||
        halir_simulation_setup_add_composition(work, "CO", "Natural", nullptr) == 0) {
      std::cout << "add_composition null argument checks failed" << std::endl;
      halir_simulation_setup_free(work);
      return 1;
    }

    if (halir_simulation_setup_add_composition(work, "CO", "Natural", &idx) != 0) {
      halir_simulation_setup_free(work);
      return 1;
    }

    if (halir_compound_set_vmr(nullptr, idx, 0.0) == 0 ||
        halir_compound_set_vmr(work, 999, 0.0) == 0 ||
        halir_compound_set_vmr(work, idx, -0.1) == 0 ||
        halir_compound_set_vmr(work, idx, 1.5) == 0) {
      std::cout << "set_vmr argument bounds checks failed" << std::endl;
      halir_simulation_setup_free(work);
      return 1;
    }

    if (halir_compound_set_concentration(nullptr, idx, 0.0, PPM) == 0 ||
        halir_compound_set_concentration(work, 999, 0.0, PPM) == 0 ||
        halir_compound_set_concentration(work, idx, -1.0, PPM) == 0) {
      std::cout << "set_concentration argument bounds checks failed" << std::endl;
      halir_simulation_setup_free(work);
      return 1;
    }

    halir_simulation_setup_free(work);
    return 0;
  }


  if (test == WORKSPACE_VALIDATION_SUCCESS) {
    halir_simulation_setup *work = halir_simulation_setup_create();
    size_t comp_idx;
    if (work == nullptr) {
      std::cout << "halir_simulation_setup_create returned NULL" << std::endl;
      return 1;
    }

    if (halir_simulation_setup_set_project(work, "Test", "/tmp/", "test") != 0 ||
        halir_simulation_setup_set_sample_env(work, 300.0, K, 1.0, ATM, 100.0, CM, 2000.0, 2100.0, 0.1, 0.0, HALIR_BOXCAR, HALIR_TRANSMISSION, "") != 0 ||
        halir_simulation_setup_add_composition(work, "CO", "Natural", &comp_idx) != 0 ||
        halir_compound_set_vmr(work, comp_idx, 0.0001) != 0) {
      std::cout << "Failed to setup valid workspace" << std::endl;
      halir_simulation_setup_free(work);
      return 1;
    }

    if (halir_simulation_setup_validate(work) != 0) {
      std::cout << "Valid workspace failed validation" << std::endl;
      halir_simulation_setup_free(work);
      return 1;
    }

    halir_simulation_setup_free(work);
    return 0;
  }

  if (test == WORKSPACE_VALIDATION_INVALID_ARGS) {
    halir_simulation_setup *work = halir_simulation_setup_create();
    if (work == nullptr) {
      std::cout << "halir_simulation_setup_create returned NULL" << std::endl;
      return 1;
    }

    if (halir_simulation_setup_validate(nullptr) == 0 ||
        halir_simulation_setup_validate(work) == 0) {
      std::cout << "validation null checks failed" << std::endl;
      halir_simulation_setup_free(work);
      return 1;
    }

    size_t comp_idx;
    if (halir_simulation_setup_add_composition(work, "CO", "Natural", &comp_idx) != 0) {
      halir_simulation_setup_free(work);
      return 1;
    }
    halir_compound_set_vmr(work, comp_idx, 0.0001);

    if (halir_simulation_setup_validate(work) == 0) {
      std::cout << "validation should fail for unset project" << std::endl;
      halir_simulation_setup_free(work);
      return 1;
    }

    if (halir_simulation_setup_set_project(work, "Test", "/tmp/", "test") != 0) {
      halir_simulation_setup_free(work);
      return 1;
    }

    if (halir_simulation_setup_validate(work) == 0) {
      std::cout << "validation should fail for unset sample env" << std::endl;
      halir_simulation_setup_free(work);
      return 1;
    }

    if (halir_simulation_setup_set_sample_env(work, 300.0, K, 1.0, ATM, 100.0, CM, 2100.0, 2000.0, 0.1, 0.0, HALIR_BOXCAR, HALIR_TRANSMISSION, "") == 0) {
      std::cout << "set_sample_env should reject invalid ROI" << std::endl;
      halir_simulation_setup_free(work);
      return 1;
    }

    if (halir_simulation_setup_set_sample_env(work, 300.0, K, 1.0, ATM, 100.0, CM, 2000.0, 2100.0, 0.1, 0.0, HALIR_BOXCAR, HALIR_TRANSMISSION, "") != 0) {
      halir_simulation_setup_free(work);
      return 1;
    }

    if (halir_simulation_setup_validate(work) != 0) {
      std::cout << "validation failed for fully set workspace" << std::endl;
      halir_simulation_setup_free(work);
      return 1;
    }

    halir_simulation_setup_free(work);
    return 0;
  }

  if (test == RESULT_CREATE_FREE_API) {
    halir_simulation_setup *work = halir_simulation_setup_create();
    halir_spectra *spectra = halir_spectra_create(16);
    halir_result *result;

    if (work == nullptr || spectra == nullptr) {
      std::cout << "result lifecycle allocation failed" << std::endl;
      halir_spectra_free(spectra);
      halir_simulation_setup_free(work);
      return 1;
    }

    if (spectra->ndatapnts != 16 || spectra->wavenum == nullptr || spectra->data == nullptr) {
      std::cout << "halir_spectra_create did not initialize buffers" << std::endl;
      halir_spectra_free(spectra);
      halir_simulation_setup_free(work);
      return 1;
    }

    result = halir_result_create(work, 2);
    if (result == nullptr || result->workspace != work || result->nspectra != 2 || result->spectra == nullptr) {
      std::cout << "halir_result_create did not initialize container" << std::endl;
      halir_spectra_free(spectra);
      halir_result_free(result);
      halir_simulation_setup_free(work);
      return 1;
    }

    halir_spectra_free(spectra);
    halir_result_free(result);
    halir_simulation_setup_free(work);
    return 0;
  }

  if (test == RESULT_API_INVALID_ARGS) {
    halir_simulation_setup *work = halir_simulation_setup_create();

    if (work == nullptr) {
      std::cout << "halir_simulation_setup_create returned NULL" << std::endl;
      return 1;
    }

    if (halir_spectra_create(0) != nullptr ||
        halir_result_create(nullptr, 1) != nullptr ||
        halir_result_create(work, 0) != nullptr) {
      std::cout << "result API invalid-argument checks failed" << std::endl;
      halir_simulation_setup_free(work);
      return 1;
    }

    halir_spectra_free(nullptr);
    halir_result_free(nullptr);
    halir_simulation_setup_free(work);
    return 0;
  }

  if (test == SPECTRA_SPC_ROUNDTRIP) {
    const char *path = "/tmp/halir_spectra_roundtrip.spc";
    halir_spectra *spectra = halir_spectra_create(5);
    halir_spectra *loaded = nullptr;

    if (spectra == nullptr) {
      std::cout << "failed to allocate spectra for spc test" << std::endl;
      return 1;
    }

    for (size_t i = 0; i < spectra->ndatapnts; i++) {
      spectra->wavenum[i] = 2000.0 + 0.5 * static_cast<double>(i);
      spectra->data[i] = 0.1 * static_cast<double>(i + 1);
    }

    if (halir_spectra_write_spc(spectra, path) != 0) {
      std::cout << "halir_spectra_write_spc failed" << std::endl;
      halir_spectra_free(spectra);
      return 1;
    }

    loaded = halir_spectra_read_spc(path);
    std::remove(path);
    if (!compare_spc_roundtrip_values(spectra, loaded)) {
      std::cout << "spectra spc roundtrip mismatch" << std::endl;
      halir_spectra_free(spectra);
      halir_spectra_free(loaded);
      return 1;
    }

    halir_spectra_free(spectra);
    halir_spectra_free(loaded);
    return 0;
  }

  if (test == RESULT_SPC_ROUNDTRIP) {
    const char *path = "/tmp/halir_result_roundtrip.spc";
    halir_simulation_setup *work = halir_simulation_setup_create();
    halir_result *result = nullptr;
    halir_result *loaded = nullptr;

    if (work == nullptr) {
      std::cout << "failed to allocate workspace for spc result test" << std::endl;
      return 1;
    }

    work->ftype = HALIR_TRANSMISSION;
    result = halir_result_create(work, 2);
    if (result == nullptr) {
      std::cout << "failed to allocate result for spc test" << std::endl;
      halir_simulation_setup_free(work);
      return 1;
    }

    for (size_t s = 0; s < result->nspectra; s++) {
      const size_t npts = 4 + s;
      result->spectra[s].ndatapnts = npts;
      result->spectra[s].wavenum = (halir_num*)calloc(npts, sizeof(halir_num));
      result->spectra[s].data = (halir_num*)calloc(npts, sizeof(halir_num));
      if ((result->spectra[s].wavenum == nullptr) || (result->spectra[s].data == nullptr)) {
        std::cout << "failed to allocate result spectra buffers" << std::endl;
        halir_result_free(result);
        halir_simulation_setup_free(work);
        return 1;
      }

      for (size_t i = 0; i < npts; i++) {
        result->spectra[s].wavenum[i] = 2100.0 + 0.25 * static_cast<double>(i) + static_cast<double>(s);
        result->spectra[s].data[i] = 1.0 + static_cast<double>(s) + 0.05 * static_cast<double>(i);
      }
    }

    if (halir_result_write_spc(result, path) != 0) {
      std::cout << "halir_result_write_spc failed" << std::endl;
      halir_result_free(result);
      halir_simulation_setup_free(work);
      return 1;
    }

    loaded = halir_result_read_spc(path);
    std::remove(path);
    if (!compare_spc_roundtrip_result(result, loaded)) {
      std::cout << "result spc roundtrip mismatch" << std::endl;
      halir_result_free(result);
      halir_result_free(loaded);
      halir_simulation_setup_free(work);
      return 1;
    }

    halir_result_free(result);
    halir_result_free(loaded);
    halir_simulation_setup_free(work);
    return 0;
  }

  if (test == RESULT_SPC_DIRECTORY_READ) {
    const char *path = "/tmp/halir_result_directory_read.spc";
    halir_result *loaded = nullptr;
    FILE *fp = fopen(path, "wb");
    SPCHDR header;
    SUBHDR subheader;
    SSFSTC dir[2];
    float x0[3] = {2100.0f, 2100.5f, 2101.0f};
    float y0[3] = {1.0f, 1.5f, 2.0f};
    float x1[2] = {2200.0f, 2200.25f};
    float y1[2] = {3.0f, 4.0f};
    long first_offset;
    long second_offset;
    long dir_offset;
    halir_num expected_x0_0;
    halir_num expected_x1_0;
    halir_num expected_y0_2;
    halir_num expected_y1_1;

    if (fp == nullptr) {
      std::cout << "failed to create spc directory test file" << std::endl;
      return 1;
    }

    std::memset(&header, 0, sizeof(header));
    header.ftflgs = static_cast<BYTE>(TMULTI | TXVALS | TXYXYS);
    header.fversn = 0x4B;
    header.fexper = SPCFTIR;
    header.fexp = static_cast<char>(0x80);
    header.fnsub = 2;
    header.fxtype = XWAVEN;
    header.fytype = YARB;
    header.ffirst = 2100.0;
    header.flast = 2200.25;
    if (fwrite(&header, sizeof(header), 1, fp) != 1) {
      fclose(fp);
      std::remove(path);
      return 1;
    }

    second_offset = std::ftell(fp);
    std::memset(&subheader, 0, sizeof(subheader));
    subheader.subexp = static_cast<char>(0x80);
    subheader.subindx = 1;
    subheader.subtime = 1.0f;
    subheader.subnext = 2.0f;
    subheader.subnpts = 2;
    fwrite(&subheader, sizeof(subheader), 1, fp);
    fwrite(x1, sizeof(float), 2, fp);
    fwrite(y1, sizeof(float), 2, fp);

    first_offset = std::ftell(fp);
    std::memset(&subheader, 0, sizeof(subheader));
    subheader.subexp = static_cast<char>(0x80);
    subheader.subindx = 0;
    subheader.subtime = 0.0f;
    subheader.subnext = 1.0f;
    subheader.subnpts = 3;
    fwrite(&subheader, sizeof(subheader), 1, fp);
    fwrite(x0, sizeof(float), 3, fp);
    fwrite(y0, sizeof(float), 3, fp);

    dir_offset = std::ftell(fp);
    std::memset(dir, 0, sizeof(dir));
    dir[0].ssfposn = static_cast<DWORD>(first_offset);
    dir[0].ssfsize = static_cast<DWORD>(sizeof(SUBHDR) + sizeof(x0) + sizeof(y0));
    dir[0].ssftime = 0.0f;
    dir[1].ssfposn = static_cast<DWORD>(second_offset);
    dir[1].ssfsize = static_cast<DWORD>(sizeof(SUBHDR) + sizeof(x1) + sizeof(y1));
    dir[1].ssftime = 1.0f;
    fwrite(dir, sizeof(SSFSTC), 2, fp);

    header.fnpts = static_cast<DWORD>(dir_offset);
    std::fseek(fp, 0, SEEK_SET);
    fwrite(&header, sizeof(header), 1, fp);
    fclose(fp);

    loaded = halir_result_read_spc(path);
    std::remove(path);
    if (loaded == nullptr || loaded->nspectra != 2) {
      std::cout << "failed to read directory-based spc result" << std::endl;
      halir_result_free(loaded);
      return 1;
    }

    if (loaded->spectra[0].ndatapnts != 3 || loaded->spectra[1].ndatapnts != 2) {
      std::cout << "directory-based spc result has wrong subfile sizes" << std::endl;
      halir_result_free(loaded);
      return 1;
    }

    expected_x0_0 = static_cast<halir_num>(x0[0]);
    expected_x1_0 = static_cast<halir_num>(x1[0]);
    expected_y0_2 = static_cast<halir_num>(y0[2]);
    expected_y1_1 = static_cast<halir_num>(y1[1]);

    if (!compare_halir_num(&loaded->spectra[0].wavenum[0], &expected_x0_0) ||
        !compare_halir_num(&loaded->spectra[1].wavenum[0], &expected_x1_0) ||
        !compare_halir_num(&loaded->spectra[0].data[2], &expected_y0_2) ||
        !compare_halir_num(&loaded->spectra[1].data[1], &expected_y1_1)) {
      std::cout << "directory-based spc result values mismatch" << std::endl;
      halir_result_free(loaded);
      return 1;
    }

    halir_result_free(loaded);
    return 0;
  }
  std::string pname = "CO_Test2";
  std::string temp = "303.15";
  std::string temp_unit = "K";
  std::string press = "1.0";
  std::string press_unit = "atm";
  std::string roi = "[2000.0, 2245.0]";

  const std::string ref1_str = build_input_json(pname, temp, temp_unit, press, press_unit, roi);
  halir_simulation_setup *work_p1;
  halir_simulation_setup *work_ref = halir_parseJSONinput(ref1_str.c_str());
  string err_msg;

  switch (test) {
    case TEMP_F_TO_K:
      temp = "86";
      temp_unit = "F";
      err_msg = "Error in temperature conversion F->K\n";
      break;
    case TEMP_C_TO_K:
      temp = "30";
      temp_unit = "C";
      err_msg =  "Error in temperature conversion C->K\n";
      break;
    case PRESS_BAR_TO_ATM:
      press = "1.0132500";
      press_unit = "bar";
      err_msg = "Error in pressure conversion bar->atm\n";
      break;
    case PRESS_MBAR_TO_ATM:
      press = "1.0132500e3";
      press_unit = "mbar";
      err_msg = "Error in pressure conversion mbar->atm\n";
      break;
    case PRESS_PA_TO_ATM:
      press = "1.0132500e5";
      press_unit = "pa";
      err_msg = "Error in pressure conversion pa->atm\n";
      break;
    case PRESS_HPA_TO_ATM:
      press = "1.0132500e3";
      press_unit = "hpa";
      err_msg = "Error in pressure conversion hpa->atm\n";
      break;
    case PRESS_MMHG_TO_ATM:
      press = "760.00000";
      press_unit = "mmhg";
      err_msg = "Error in pressure conversion mmHg->atm\n";
      //if (!compare_halir_workspace(work_ref, work_p1)) {
        //std::cout.precision(13);
        //std::cout << std::scientific;
        //std::cout << "ref: " << work_ref->press << " inp: " << work_p1->press << std::endl;
        //return 1;
      //}
      break;
    case ROI_TOO_MANY_VALUES:
      roi = "[2000.0, 2245.0, 2300.0]";
      break;
    case ROI_NON_NUMERIC:
      roi = "[2000.0, \"bad\"]";
      break;
    case PNAME_TOO_LONG:
      pname = std::string(PATH_LEN + 32, 'a');
      break;
  }
  if (work_ref == NULL) {
    std::cout << "Reference parser input should be valid" << std::endl;
    return 1;
  }

  std::string inp2_str = build_input_json(pname, temp, temp_unit, press, press_unit, roi);
  if (inp2_str.empty()) {
    std::cout << "Failed to prepare parser input" << std::endl;
    halir_simulation_setup_free(work_ref);
    return 1;
  }

  work_p1 = halir_parseJSONinput(inp2_str.c_str());

  if (test == ROI_TOO_MANY_VALUES || test == ROI_NON_NUMERIC || test == PNAME_TOO_LONG) {
    if (work_p1 != NULL) {
      std::cout << "Expected parser to reject invalid input" << std::endl;
      halir_simulation_setup_free(work_p1);
      halir_simulation_setup_free(work_ref);
      return 1;
    }
    halir_simulation_setup_free(work_ref);
    return 0;
  }

  if (!compare_halir_workspace(work_ref, work_p1)) {
    std::cout << err_msg << std::endl;
    halir_simulation_setup_free(work_p1);
    halir_simulation_setup_free(work_ref);
    return 1;
  }

  halir_simulation_setup_free(work_p1);
  halir_simulation_setup_free(work_ref);

  //halir_print_simulation_setup(work_ref);
  return 0;
}
