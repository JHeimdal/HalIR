#include <cstdio>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <vector>
#include <string>
#include <boost/json/src.hpp>

extern "C" {
#include "HalIR/halir.h"
}

namespace json = boost::json;

namespace {

// Helper to build minimal workspace with valid sample environment for prmfile loading tests
halir_workspace* build_test_workspace()
{
  halir_workspace *work = halir_workspace_create();
  if (work == nullptr) {
    return nullptr;
  }

  if (halir_workspace_set_project(work, "TestProject", "/tmp/", "") != 0 ||
      halir_workspace_set_sample_env(work, 303.15, K, 1.0, ATM, 300.0, M, 2000.0, 2245.0, 0.1, 0.1, HALIR_BOXCAR, HALIR_TRANSMISSION, "") != 0) {
    halir_workspace_free(work);
    return nullptr;
  }

  return work;
}

bool write_truncated_prmfile(const std::filesystem::path &file_path)
{
  std::ofstream out(file_path, std::ios::binary | std::ios::trunc);
  if (!out.is_open()) {
    return false;
  }

  const char truncated_bytes[] = {0x01, 0x02, 0x03, 0x04};
  out.write(truncated_bytes, sizeof(truncated_bytes));
  return out.good();
}

bool write_header_only_prmfile(const std::filesystem::path &file_path)
{
  halir_HitranHead head = {};
  std::ofstream out(file_path, std::ios::binary | std::ios::trunc);
  if (!out.is_open()) {
    return false;
  }

  head.nisotp = 1;
  head.ndatapnts = 1;
  out.write(reinterpret_cast<const char*>(&head), sizeof(head));
  return out.good();
}

bool write_missing_line_data_prmfile(const std::filesystem::path &file_path)
{
  halir_HitranHead head = {};
  const int molecule_id = 1;
  std::ofstream out(file_path, std::ios::binary | std::ios::trunc);
  if (!out.is_open()) {
    return false;
  }

  head.nisotp = 1;
  head.ndatapnts = 1;
  out.write(reinterpret_cast<const char*>(&head), sizeof(head));
  out.write(reinterpret_cast<const char*>(&molecule_id), sizeof(molecule_id));
  return out.good();
}

} // namespace

int main(int argc, char **argv)
{
  enum TEST {
    VALID_HPAR,
    MISSING_HPAR,
    TRUNCATED_HPAR,
    HEADER_ONLY_HPAR,
    MISSING_LINE_DATA_HPAR,
  };

  using namespace std;

  vector<string> args(argv + 1, argv + argc);
  TEST test = (TEST)stoi(args[0]);

  const std::string valid_prmfile = "/home/jimmy/Programs/HalIR/tests/CO_test/CO.hpar";
  halir_workspace *workspace = build_test_workspace();
  std::filesystem::path truncated_path;
  bool temporary_fixture_created = false;
  std::string test_prmfile = valid_prmfile;

  if (workspace == nullptr) {
    std::cout << "Failed to build test workspace" << std::endl;
    return 1;
  }

  size_t comp_index;
  if (halir_workspace_add_composition(workspace, "CO", "Natural", &comp_index) != 0) {
    std::cout << "Failed to add composition" << std::endl;
    halir_workspace_free(workspace);
    return 1;
  }

  switch (test) {
    case VALID_HPAR:
      test_prmfile = valid_prmfile;
      break;
    case MISSING_HPAR:
      test_prmfile = "/home/jimmy/Programs/HalIR/tests/CO_test/does_not_exist.hpar";
      break;
    case TRUNCATED_HPAR:
      truncated_path = std::filesystem::temp_directory_path() / "halir_truncated_test.hpar";
      if (!write_truncated_prmfile(truncated_path)) {
        std::cout << "Could not create truncated prmfile fixture" << std::endl;
        halir_workspace_free(workspace);
        return 1;
      }
      temporary_fixture_created = true;
      test_prmfile = truncated_path.string();
      break;
    case HEADER_ONLY_HPAR:
      truncated_path = std::filesystem::temp_directory_path() / "halir_header_only_test.hpar";
      if (!write_header_only_prmfile(truncated_path)) {
        std::cout << "Could not create header-only prmfile fixture" << std::endl;
        halir_workspace_free(workspace);
        return 1;
      }
      temporary_fixture_created = true;
      test_prmfile = truncated_path.string();
      break;
    case MISSING_LINE_DATA_HPAR:
      truncated_path = std::filesystem::temp_directory_path() / "halir_missing_line_data_test.hpar";
      if (!write_missing_line_data_prmfile(truncated_path)) {
        std::cout << "Could not create missing-line-data prmfile fixture" << std::endl;
        halir_workspace_free(workspace);
        return 1;
      }
      temporary_fixture_created = true;
      test_prmfile = truncated_path.string();
      break;
  }

  int load_result = halir_compound_load_prmfile(workspace, comp_index, test_prmfile.c_str());

  if (test == VALID_HPAR) {
    if (load_result != 0) {
      std::cout << "Expected valid prmfile to load" << std::endl;
      halir_workspace_free(workspace);
      return 1;
    }
    if (workspace->composition[comp_index]->hitran_prms == nullptr) {
      std::cout << "Expected prmfile spectral lines to load" << std::endl;
      halir_workspace_free(workspace);
      return 1;
    }
    if (workspace->composition[comp_index]->hitran_head.ndatapnts <= 0) {
      std::cout << "Expected prmfile to expose positive line count" << std::endl;
      halir_workspace_free(workspace);
      return 1;
    }
    halir_workspace_free(workspace);
    return 0;
  }

  if (temporary_fixture_created) {
    std::filesystem::remove(truncated_path);
  }

  if (load_result == 0) {
    std::cout << "Expected prmfile failure to return non-zero" << std::endl;
    halir_workspace_free(workspace);
    return 1;
  }

  halir_workspace_free(workspace);
  return 0;
}
