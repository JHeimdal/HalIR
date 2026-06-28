#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

extern "C" {
#include "HalIR/halir.h"
}

#ifndef HALIR_TESTS_SOURCE_DIR
#define HALIR_TESTS_SOURCE_DIR "."
#endif

namespace {

struct Golden {
  size_t n;
  double w0;
  double wN;
  double sum;
  double max;
  double mid;
};

bool write_single_line_prmfile(const std::filesystem::path &file_path,
                               int molec_num,
                               int isotp_num,
                               double trans_mu,
                               double line_I,
                               double air_B,
                               double self_B,
                               double low_state_en,
                               double temp_air_B,
                               double pressure_S,
                               double molec_mass)
{
  halir_HitranHead head = {};
  halir_HitranLine line = {};
  const int molecule_id = molec_num;
  std::ofstream out(file_path, std::ios::binary | std::ios::trunc);

  if (!out.is_open()) {
    return false;
  }

  head.nisotp = 1;
  head.ndatapnts = 1;
  head.roi_low = 2098.0f;
  head.roi_high = 2102.0f;

  line.molec_num = molec_num;
  line.isotp_num = isotp_num;
  line.trans_mu = trans_mu;
  line.line_I = line_I;
  line.air_B = air_B;
  line.self_B = self_B;
  line.low_state_en = low_state_en;
  line.temp_air_B = temp_air_B;
  line.pressure_S = pressure_S;
  line.molecMass = molec_mass;

  out.write(reinterpret_cast<const char*>(&head), sizeof(head));
  out.write(reinterpret_cast<const char*>(&molecule_id), sizeof(molecule_id));
  out.write(reinterpret_cast<const char*>(&line), sizeof(line));
  out.close();
  return !out.fail();
}

bool is_monotonic_increasing(const halir_num *arr, size_t n)
{
  if (arr == nullptr || n < 2) {
    return false;
  }
  for (size_t i = 1; i < n; i++) {
    if (!(arr[i] > arr[i - 1])) {
      return false;
    }
  }
  return true;
}

bool all_finite(const halir_num *arr, size_t n)
{
  if (arr == nullptr) {
    return false;
  }
  for (size_t i = 0; i < n; i++) {
    if (!std::isfinite(arr[i])) {
      return false;
    }
  }
  return true;
}

double relative_or_absolute_error(double got, double expected)
{
  const double abs_err = std::abs(got - expected);
  const double denom = std::max(1.0, std::abs(expected));
  return abs_err / denom;
}

double spectrum_sum(const halir_spectra &s)
{
  double acc = 0.0;
  for (size_t i = 0; i < s.ndatapnts; i++) {
    acc += s.data[i];
  }
  return acc;
}

double spectrum_max(const halir_spectra &s)
{
  double m = -std::numeric_limits<double>::infinity();
  for (size_t i = 0; i < s.ndatapnts; i++) {
    if (s.data[i] > m) {
      m = s.data[i];
    }
  }
  return m;
}

void dump_metrics(const halir_result *result)
{
  for (size_t i = 0; i < result->nspectra; i++) {
    const halir_spectra &s = result->spectra[i];
    std::cout << "spec=" << i
              << " n=" << s.ndatapnts
              << " w0=" << s.wavenum[0]
              << " wN=" << s.wavenum[s.ndatapnts - 1]
              << " sum=" << spectrum_sum(s)
              << " max=" << spectrum_max(s)
              << " mid=" << s.data[s.ndatapnts / 2]
              << std::endl;
  }
}

bool load_golden_metrics(const std::filesystem::path &fixture_path,
                         std::vector<Golden> *metrics_out)
{
  std::ifstream in(fixture_path);
  Golden m = {};

  if ((metrics_out == nullptr) || !in.is_open()) {
    return false;
  }

  metrics_out->clear();
  while (in >> m.n >> m.w0 >> m.wN >> m.sum >> m.max >> m.mid) {
    metrics_out->push_back(m);
  }

  if (!in.eof()) {
    return false;
  }

  return !metrics_out->empty();
}

} // namespace

int main(int argc, char **argv)
{
  const bool dump_mode = (argc > 1) && (std::string(argv[1]) == "--dump");
  const std::filesystem::path prm1 = std::filesystem::temp_directory_path() / "halir_calc_reg_1.hpar";
  const std::filesystem::path prm2 = std::filesystem::temp_directory_path() / "halir_calc_reg_2.hpar";
  const std::filesystem::path fixture_path =
      std::filesystem::path(HALIR_TESTS_SOURCE_DIR) / "fixtures" / "calc_result_metrics.txt";
  const double rel_tol = 5e-3;
  const double abs_tol = 1e-14;
  std::vector<Golden> golden;

  halir_simulation_setup *work = nullptr;
  halir_result *result = nullptr;
  const halir_spectra *s0 = nullptr;
  const halir_spectra *s1 = nullptr;
  const Golden *g0 = nullptr;
  const Golden *g1 = nullptr;
  size_t comp_a = 0;
  size_t comp_b = 0;
  int exit_code = 1;

  if (!write_single_line_prmfile(prm1, 1, 1, 2100.0, 1.0e-20, 0.07, 0.10, 50.0, 0.70, 0.0, 18.0) ||
      !write_single_line_prmfile(prm2, 1, 1, 2100.5, 1.2e-20, 0.08, 0.11, 55.0, 0.70, 0.0, 18.0)) {
    std::cout << "Failed to create prmfile fixtures" << std::endl;
    goto cleanup;
  }

  work = halir_simulation_setup_create();
  if (work == nullptr) {
    std::cout << "halir_simulation_setup_create failed" << std::endl;
    goto cleanup;
  }

  if (halir_simulation_setup_set_project(work, "CalcRegression", "/tmp", "") != 0 ||
      halir_simulation_setup_set_sample_env(work,
                                     296.0, K,
                                     1.0, ATM,
                                     100.0, CM,
                                     2098.0, 2102.0,
                                     0.1, 0.1,
                                     HALIR_BOXCAR,
                                     HALIR_TRANSMISSION,
                                     "") != 0) {
    std::cout << "Failed to initialize workspace metadata" << std::endl;
    goto cleanup;
  }

  if (halir_simulation_setup_add_composition(work, "H2O", "Natural", &comp_a) != 0 ||
      halir_simulation_setup_add_composition(work, "H2O", "Natural", &comp_b) != 0 ||
      halir_compound_set_vmr(work, comp_a, 8.0e-4) != 0 ||
      halir_compound_set_vmr(work, comp_b, 6.0e-4) != 0 ||
      halir_compound_load_prmfile(work, comp_a, prm1.string().c_str()) != 0 ||
      halir_compound_load_prmfile(work, comp_b, prm2.string().c_str()) != 0) {
    std::cout << "Failed to prepare composition inputs" << std::endl;
    goto cleanup;
  }

  result = halir_calculate_result(work);
  if (result == nullptr) {
    std::cout << "halir_calculate_result returned NULL" << std::endl;
    goto cleanup;
  }

  if (result->nspectra != 2 || result->spectra == nullptr) {
    std::cout << "Unexpected result container shape" << std::endl;
    goto cleanup;
  }

  for (size_t i = 0; i < result->nspectra; i++) {
    const halir_spectra &s = result->spectra[i];
    if (s.ndatapnts < 8 || s.wavenum == nullptr || s.data == nullptr) {
      std::cout << "Spectrum buffer not initialized" << std::endl;
      goto cleanup;
    }
    if (!is_monotonic_increasing(s.wavenum, s.ndatapnts)) {
      std::cout << "Wavenumber grid is not strictly increasing" << std::endl;
      goto cleanup;
    }
    if (!all_finite(s.wavenum, s.ndatapnts) || !all_finite(s.data, s.ndatapnts)) {
      std::cout << "Non-finite values in spectrum outputs" << std::endl;
      goto cleanup;
    }
  }

  if (dump_mode) {
    dump_metrics(result);
    exit_code = 0;
    goto cleanup;
  }

  if (!load_golden_metrics(fixture_path, &golden)) {
    std::cout << "Failed to load golden fixture: " << fixture_path << std::endl;
    goto cleanup;
  }
  if (golden.size() != result->nspectra) {
    std::cout << "Golden fixture spectrum count mismatch" << std::endl;
    goto cleanup;
  }

  s0 = &result->spectra[0];
  s1 = &result->spectra[1];

  g0 = &golden[0];
  g1 = &golden[1];

  if (s0->ndatapnts != g0->n || s1->ndatapnts != g1->n) {
    std::cout << "Golden datapoint counts do not match" << std::endl;
    goto cleanup;
  }

  if (std::abs(s0->wavenum[0] - g0->w0) > abs_tol ||
      std::abs(s0->wavenum[s0->ndatapnts - 1] - g0->wN) > abs_tol ||
      relative_or_absolute_error(spectrum_sum(*s0), g0->sum) > rel_tol ||
      relative_or_absolute_error(spectrum_max(*s0), g0->max) > rel_tol ||
      relative_or_absolute_error(s0->data[s0->ndatapnts / 2], g0->mid) > rel_tol) {
    std::cout << "Golden check failed for spectrum 0" << std::endl;
    goto cleanup;
  }

  if (std::abs(s1->wavenum[0] - g1->w0) > abs_tol ||
      std::abs(s1->wavenum[s1->ndatapnts - 1] - g1->wN) > abs_tol ||
      relative_or_absolute_error(spectrum_sum(*s1), g1->sum) > rel_tol ||
      relative_or_absolute_error(spectrum_max(*s1), g1->max) > rel_tol ||
      relative_or_absolute_error(s1->data[s1->ndatapnts / 2], g1->mid) > rel_tol) {
    std::cout << "Golden check failed for spectrum 1" << std::endl;
    goto cleanup;
  }

  exit_code = 0;

cleanup:
  halir_result_free(result);
  halir_simulation_setup_free(work);
  std::filesystem::remove(prm1);
  std::filesystem::remove(prm2);
  return exit_code;
}
