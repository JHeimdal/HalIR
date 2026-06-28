#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

extern "C" {
#include "HalIR/halir.h"
#include "HalIR/halir_cls.h"
}

static bool approx(halir_num a, halir_num b, halir_num tol = 1e-9) {
  return std::abs(a - b) <= tol * std::max((halir_num)1.0, std::abs(b));
}

int main(int argc, char **argv)
{
  enum TEST {
    MATRIX_LIFECYCLE,
    GRID_BUILD,
    RESAMPLE_LINEAR,
    SYNTHESIZE_A,
    CALIBRATE_RECOVER,
    PREDICT_RECOVER,
    INVALID_ARGS,
  };

  using namespace std;

  vector<string> args(argv + 1, argv + argc);
  if (args.empty()) {
    cout << "missing test selector" << endl;
    return 1;
  }
  TEST test = (TEST)stoi(args[0]);

  if (test == MATRIX_LIFECYCLE) {
    halir_matrix *m = halir_matrix_create(2, 3);
    if (m == nullptr || m->rows != 2 || m->cols != 3) {
      cout << "matrix create failed" << endl;
      return 1;
    }
    halir_matrix_set(m, 1, 2, 7.5);
    if (!approx(halir_matrix_get(m, 1, 2), 7.5) || !approx(halir_matrix_get(m, 0, 0), 0.0)) {
      cout << "matrix get/set failed" << endl;
      halir_matrix_free(m);
      return 1;
    }
    halir_matrix_free(m);
    halir_matrix_free(nullptr);
    if (halir_matrix_create(0, 3) != nullptr) {
      cout << "matrix create should reject zero dim" << endl;
      return 1;
    }
    return 0;
  }

  if (test == GRID_BUILD) {
    halir_cls_window win[2] = {{10.0, 12.0}, {20.0, 21.0}};
    halir_cls_grid *grid = halir_cls_build_grid(win, 2, 1.0);
    if (grid == nullptr) {
      cout << "build_grid returned NULL" << endl;
      return 1;
    }
    // window 0: 10,11,12 (3 pts); window 1: 20,21 (2 pts) => 5
    if (grid->n != 5) {
      cout << "build_grid wrong size: " << grid->n << endl;
      halir_cls_grid_free(grid);
      return 1;
    }
    if (!approx(grid->wavenum[0], 10.0) || !approx(grid->wavenum[2], 12.0) ||
        !approx(grid->wavenum[3], 20.0) || !approx(grid->wavenum[4], 21.0)) {
      cout << "build_grid wrong values" << endl;
      halir_cls_grid_free(grid);
      return 1;
    }
    halir_cls_grid_free(grid);
    return 0;
  }

  if (test == RESAMPLE_LINEAR) {
    // Source: y = 2x on [0,10]
    vector<halir_num> sw, sd;
    for (int i = 0; i <= 10; i++) { sw.push_back(i); sd.push_back(2.0 * i); }
    halir_cls_window win = {2.0, 8.0};
    halir_cls_grid *grid = halir_cls_build_grid(&win, 1, 0.5);
    if (grid == nullptr) { cout << "grid build failed" << endl; return 1; }
    vector<halir_num> out(grid->n, -1.0);
    if (halir_cls_resample(sw.data(), sd.data(), sw.size(), grid, out.data()) != 0) {
      cout << "resample failed" << endl;
      halir_cls_grid_free(grid);
      return 1;
    }
    for (size_t i = 0; i < grid->n; i++) {
      if (!approx(out[i], 2.0 * grid->wavenum[i], 1e-9)) {
        cout << "resample linear mismatch at " << i << endl;
        halir_cls_grid_free(grid);
        return 1;
      }
    }
    halir_cls_grid_free(grid);
    // Out-of-range -> zero
    halir_cls_grid *g2 = halir_cls_grid_create(2);
    g2->wavenum[0] = -5.0; g2->wavenum[1] = 100.0;
    vector<halir_num> o2(2, 9.0);
    halir_cls_resample(sw.data(), sd.data(), sw.size(), g2, o2.data());
    if (!approx(o2[0], 0.0) || !approx(o2[1], 0.0)) {
      cout << "resample out-of-range not zeroed" << endl;
      halir_cls_grid_free(g2);
      return 1;
    }
    halir_cls_grid_free(g2);
    return 0;
  }

  if (test == SYNTHESIZE_A) {
    // C (2x2), R (2x3); A = C*R
    halir_matrix *C = halir_matrix_create(2, 2);
    halir_matrix *R = halir_matrix_create(2, 3);
    halir_num cvals[4] = {1, 2, 3, 4};
    halir_num rvals[6] = {1, 0, 1, 0, 1, 1};
    for (int i = 0; i < 4; i++) C->data[i] = cvals[i];
    for (int i = 0; i < 6; i++) R->data[i] = rvals[i];
    halir_matrix *A = halir_cls_synthesize_A(C, R);
    if (A == nullptr || A->rows != 2 || A->cols != 3) {
      cout << "synthesize_A failed" << endl;
      return 1;
    }
    // Row0 = 1*[1,0,1] + 2*[0,1,1] = [1,2,3]
    // Row1 = 3*[1,0,1] + 4*[0,1,1] = [3,4,7]
    halir_num expect[6] = {1, 2, 3, 3, 4, 7};
    for (int i = 0; i < 6; i++) {
      if (!approx(A->data[i], expect[i])) {
        cout << "synthesize_A wrong value at " << i << endl;
        halir_matrix_free(A); halir_matrix_free(C); halir_matrix_free(R);
        return 1;
      }
    }
    halir_matrix_free(A); halir_matrix_free(C); halir_matrix_free(R);
    return 0;
  }

  if (test == CALIBRATE_RECOVER) {
    const size_t nc = 2, nw = 4, ns = 3;
    halir_matrix *R = halir_matrix_create(nc, nw);
    halir_num rvals[8] = {1, 2, 3, 4, 0.5, 0.0, 1.0, 2.0};
    for (size_t i = 0; i < nc * nw; i++) R->data[i] = rvals[i];
    halir_matrix *C = halir_matrix_create(ns, nc);
    halir_num cvals[6] = {1, 0, 0, 1, 1, 1};
    for (size_t i = 0; i < ns * nc; i++) C->data[i] = cvals[i];
    halir_matrix *A = halir_cls_synthesize_A(C, R);
    halir_cls_grid *grid = halir_cls_grid_create(nw);
    for (size_t i = 0; i < nw; i++) grid->wavenum[i] = 10.0 + i;

    halir_cls_calibration *cal = halir_cls_calibrate(A, C, grid);
    if (cal == nullptr) {
      cout << "calibrate returned NULL" << endl;
      return 1;
    }
    int rc = 0;
    for (size_t i = 0; i < nc * nw; i++) {
      if (!approx(cal->K->data[i], R->data[i], 1e-7)) {
        cout << "calibrate did not recover R at " << i
             << " got " << cal->K->data[i] << " want " << R->data[i] << endl;
        rc = 1; break;
      }
    }
    halir_cls_calibration_free(cal);
    halir_matrix_free(A); halir_matrix_free(C); halir_matrix_free(R);
    halir_cls_grid_free(grid);
    return rc;
  }

  if (test == PREDICT_RECOVER) {
    const size_t nc = 2, nw = 4, ns = 3;
    halir_matrix *R = halir_matrix_create(nc, nw);
    halir_num rvals[8] = {1, 2, 3, 4, 0.5, 0.0, 1.0, 2.0};
    for (size_t i = 0; i < nc * nw; i++) R->data[i] = rvals[i];
    halir_matrix *C = halir_matrix_create(ns, nc);
    halir_num cvals[6] = {1, 0, 0, 1, 1, 1};
    for (size_t i = 0; i < ns * nc; i++) C->data[i] = cvals[i];
    halir_matrix *A = halir_cls_synthesize_A(C, R);
    halir_cls_grid *grid = halir_cls_grid_create(nw);
    for (size_t i = 0; i < nw; i++) grid->wavenum[i] = 10.0 + i;
    halir_cls_calibration *cal = halir_cls_calibrate(A, C, grid);
    if (cal == nullptr) { cout << "calibrate failed" << endl; return 1; }

    // Known concentrations; build sample = conc * K on the grid.
    halir_num conc_true[2] = {2.0, 3.0};
    vector<halir_num> sdata(nw, 0.0), swave(nw, 0.0);
    for (size_t k = 0; k < nw; k++) {
      swave[k] = grid->wavenum[k];
      for (size_t r = 0; r < nc; r++) sdata[k] += conc_true[r] * R->data[r * nw + k];
    }
    halir_spectra sample;
    std::memset(&sample, 0, sizeof(sample));
    sample.ndatapnts = nw;
    sample.wavenum = swave.data();
    sample.data = sdata.data();

    halir_cls_prediction *pred = halir_cls_predict(cal, &sample);
    int rc = 0;
    if (pred == nullptr) {
      cout << "predict returned NULL" << endl;
      rc = 1;
    } else {
      if (!approx(pred->conc[0], 2.0, 1e-7) || !approx(pred->conc[1], 3.0, 1e-7)) {
        cout << "predict conc wrong: " << pred->conc[0] << ", " << pred->conc[1] << endl;
        rc = 1;
      }
      for (size_t k = 0; k < nw && rc == 0; k++) {
        if (std::abs(pred->residual.data[k]) > 1e-7) {
          cout << "predict residual nonzero at " << k << ": " << pred->residual.data[k] << endl;
          rc = 1;
        }
        if (!approx(pred->fitted.data[k], sdata[k], 1e-7)) {
          cout << "predict fitted mismatch at " << k << endl;
          rc = 1;
        }
      }
      halir_cls_prediction_free(pred);
    }
    halir_cls_calibration_free(cal);
    halir_matrix_free(A); halir_matrix_free(C); halir_matrix_free(R);
    halir_cls_grid_free(grid);
    return rc;
  }

  if (test == INVALID_ARGS) {
    if (halir_cls_build_grid(nullptr, 1, 1.0) != nullptr) { cout << "build_grid null" << endl; return 1; }
    halir_cls_window w = {1.0, 2.0};
    if (halir_cls_build_grid(&w, 1, 0.0) != nullptr) { cout << "build_grid step0" << endl; return 1; }
    if (halir_cls_synthesize_A(nullptr, nullptr) != nullptr) { cout << "synth null" << endl; return 1; }
    if (halir_cls_calibrate(nullptr, nullptr, nullptr) != nullptr) { cout << "calibrate null" << endl; return 1; }
    if (halir_cls_predict(nullptr, nullptr) != nullptr) { cout << "predict null" << endl; return 1; }
    // dimension mismatch in synthesize: C cols != R rows
    halir_matrix *C = halir_matrix_create(2, 3);
    halir_matrix *R = halir_matrix_create(2, 2);
    if (halir_cls_synthesize_A(C, R) != nullptr) {
      cout << "synth dim mismatch not caught" << endl;
      halir_matrix_free(C); halir_matrix_free(R);
      return 1;
    }
    halir_matrix_free(C); halir_matrix_free(R);
    // NULL-safe frees
    halir_cls_grid_free(nullptr);
    halir_cls_calibration_free(nullptr);
    halir_cls_prediction_free(nullptr);
    return 0;
  }

  cout << "unknown test selector" << endl;
  return 1;
}
