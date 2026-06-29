#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>

#include <HalIR/halir.h>
#include <HalIR/halir_cls.h>

/* ----------------------------------------------------------------------- */
/* Matrix lifecycle                                                        */
/* ----------------------------------------------------------------------- */

halir_matrix *halir_matrix_create(size_t rows, size_t cols)
{
  halir_matrix *m;

  if ((rows == 0) || (cols == 0)) {
    return NULL;
  }

  m = (halir_matrix*)calloc(1, sizeof(halir_matrix));
  if (m == NULL) {
    return NULL;
  }

  m->data = (halir_num*)calloc(rows * cols, sizeof(halir_num));
  if (m->data == NULL) {
    free(m);
    return NULL;
  }

  m->rows = rows;
  m->cols = cols;
  return m;
}

void halir_matrix_free(halir_matrix *m)
{
  if (m == NULL) {
    return;
  }
  free(m->data);
  free(m);
}

/* ----------------------------------------------------------------------- */
/* Grid lifecycle and resampling                                           */
/* ----------------------------------------------------------------------- */

halir_cls_grid *halir_cls_grid_create(size_t n)
{
  halir_cls_grid *grid;

  if (n == 0) {
    return NULL;
  }

  grid = (halir_cls_grid*)calloc(1, sizeof(halir_cls_grid));
  if (grid == NULL) {
    return NULL;
  }

  grid->wavenum = (halir_num*)calloc(n, sizeof(halir_num));
  if (grid->wavenum == NULL) {
    free(grid);
    return NULL;
  }

  grid->n = n;
  return grid;
}

void halir_cls_grid_free(halir_cls_grid *grid)
{
  if (grid == NULL) {
    return;
  }
  free(grid->wavenum);
  free(grid);
}

halir_cls_grid *halir_cls_build_grid(const halir_cls_window *win, size_t nwin, halir_num step)
{
  halir_cls_grid *grid;
  size_t total = 0;
  size_t idx = 0;

  if ((win == NULL) || (nwin == 0) || (!isfinite(step)) || (step <= 0.0)) {
    fprintf(stderr, "halir_cls_build_grid: invalid arguments\n");
    return NULL;
  }

  for (size_t w = 0; w < nwin; w++) {
    if ((!isfinite(win[w].low)) || (!isfinite(win[w].high)) || (win[w].high <= win[w].low)) {
      fprintf(stderr, "halir_cls_build_grid: invalid window %zu\n", w);
      return NULL;
    }
    total += (size_t)floor((win[w].high - win[w].low) / step + 1e-9) + 1;
  }

  grid = halir_cls_grid_create(total);
  if (grid == NULL) {
    fprintf(stderr, "halir_cls_build_grid: allocation failed\n");
    return NULL;
  }

  for (size_t w = 0; w < nwin; w++) {
    size_t npts = (size_t)floor((win[w].high - win[w].low) / step + 1e-9) + 1;
    for (size_t k = 0; k < npts; k++) {
      grid->wavenum[idx++] = win[w].low + (halir_num)k * step;
    }
  }

  return grid;
}

/* Binary search: largest index i with src[i] <= val, or SIZE_MAX if none. */
static size_t lower_bracket(const halir_num *src, size_t n, halir_num val)
{
  size_t lo = 0;
  size_t hi = n;            /* search in [lo, hi) for first src[i] > val */

  while (lo < hi) {
    size_t mid = lo + (hi - lo) / 2;
    if (src[mid] <= val) {
      lo = mid + 1;
    } else {
      hi = mid;
    }
  }
  if (lo == 0) {
    return (size_t)-1;      /* val < src[0] */
  }
  return lo - 1;
}

static int validate_strictly_increasing(const halir_num *src, const halir_num *data, size_t n)
{
  if ((src == NULL) || (data == NULL) || (n == 0)) {
    return 1;
  }

  for (size_t i = 0; i < n; i++) {
    if ((!isfinite(src[i])) || (!isfinite(data[i]))) {
      fprintf(stderr, "halir_cls_resample: source contains non-finite values\n");
      return 1;
    }
    if ((i > 0) && (src[i] <= src[i - 1])) {
      fprintf(stderr, "halir_cls_resample: source wavenumbers must be strictly increasing\n");
      return 1;
    }
  }

  return 0;
}

static halir_cls_calibration *halir_cls_calibration_create(size_t nc,
                                                           size_t nw,
                                                           const halir_cls_grid *grid)
{
  halir_cls_calibration *cal;

  cal = (halir_cls_calibration*)calloc(1, sizeof(halir_cls_calibration));
  if (cal == NULL) {
    return NULL;
  }

  cal->ncomp = nc;
  cal->grid.wavenum = (halir_num*)calloc(nw, sizeof(halir_num));
  cal->K = halir_matrix_create(nc, nw);
  if ((cal->grid.wavenum == NULL) || (cal->K == NULL)) {
    halir_cls_calibration_free(cal);
    return NULL;
  }

  cal->grid.n = nw;
  memcpy(cal->grid.wavenum, grid->wavenum, nw * sizeof(halir_num));
  return cal;
}

static int qr_has_full_column_rank(const gsl_matrix *QR, size_t rows, size_t cols)
{
  double max_diag = 0.0;
  double tol;

  if ((QR == NULL) || (rows == 0) || (cols == 0)) {
    return 0;
  }

  for (size_t i = 0; i < cols; i++) {
    double diag = fabs(gsl_matrix_get(QR, i, i));
    if (diag > max_diag) {
      max_diag = diag;
    }
  }

  if (max_diag == 0.0) {
    return 0;
  }

  tol = DBL_EPSILON * (double)((rows > cols) ? rows : cols) * max_diag * 16.0;
  for (size_t i = 0; i < cols; i++) {
    if (fabs(gsl_matrix_get(QR, i, i)) <= tol) {
      return 0;
    }
  }

  return 1;
}

static halir_cls_calibration *halir_cls_calibrate_weighted_impl(const halir_matrix *A,
                                                                const halir_matrix *C,
                                                                const halir_matrix *W,
                                                                const halir_cls_grid *grid)
{
  halir_cls_calibration *cal = NULL;
  gsl_matrix *QR = NULL;
  gsl_vector *tau = NULL;
  gsl_vector *b = NULL;
  gsl_vector *x = NULL;
  gsl_vector *residual = NULL;
  gsl_error_handler_t *old_handler;
  size_t ns, nc, nw;

  if ((A == NULL) || (C == NULL) || (grid == NULL) ||
      (A->data == NULL) || (C->data == NULL) || (grid->wavenum == NULL)) {
    fprintf(stderr, "halir_cls_calibrate: invalid arguments\n");
    return NULL;
  }

  ns = A->rows;
  nc = C->cols;
  nw = A->cols;

  if (C->rows != ns) {
    fprintf(stderr, "halir_cls_calibrate: A rows (%zu) != C rows (%zu)\n", ns, C->rows);
    return NULL;
  }
  if (grid->n != nw) {
    fprintf(stderr, "halir_cls_calibrate: grid size (%zu) != A cols (%zu)\n", grid->n, nw);
    return NULL;
  }
  if (ns < nc) {
    fprintf(stderr, "halir_cls_calibrate: need samples (%zu) >= components (%zu)\n", ns, nc);
    return NULL;
  }
  if (W != NULL) {
    if ((W->data == NULL) || (W->rows != ns) || (W->cols != nw)) {
      fprintf(stderr, "halir_cls_calibrate_weighted: weights must match A dimensions\n");
      return NULL;
    }
  }

  cal = halir_cls_calibration_create(nc, nw, grid);
  if (cal == NULL) {
    fprintf(stderr, "halir_cls_calibrate: allocation failed\n");
    return NULL;
  }

  QR = gsl_matrix_alloc(ns, nc);
  tau = gsl_vector_alloc(nc);
  b = gsl_vector_alloc(ns);
  x = gsl_vector_alloc(nc);
  residual = gsl_vector_alloc(ns);
  if ((QR == NULL) || (tau == NULL) || (b == NULL) || (x == NULL) || (residual == NULL)) {
    fprintf(stderr, "halir_cls_calibrate: GSL allocation failed\n");
    goto fail;
  }

  old_handler = gsl_set_error_handler_off();

  if (W == NULL) {
    gsl_matrix_const_view cv = gsl_matrix_const_view_array(C->data, ns, nc);
    if (gsl_matrix_memcpy(QR, &cv.matrix) != 0) {
      fprintf(stderr, "halir_cls_calibrate: matrix copy failed\n");
      gsl_set_error_handler(old_handler);
      goto fail;
    }
    if (gsl_linalg_QR_decomp(QR, tau) != 0) {
      gsl_set_error_handler(old_handler);
      fprintf(stderr, "halir_cls_calibrate: QR decomposition failed\n");
      goto fail;
    }
    if (!qr_has_full_column_rank(QR, ns, nc)) {
      gsl_set_error_handler(old_handler);
      fprintf(stderr, "halir_cls_calibrate: calibration design is rank deficient\n");
      goto fail;
    }
  }

  for (size_t j = 0; j < nw; j++) {
    if (W != NULL) {
      for (size_t i = 0; i < ns; i++) {
        halir_num weight = W->data[i * nw + j];
        halir_num scale;

        if ((!isfinite(weight)) || (weight < 0.0)) {
          gsl_set_error_handler(old_handler);
          fprintf(stderr, "halir_cls_calibrate_weighted: invalid weight at sample %zu, frequency %zu\n", i, j);
          goto fail;
        }

        scale = sqrt(weight);
        gsl_vector_set(b, i, scale * A->data[i * nw + j]);
        for (size_t r = 0; r < nc; r++) {
          gsl_matrix_set(QR, i, r, scale * C->data[i * nc + r]);
        }
      }
      if (gsl_linalg_QR_decomp(QR, tau) != 0) {
        gsl_set_error_handler(old_handler);
        fprintf(stderr, "halir_cls_calibrate_weighted: QR decomposition failed (column %zu)\n", j);
        goto fail;
      }
      if (!qr_has_full_column_rank(QR, ns, nc)) {
        gsl_set_error_handler(old_handler);
        fprintf(stderr, "halir_cls_calibrate_weighted: calibration design is rank deficient at column %zu\n", j);
        goto fail;
      }
    } else {
      for (size_t i = 0; i < ns; i++) {
        gsl_vector_set(b, i, A->data[i * nw + j]);
      }
    }

    if (gsl_linalg_QR_lssolve(QR, tau, b, x, residual) != 0) {
      gsl_set_error_handler(old_handler);
      fprintf(stderr, "halir_cls_calibrate: least-squares solve failed (column %zu)\n", j);
      goto fail;
    }
    for (size_t i = 0; i < nc; i++) {
      cal->K->data[i * nw + j] = gsl_vector_get(x, i);
    }
  }

  gsl_set_error_handler(old_handler);
  gsl_matrix_free(QR);
  gsl_vector_free(tau);
  gsl_vector_free(b);
  gsl_vector_free(x);
  gsl_vector_free(residual);
  return cal;

fail:
  gsl_matrix_free(QR);
  gsl_vector_free(tau);
  gsl_vector_free(b);
  gsl_vector_free(x);
  gsl_vector_free(residual);
  halir_cls_calibration_free(cal);
  return NULL;
}

int halir_cls_resample(const halir_num *src_wavenum, const halir_num *src_data,
                       size_t src_n, const halir_cls_grid *grid, halir_num *out)
{
  if ((src_wavenum == NULL) || (src_data == NULL) || (src_n == 0) ||
      (grid == NULL) || (grid->wavenum == NULL) || (out == NULL)) {
    fprintf(stderr, "halir_cls_resample: invalid arguments\n");
    return 1;
  }

  if (validate_strictly_increasing(src_wavenum, src_data, src_n) != 0) {
    return 1;
  }

  for (size_t i = 0; i < grid->n; i++) {
    if (!isfinite(grid->wavenum[i])) {
      fprintf(stderr, "halir_cls_resample: grid contains non-finite values\n");
      return 1;
    }
  }

  for (size_t i = 0; i < grid->n; i++) {
    halir_num x = grid->wavenum[i];

    if ((x < src_wavenum[0]) || (x > src_wavenum[src_n - 1])) {
      out[i] = 0.0;
      continue;
    }

    size_t lo = lower_bracket(src_wavenum, src_n, x);
    if (lo == (size_t)-1) {
      out[i] = 0.0;
      continue;
    }
    if (lo >= src_n - 1) {
      /* x == last sample */
      out[i] = src_data[src_n - 1];
      continue;
    }

    halir_num x0 = src_wavenum[lo];
    halir_num x1 = src_wavenum[lo + 1];
    halir_num dx = x1 - x0;
    if (dx <= 0.0) {
      out[i] = src_data[lo];
    } else {
      halir_num t = (x - x0) / dx;
      out[i] = src_data[lo] + t * (src_data[lo + 1] - src_data[lo]);
    }
  }

  return 0;
}

/* ----------------------------------------------------------------------- */
/* Calibration set generation                                              */
/* ----------------------------------------------------------------------- */

halir_matrix *halir_cls_component_responses(const halir_result *res,
                                            const halir_cls_grid *grid)
{
  halir_matrix *R;

  if ((res == NULL) || (res->spectra == NULL) || (res->nspectra == 0) ||
      (grid == NULL) || (grid->wavenum == NULL) || (grid->n == 0)) {
    fprintf(stderr, "halir_cls_component_responses: invalid arguments\n");
    return NULL;
  }

  R = halir_matrix_create(res->nspectra, grid->n);
  if (R == NULL) {
    fprintf(stderr, "halir_cls_component_responses: allocation failed\n");
    return NULL;
  }

  for (size_t s = 0; s < res->nspectra; s++) {
    const halir_spectra *sp = &res->spectra[s];
    halir_num *row = &R->data[s * R->cols];

    if ((sp->wavenum == NULL) || (sp->data == NULL) || (sp->ndatapnts == 0)) {
      fprintf(stderr, "halir_cls_component_responses: empty spectrum %zu\n", s);
      halir_matrix_free(R);
      return NULL;
    }

    if (halir_cls_resample(sp->wavenum, sp->data, sp->ndatapnts, grid, row) != 0) {
      halir_matrix_free(R);
      return NULL;
    }

    halir_num conc = sp->composition.vmr;
    if (isfinite(conc) && (conc > 0.0)) {
      for (size_t k = 0; k < grid->n; k++) {
        row[k] /= conc;
      }
    }
  }

  return R;
}

halir_matrix *halir_cls_synthesize_A(const halir_matrix *C, const halir_matrix *R)
{
  halir_matrix *A;

  if ((C == NULL) || (R == NULL) || (C->data == NULL) || (R->data == NULL)) {
    fprintf(stderr, "halir_cls_synthesize_A: invalid arguments\n");
    return NULL;
  }
  if (C->cols != R->rows) {
    fprintf(stderr, "halir_cls_synthesize_A: dimension mismatch (C cols %zu != R rows %zu)\n",
            C->cols, R->rows);
    return NULL;
  }

  A = halir_matrix_create(C->rows, R->cols);
  if (A == NULL) {
    fprintf(stderr, "halir_cls_synthesize_A: allocation failed\n");
    return NULL;
  }

  {
    gsl_matrix_const_view cv = gsl_matrix_const_view_array(C->data, C->rows, C->cols);
    gsl_matrix_const_view rv = gsl_matrix_const_view_array(R->data, R->rows, R->cols);
    gsl_matrix_view av = gsl_matrix_view_array(A->data, A->rows, A->cols);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &cv.matrix, &rv.matrix, 0.0, &av.matrix);
  }

  return A;
}

halir_matrix *halir_cls_matrix_scale_rows(const halir_matrix *M,
                                          const halir_num *row_scale)
{
  halir_matrix *scaled;

  if ((M == NULL) || (row_scale == NULL) || (M->data == NULL) || (M->rows == 0) || (M->cols == 0)) {
    fprintf(stderr, "halir_cls_matrix_scale_rows: invalid arguments\n");
    return NULL;
  }

  scaled = halir_matrix_create(M->rows, M->cols);
  if (scaled == NULL) {
    fprintf(stderr, "halir_cls_matrix_scale_rows: allocation failed\n");
    return NULL;
  }

  for (size_t i = 0; i < M->rows; i++) {
    halir_num factor = row_scale[i];
    if (!isfinite(factor)) {
      fprintf(stderr, "halir_cls_matrix_scale_rows: non-finite scale factor at row %zu\n", i);
      halir_matrix_free(scaled);
      return NULL;
    }
    for (size_t j = 0; j < M->cols; j++) {
      scaled->data[i * M->cols + j] = factor * M->data[i * M->cols + j];
    }
  }

  return scaled;
}

halir_matrix *halir_cls_design_augment_column(const halir_matrix *C,
                                              const halir_num *column)
{
  halir_matrix *augmented;

  if ((C == NULL) || (column == NULL) || (C->data == NULL) || (C->rows == 0) || (C->cols == 0)) {
    fprintf(stderr, "halir_cls_design_augment_column: invalid arguments\n");
    return NULL;
  }

  augmented = halir_matrix_create(C->rows, C->cols + 1);
  if (augmented == NULL) {
    fprintf(stderr, "halir_cls_design_augment_column: allocation failed\n");
    return NULL;
  }

  for (size_t i = 0; i < C->rows; i++) {
    if (!isfinite(column[i])) {
      fprintf(stderr, "halir_cls_design_augment_column: non-finite column value at row %zu\n", i);
      halir_matrix_free(augmented);
      return NULL;
    }

    memcpy(&augmented->data[i * augmented->cols],
           &C->data[i * C->cols],
           C->cols * sizeof(halir_num));
    augmented->data[i * augmented->cols + C->cols] = column[i];
  }

  return augmented;
}

halir_matrix *halir_cls_design_augment_inverse_pathlength(const halir_matrix *C,
                                                          const halir_num *pathlength)
{
  halir_num *inverse_pathlength;
  halir_matrix *augmented;

  if ((C == NULL) || (pathlength == NULL) || (C->data == NULL) || (C->rows == 0) || (C->cols == 0)) {
    fprintf(stderr, "halir_cls_design_augment_inverse_pathlength: invalid arguments\n");
    return NULL;
  }

  inverse_pathlength = (halir_num*)calloc(C->rows, sizeof(halir_num));
  if (inverse_pathlength == NULL) {
    fprintf(stderr, "halir_cls_design_augment_inverse_pathlength: allocation failed\n");
    return NULL;
  }

  for (size_t i = 0; i < C->rows; i++) {
    if ((!isfinite(pathlength[i])) || (pathlength[i] <= 0.0)) {
      fprintf(stderr, "halir_cls_design_augment_inverse_pathlength: invalid pathlength at row %zu\n", i);
      free(inverse_pathlength);
      return NULL;
    }
    inverse_pathlength[i] = 1.0 / pathlength[i];
  }

  augmented = halir_cls_design_augment_column(C, inverse_pathlength);
  free(inverse_pathlength);
  return augmented;
}

halir_matrix *halir_cls_design_identity(size_t ncomp)
{
  halir_matrix *C = halir_matrix_create(ncomp, ncomp);
  if (C == NULL) {
    return NULL;
  }
  for (size_t i = 0; i < ncomp; i++) {
    C->data[i * ncomp + i] = 1.0;
  }
  return C;
}

/* ----------------------------------------------------------------------- */
/* Calibration and prediction                                              */
/* ----------------------------------------------------------------------- */

halir_cls_calibration *halir_cls_calibrate(const halir_matrix *A,
                                           const halir_matrix *C,
                                           const halir_cls_grid *grid)
{
  return halir_cls_calibrate_weighted_impl(A, C, NULL, grid);
}

halir_cls_calibration *halir_cls_calibrate_weighted(const halir_matrix *A,
                                                    const halir_matrix *C,
                                                    const halir_matrix *W,
                                                    const halir_cls_grid *grid)
{
  if (W == NULL) {
    fprintf(stderr, "halir_cls_calibrate_weighted: invalid arguments\n");
    return NULL;
  }

  return halir_cls_calibrate_weighted_impl(A, C, W, grid);
}

void halir_cls_calibration_free(halir_cls_calibration *cal)
{
  if (cal == NULL) {
    return;
  }
  free(cal->grid.wavenum);
  halir_matrix_free(cal->K);
  free(cal);
}

halir_cls_prediction *halir_cls_predict(const halir_cls_calibration *cal,
                                        const halir_spectra *sample)
{
  halir_cls_prediction *pred = NULL;
  halir_num *sample_grid = NULL;
  gsl_matrix *QRt = NULL;
  gsl_vector *tau = NULL;
  gsl_vector *b = NULL;
  gsl_vector *x = NULL;
  gsl_vector *residual = NULL;
  gsl_error_handler_t *old_handler;
  size_t nc, nw;

  if ((cal == NULL) || (cal->K == NULL) || (cal->grid.wavenum == NULL) ||
      (sample == NULL) || (sample->wavenum == NULL) || (sample->data == NULL) ||
      (sample->ndatapnts == 0)) {
    fprintf(stderr, "halir_cls_predict: invalid arguments\n");
    return NULL;
  }

  nc = cal->K->rows;
  nw = cal->K->cols;

  if (nw < nc) {
    fprintf(stderr, "halir_cls_predict: need grid points (%zu) >= components (%zu)\n", nw, nc);
    return NULL;
  }

  sample_grid = (halir_num*)calloc(nw, sizeof(halir_num));
  if (sample_grid == NULL) {
    return NULL;
  }
  if (halir_cls_resample(sample->wavenum, sample->data, sample->ndatapnts,
                         &cal->grid, sample_grid) != 0) {
    free(sample_grid);
    return NULL;
  }

  pred = (halir_cls_prediction*)calloc(1, sizeof(halir_cls_prediction));
  if (pred == NULL) {
    free(sample_grid);
    return NULL;
  }
  pred->ncomp = nc;
  pred->conc = (halir_num*)calloc(nc, sizeof(halir_num));
  pred->fitted.wavenum = (halir_num*)calloc(nw, sizeof(halir_num));
  pred->fitted.data = (halir_num*)calloc(nw, sizeof(halir_num));
  pred->residual.wavenum = (halir_num*)calloc(nw, sizeof(halir_num));
  pred->residual.data = (halir_num*)calloc(nw, sizeof(halir_num));
  if ((pred->conc == NULL) || (pred->fitted.wavenum == NULL) || (pred->fitted.data == NULL) ||
      (pred->residual.wavenum == NULL) || (pred->residual.data == NULL)) {
    fprintf(stderr, "halir_cls_predict: allocation failed\n");
    goto fail;
  }
  pred->fitted.ndatapnts = nw;
  pred->residual.ndatapnts = nw;

  /* Solve K^T * conc = sample for conc (least squares, K^T is nw x nc). */
  QRt = gsl_matrix_alloc(nw, nc);
  tau = gsl_vector_alloc(nc);
  b = gsl_vector_alloc(nw);
  x = gsl_vector_alloc(nc);
  residual = gsl_vector_alloc(nw);
  if ((QRt == NULL) || (tau == NULL) || (b == NULL) || (x == NULL) || (residual == NULL)) {
    fprintf(stderr, "halir_cls_predict: GSL allocation failed\n");
    goto fail;
  }

  for (size_t i = 0; i < nw; i++) {
    for (size_t r = 0; r < nc; r++) {
      gsl_matrix_set(QRt, i, r, cal->K->data[r * nw + i]); /* (K^T)[i,r] = K[r,i] */
    }
    gsl_vector_set(b, i, sample_grid[i]);
  }

  old_handler = gsl_set_error_handler_off();
  if (gsl_linalg_QR_decomp(QRt, tau) != 0) {
    gsl_set_error_handler(old_handler);
    fprintf(stderr, "halir_cls_predict: QR decomposition failed\n");
    goto fail;
  }
  if (!qr_has_full_column_rank(QRt, nw, nc)) {
    gsl_set_error_handler(old_handler);
    fprintf(stderr, "halir_cls_predict: calibration matrix is rank deficient\n");
    goto fail;
  }
  if (gsl_linalg_QR_lssolve(QRt, tau, b, x, residual) != 0) {
    gsl_set_error_handler(old_handler);
    fprintf(stderr, "halir_cls_predict: least-squares solve failed\n");
    goto fail;
  }
  gsl_set_error_handler(old_handler);

  for (size_t r = 0; r < nc; r++) {
    pred->conc[r] = gsl_vector_get(x, r);
  }

  /* fitted = conc * K ; residual = sample - fitted */
  for (size_t k = 0; k < nw; k++) {
    halir_num acc = 0.0;
    for (size_t r = 0; r < nc; r++) {
      acc += pred->conc[r] * cal->K->data[r * nw + k];
    }
    pred->fitted.wavenum[k] = cal->grid.wavenum[k];
    pred->fitted.data[k] = acc;
    pred->residual.wavenum[k] = cal->grid.wavenum[k];
    pred->residual.data[k] = sample_grid[k] - acc;
  }

  gsl_matrix_free(QRt);
  gsl_vector_free(tau);
  gsl_vector_free(b);
  gsl_vector_free(x);
  gsl_vector_free(residual);
  free(sample_grid);
  return pred;

fail:
  gsl_matrix_free(QRt);
  gsl_vector_free(tau);
  gsl_vector_free(b);
  gsl_vector_free(x);
  gsl_vector_free(residual);
  free(sample_grid);
  halir_cls_prediction_free(pred);
  return NULL;
}

void halir_cls_prediction_free(halir_cls_prediction *pred)
{
  if (pred == NULL) {
    return;
  }
  free(pred->conc);
  free(pred->fitted.wavenum);
  free(pred->fitted.data);
  free(pred->residual.wavenum);
  free(pred->residual.data);
  free(pred);
}
