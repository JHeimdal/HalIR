#ifndef HALIR_CLS_H_
#define HALIR_CLS_H_

#include <stddef.h>
#include <HalIR/halir.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @file halir_cls.h
 * @brief Classical Least Squares (CLS) calibration and prediction for HalIR.
 *
 * Pipeline (see project workflow):
 *   1. Generate synthetic component responses on a common wavenumber grid.
 *   2. Build a calibration set A = C * R from a concentration design C.
 *   3. Calibrate: solve K (Ncomp x Nwave) from A and C.
 *   4. Predict: recover sample concentrations and fitted/residual spectra.
 *
 * Matrices are dense and row-major. All numerical values use halir_num.
 */

/**
 * @brief Dense, row-major matrix of halir_num values.
 */
typedef struct {
  size_t rows;       /**< Number of rows. */
  size_t cols;       /**< Number of columns. */
  halir_num *data;   /**< rows*cols values, row-major. */
} halir_matrix;

/**
 * @brief A single calibration window [low, high] in wavenumber (cm^-1).
 */
typedef struct {
  halir_num low;     /**< Lower bound (cm^-1), inclusive. */
  halir_num high;    /**< Upper bound (cm^-1), inclusive. */
} halir_cls_window;

/**
 * @brief A common wavenumber grid used by the CLS matrices.
 */
typedef struct {
  size_t n;            /**< Number of grid points. */
  halir_num *wavenum;  /**< Monotonically increasing grid (length n). */
} halir_cls_grid;

/**
 * @brief Output of CLS calibration: the K (Ncomp x Nwave) matrix on a grid.
 */
typedef struct {
  halir_cls_grid grid; /**< Wavenumber grid the columns of K refer to. */
  size_t ncomp;        /**< Number of components (rows of K). */
  halir_matrix *K;     /**< Calibration matrix, ncomp x grid.n. */
} halir_cls_calibration;

/**
 * @brief Output of CLS prediction for a single sample spectrum.
 */
typedef struct {
  size_t ncomp;            /**< Number of predicted concentrations. */
  halir_num *conc;         /**< Predicted concentrations (length ncomp). */
  halir_spectra fitted;    /**< Fitted spectrum on the calibration grid. */
  halir_spectra residual;  /**< Residual (sample - fitted) on the grid. */
} halir_cls_prediction;

/* ----------------------------------------------------------------------- */
/* Matrix lifecycle                                                        */
/* ----------------------------------------------------------------------- */

/**
 * @brief Allocate a zero-initialized matrix.
 * @param rows Number of rows.
 * @param cols Number of columns.
 * @return New matrix or NULL on invalid args / allocation failure.
 */
halir_matrix *halir_matrix_create(size_t rows, size_t cols);

/**
 * @brief Free a matrix (NULL-safe).
 * @param m Matrix to free.
 */
void halir_matrix_free(halir_matrix *m);

/** @brief Element accessor (no bounds checking). */
static inline halir_num halir_matrix_get(const halir_matrix *m, size_t r, size_t c)
{
  return m->data[r * m->cols + c];
}

/** @brief Element setter (no bounds checking). */
static inline void halir_matrix_set(halir_matrix *m, size_t r, size_t c, halir_num v)
{
  m->data[r * m->cols + c] = v;
}

/* ----------------------------------------------------------------------- */
/* Grid lifecycle and resampling                                           */
/* ----------------------------------------------------------------------- */

/**
 * @brief Allocate a grid of n points (wavenum zero-initialized).
 * @param n Number of grid points.
 * @return New grid or NULL on failure.
 */
halir_cls_grid *halir_cls_grid_create(size_t n);

/**
 * @brief Free a grid (NULL-safe).
 * @param grid Grid to free.
 */
void halir_cls_grid_free(halir_cls_grid *grid);

/**
 * @brief Build a uniform grid spanning the union of calibration windows.
 *
 * Each window contributes a uniform segment with the given step; the segments
 * are concatenated in window order. Windows must be valid (high > low) and
 * step must be positive and finite.
 *
 * @param win   Array of nwin windows.
 * @param nwin  Number of windows (>= 1).
 * @param step  Grid spacing (cm^-1), > 0.
 * @return New grid or NULL on invalid args / allocation failure.
 */
halir_cls_grid *halir_cls_build_grid(const halir_cls_window *win, size_t nwin, halir_num step);

/**
 * @brief Linearly resample a spectrum onto a target grid.
 *
 * Points outside the source range are set to 0. The source wavenum array must
 * be monotonically increasing. Writes grid->n values into out.
 *
 * @param src_wavenum Source wavenumber array.
 * @param src_data Source data array.
 * @param src_n Number of source points.
 * @param grid Target grid.
 * @param out Output array with at least grid->n elements.
 * @return 0 on success, non-zero on invalid args.
 */
int halir_cls_resample(const halir_num *src_wavenum, const halir_num *src_data,
                       size_t src_n, const halir_cls_grid *grid, halir_num *out);

/* ----------------------------------------------------------------------- */
/* Calibration set generation                                              */
/* ----------------------------------------------------------------------- */

/**
 * @brief Build the component response matrix R (Ncomp x Nwave).
 *
 * Each per-component spectrum in @p res is resampled onto @p grid and divided
 * by that component's concentration (vmr) to give a unit-concentration
 * response. Components with non-positive vmr are left as the raw resampled
 * spectrum (treated as already unit-scaled).
 *
 * @param res Result spectra from halir_calculate_result.
 * @param grid Common CLS wavenumber grid.
 * @return New matrix (res->nspectra x grid->n) or NULL on failure.
 */
halir_matrix *halir_cls_component_responses(const halir_result *res,
                                            const halir_cls_grid *grid);

/**
 * @brief Synthesize calibration spectra A = C * R.
 *
 * @param C  Concentration design, Nsamples x Ncomp.
 * @param R  Component responses, Ncomp x Nwave.
 * @return New matrix Nsamples x Nwave or NULL on invalid args / failure.
 */
halir_matrix *halir_cls_synthesize_A(const halir_matrix *C, const halir_matrix *R);

/**
 * @brief Return a copy of a matrix with each row scaled by a caller-supplied factor.
 *
 * This helper is intended for caller-side preprocessing of calibration or
 * sample matrices, including inverse-pathlength normalization.
 *
 * @param M            Source matrix, Nrows x Ncols.
 * @param row_scale    One multiplicative scale factor per row, length Nrows.
 * @return New matrix Nrows x Ncols or NULL on invalid args / failure.
 */
halir_matrix *halir_cls_matrix_scale_rows(const halir_matrix *M,
                                          const halir_num *row_scale);

/**
 * @brief Append one explicit design column to a concentration design matrix.
 *
 * This helper is intended for caller-controlled augmented CLS designs, such as
 * constant intercept terms or inverse-pathlength columns. The appended column
 * must contain one value per sample (row of @p C).
 *
 * @param C       Base concentration design, Nsamples x Ncomp.
 * @param column  Column values to append, length Nsamples.
 * @return New matrix Nsamples x (Ncomp + 1) or NULL on invalid args / failure.
 */
halir_matrix *halir_cls_design_augment_column(const halir_matrix *C,
                                              const halir_num *column);

/**
 * @brief Append an inverse-pathlength column to a concentration design matrix.
 *
 * Each pathlength value must be finite and strictly positive. The appended
 * column contains 1/pathlength for each sample row and is intended for the
 * explicit intercept/pathlength workflow described in the CLS paper.
 *
 * @param C            Base concentration design, Nsamples x Ncomp.
 * @param pathlength   Pathlength values, length Nsamples.
 * @return New matrix Nsamples x (Ncomp + 1) or NULL on invalid args / failure.
 */
halir_matrix *halir_cls_design_augment_inverse_pathlength(const halir_matrix *C,
                                                          const halir_num *pathlength);

/**
 * @brief Build an identity concentration design (one pure sample per component).
 * @param ncomp Number of components.
 * @return New Ncomp x Ncomp identity matrix or NULL on failure.
 */
halir_matrix *halir_cls_design_identity(size_t ncomp);

/* ----------------------------------------------------------------------- */
/* Calibration and prediction                                              */
/* ----------------------------------------------------------------------- */

/**
 * @brief Solve the CLS calibration matrix K from A and C.
 *
 * Solves the least-squares system A = C * K for K (Ncomp x Nwave) using a QR
 * factorization of C. Requires Nsamples >= Ncomp.
 *
 * @param A     Calibration spectra, Nsamples x Nwave.
 * @param C     Concentrations, Nsamples x Ncomp.
 * @param grid  Wavenumber grid (grid->n must equal A->cols); copied into output.
 * @return New calibration object or NULL on invalid args / failure.
 */
halir_cls_calibration *halir_cls_calibrate(const halir_matrix *A,
                                           const halir_matrix *C,
                                           const halir_cls_grid *grid);

/**
 * @brief Solve the weighted CLS calibration matrix K from A, C, and W.
 *
 * Solves the weighted least-squares system A = C * K for K (Ncomp x Nwave)
 * using per-frequency diagonal weights. The weight matrix W must have the
 * same shape as A, where W[i,j] is the weight applied to sample i at
 * frequency j. Each spectral column is solved independently after scaling the
 * calibration equations by sqrt(W[:,j]). Requires Nsamples >= Ncomp.
 *
 * @param A        Calibration spectra, Nsamples x Nwave.
 * @param C        Concentrations, Nsamples x Ncomp.
 * @param W        Weights, Nsamples x Nwave.
 * @param grid     Wavenumber grid (grid->n must equal A->cols); copied into output.
 * @return New calibration object or NULL on invalid args / failure.
 */
halir_cls_calibration *halir_cls_calibrate_weighted(const halir_matrix *A,
                                                    const halir_matrix *C,
                                                    const halir_matrix *W,
                                                    const halir_cls_grid *grid);

/**
 * @brief Free a calibration object (NULL-safe).
 * @param cal Calibration object to free.
 */
void halir_cls_calibration_free(halir_cls_calibration *cal);

/**
 * @brief Predict concentrations and fitted/residual spectra for a sample.
 *
 * The sample spectrum is resampled onto the calibration grid, then the system
 * sample = conc * K is solved for conc (length Ncomp) using a QR factorization
 * of K^T. Requires Nwave >= Ncomp.
 *
 * @param cal     Calibration produced by halir_cls_calibrate.
 * @param sample  Measured sample spectrum (own grid).
 * @return New prediction object or NULL on invalid args / failure.
 */
halir_cls_prediction *halir_cls_predict(const halir_cls_calibration *cal,
                                        const halir_spectra *sample);

/**
 * @brief Free a prediction object (NULL-safe).
 * @param pred Prediction object to free.
 */
void halir_cls_prediction_free(halir_cls_prediction *pred);

#ifdef __cplusplus
}
#endif

#endif /* HALIR_CLS_H_ */
