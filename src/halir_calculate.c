#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_vector_float.h>
#include <gsl/gsl_math.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include <HalIR/halir.h>
#include <HalIR/tips.h>
#include <cerf.h>

#ifdef HALIR_HAVE_OPENMP
#include <omp.h>
#endif

static void
copy_compound_metadata(const halir_compound *src, halir_compound *dst)
{
  if ((src == NULL) || (dst == NULL)) {
    return;
  }

  memset(dst, 0, sizeof(*dst));
  memcpy(dst->molec, src->molec, sizeof(dst->molec));
  memcpy(dst->isotop, src->isotop, sizeof(dst->isotop));
  dst->vmr = src->vmr;
  dst->conc = src->conc;
  dst->concU = src->concU;
  memcpy(dst->prmfile, src->prmfile, sizeof(dst->prmfile));
  dst->hitran_head = src->hitran_head;
  dst->hitran_head.molecs = NULL;
  dst->hitran_head.nisotp = 0;
  dst->hitran_head.ndatapnts = 0;
  dst->hitran_prms = NULL;
}

size_t find_nearest_index(gsl_vector_float *v, float val)
{
  gsl_vector_float *tmp = gsl_vector_float_alloc(v->size);
  gsl_vector_float_memcpy(tmp, v);
  gsl_vector_float_add_constant(tmp, -1*val);
  gsl_vector_float_mul(tmp, tmp);
  size_t idx = gsl_vector_float_min_index(tmp);
  gsl_vector_float_free(tmp);
  return idx;
}

halir_result *halir_calculate_result(halir_simulation_setup *work)
{
  double q296, qT;
  double T;
  double P;
  double vc, S, q, tfac, mu_step, mu_off1;
  double sig_v = 0.1;
  float alphaD_max, alphaD_min, alphaL_max,  alphaL_min;
  halir_result *result = NULL;

  const double kb_si = 1.3806503e-23; // J/K
  const double kb_erg = 1.3806503e-16; // erg/K
  const double invkb_si = 1./kb_si; // K/J
  const double c_si = 299792458.; // m/s - SI speed of light; used in MM (molecular speed term)
  const double c_erg = 29979245800.; // cm/s - CGS speed of light; used via hc_erg
  const double atmmass_si=1.6605e-27; // kg/amu - atomic mass unit; used in MM
  const double hc_erg = 6.62607015e-27*c_erg; // h*c in erg*cm
  const double c1_erg = hc_erg/kb_erg; // second radiation constant hc/k in cm*K
  const double pi = 3.14159265358979;
  const double ln2 = log(2);
  const double sqrt_pi = sqrt(pi);
  const double sqrt_ln2 = sqrt(ln2);

  if ((work == NULL) || (work->composition == NULL) || (work->composition_length == 0)) {
    fprintf(stderr, "Invalid workspace for calculation\n");
    return NULL;
  }

  T = work->temp;
  P = work->press;

  if ((!isfinite(T)) || (!isfinite(P)) || (T <= 0.0) || (P <= 0.0)) {
    fprintf(stderr, "Invalid temperature or pressure in workspace\n");
    return NULL;
  }
  if ((!isfinite(work->ROI[0])) || (!isfinite(work->ROI[1])) || (work->ROI[1] <= work->ROI[0])) {
    fprintf(stderr, "Invalid ROI range in workspace\n");
    return NULL;
  }

  result = halir_result_create(work, work->composition_length);
  if (result == NULL) {
    fprintf(stderr, "Could not allocate result container\n");
    return NULL;
  }

  for (size_t comp = 0; comp < work->composition_length; comp++) {
    halir_HitranHead *head = &work->composition[comp]->hitran_head;
    halir_HitranLine *prms = work->composition[comp]->hitran_prms;
    gsl_vector_float *v0 = NULL;
    gsl_vector_float *p_S = NULL;
    gsl_vector_float *a_B = NULL;
    gsl_vector_float *MM = NULL;
    gsl_vector_float *taB = NULL;
    gsl_vector_float *alphaD = NULL;
    gsl_vector_float *alphaL = NULL;
    gsl_vector_float *y = NULL;
    double *line_vc = NULL;
    double *line_S = NULL;
    double *line_scale = NULL;
    double *line_aL = NULL;
    size_t *line_start = NULL;
    size_t *line_end = NULL;
    double *thread_buf = NULL;
    size_t mu_size = 0;

    if ((prms == NULL) || (head->ndatapnts <= 0)) {
      fprintf(stderr, "Missing spectral parameters for composition index %zu\n", comp);
      goto calc_error;
    }

    copy_compound_metadata(work->composition[comp], &result->spectra[comp].composition);

    q = work->composition[comp]->vmr;
    tfac = sqrt(2*ln2*kb_si*T);
    // Allocate vectors used in the calculations
    v0 = gsl_vector_float_alloc(head->ndatapnts);
    p_S = gsl_vector_float_alloc(head->ndatapnts);
    a_B = gsl_vector_float_alloc(head->ndatapnts);
    MM = gsl_vector_float_alloc(head->ndatapnts);
    taB = gsl_vector_float_alloc(head->ndatapnts);
    alphaD = gsl_vector_float_alloc(head->ndatapnts);
    alphaL = gsl_vector_float_alloc(head->ndatapnts);
    if ((v0 == NULL) || (p_S == NULL) || (a_B == NULL) || (MM == NULL) ||
        (taB == NULL) || (alphaD == NULL) || (alphaL == NULL)) {
      fprintf(stderr, "Failed to allocate temporary vectors for composition index %zu\n", comp);
      goto comp_error;
    }

    // Populate the vectors with numbers
    for (int i=0; i < head->ndatapnts; i++) {
      gsl_vector_float_set(v0, i, prms[i].trans_mu);
      gsl_vector_float_set(p_S , i, prms[i].pressure_S);
      gsl_vector_float_set(a_B , i, prms[i].air_B );
      gsl_vector_float_set(alphaL , i, prms[i].self_B );
      gsl_vector_float_set(MM  , i, sqrtf(c_si*c_si*atmmass_si*prms[i].molecMass));
      gsl_vector_float_set(taB , i, powf( 296/T, prms[i].temp_air_B) );
    }
    // Correct for pressure shift
    gsl_vector_float_axpby(P, p_S, 1., v0); // P should be P/P_0 but internali alwas use atm pressure units then P_0 = 1
    // Calculate alphaD values
    gsl_vector_float_memcpy(alphaD, v0);
    gsl_vector_float_scale(alphaD, tfac);
    /*gsl_vector_float_scale(MM, c);*/
    gsl_vector_float_div(alphaD, MM);
    // Calculate alphaL values
    gsl_vector_float_axpby((1-q), a_B, q, alphaL);
    gsl_vector_float_mul(alphaL, taB);
    gsl_vector_float_scale(alphaL, P);
    // Get min/max of alphaL and alphaD
    gsl_vector_float_minmax(alphaD, &alphaD_min, &alphaD_max);
    gsl_vector_float_minmax(alphaL, &alphaL_min, &alphaL_max);

    mu_step = sig_v * (alphaD_min + alphaL_min);
    if ((!isfinite(mu_step)) || (mu_step <= 0.0)) {
      fprintf(stderr, "Invalid line widths produced non-positive sampling step\n");
      goto comp_error;
    }
    mu_off1 = ceil(50*GSL_MAX(alphaD_max, alphaL_max));
    if ((!isfinite(mu_off1)) || (mu_off1 < 0.0)) {
      fprintf(stderr, "Invalid line widths produced invalid sampling offset\n");
      goto comp_error;
    }

    /*printf("mu_step: %f mu_off1: %f\n", mu_step, mu_off1);*/
    mu_size = (size_t)ceil(((work->ROI[1]-work->ROI[0])+2*mu_off1)/mu_step);
    if (mu_size < 2) {
      mu_size = 2;
    }
    mu_step = ((work->ROI[1]-work->ROI[0])+2*mu_off1)/(mu_size-1);
    if ((!isfinite(mu_step)) || (mu_step <= 0.0)) {
      fprintf(stderr, "Invalid final sampling step\n");
      goto comp_error;
    }
    double start_value = work->ROI[0]-mu_off1;
    size_t off1 = (size_t)ceil(mu_off1/mu_step);
    /*printf("mu_step: %f mu_size: %ld\n", mu_step, mu_size);*/

    y = gsl_vector_float_calloc(mu_size);
    if (y == NULL) {
      fprintf(stderr, "Failed to allocate output buffer for composition index %zu\n", comp);
      goto comp_error;
    }
    // Print some information (possible in later callback?)
    /*printf("alphaL_min: %f alphaL_max %f\n", alphaL_min, alphaL_max);*/
    /*printf("alphaD_min: %f alphaD_max %f\n", alphaD_min, alphaD_max);*/

    // Per-line parameter arrays for the line-profile phase
    line_vc = (double*)malloc((size_t)head->ndatapnts * sizeof(double));
    line_S = (double*)malloc((size_t)head->ndatapnts * sizeof(double));
    line_scale = (double*)malloc((size_t)head->ndatapnts * sizeof(double));
    line_aL = (double*)malloc((size_t)head->ndatapnts * sizeof(double));
    line_start = (size_t*)malloc((size_t)head->ndatapnts * sizeof(size_t));
    line_end = (size_t*)malloc((size_t)head->ndatapnts * sizeof(size_t));
    if ((line_vc == NULL) || (line_S == NULL) || (line_scale == NULL) ||
        (line_aL == NULL) || (line_start == NULL) || (line_end == NULL)) {
      fprintf(stderr, "Failed to allocate per-line arrays for composition index %zu\n", comp);
      goto comp_error;
    }

    // Serial analysis phase: compute line strengths (incl. Fortran TIPS),
    // widths, amplitude prefactor and grid windows for every line. This
    // keeps the non-reentrant tips_2020() out of the parallel region.
    for (int mm = 0; mm < head->ndatapnts; mm++) {
      vc = gsl_vector_float_get(v0, mm);
      // analytic nearest index on the uniform grid (replaces find_nearest_index)
      double fidx = (vc - start_value) / mu_step;
      long lidx = lround(fidx);
      if (lidx < 0) {
        lidx = 0;
      } else if ((size_t)lidx >= mu_size) {
        lidx = (long)mu_size - 1;
      }
      size_t idx = (size_t)lidx;

      q296 = tips_2020(prms[mm].molec_num, prms[mm].isotp_num, 296);
      qT = tips_2020(prms[mm].molec_num, prms[mm].isotp_num, T);

      double a_D = (double)gsl_vector_float_get(alphaD, mm);
      double a_L = (double)gsl_vector_float_get(alphaL, mm);

      if (qT == 0.0) {
        line_start[mm] = 0;
        line_end[mm] = 0;
        continue;
      }
      S = prms[mm].line_I * q296/qT * exp(c1_erg*prms[mm].low_state_en/T)/exp(c1_erg*prms[mm].low_state_en/296)*((1-exp(-c1_erg*vc/T))/(1-exp(-c1_erg*vc/296)));

      if (a_D <= 0.0f || !isfinite(a_D) || !isfinite(a_L) || !isfinite(S)) {
        line_start[mm] = 0;
        line_end[mm] = 0;
        continue;
      }

      size_t start_i = (idx > off1) ? (idx - off1) : 0;
      size_t end_i = idx + off1;
      if (end_i > mu_size) {
        end_i = mu_size;
      }

      line_vc[mm] = vc;
      // amplitude prefactor: keep identical left-to-right association as the
      // original loop so the per-point value is bit-for-bit unchanged
      line_S[mm] = sqrt_ln2*S/(sqrt_pi*a_D) * q * work->pathL/100 * P*101325/1e4 * invkb_si/T;
      // profile scale sqrt_ln2/a_D hoisted out of the per-point inner loop
      line_scale[mm] = sqrt_ln2/a_D;
      line_aL[mm] = sqrt_ln2*a_L/a_D;
      line_start[mm] = start_i;
      line_end[mm] = end_i;
    }

    // Parallel line-profile phase: each thread accumulates its lines into a
    // private slice of thread_buf; slices are merged after the region. The
    // component loop itself stays serial (compute load is in the profile).
    {
      int nthreads = 1;
#ifdef HALIR_HAVE_OPENMP
      nthreads = omp_get_max_threads();
      if (nthreads < 1) {
        nthreads = 1;
      }
#endif
      thread_buf = (double*)calloc((size_t)nthreads * mu_size, sizeof(double));
      if (thread_buf == NULL) {
        fprintf(stderr, "Failed to allocate thread buffers for composition index %zu\n", comp);
        goto comp_error;
      }

#ifdef HALIR_HAVE_OPENMP
      #pragma omp parallel
      {
        int tid = omp_get_thread_num();
        double *buf = thread_buf + (size_t)tid * mu_size;
        #pragma omp for schedule(static)
        for (int mm = 0; mm < head->ndatapnts; mm++) {
          double vc_m = line_vc[mm];
          double amp = line_S[mm];
          double scale = line_scale[mm];
          double yterm = line_aL[mm];
          for (size_t i = line_start[mm]; i < line_end[mm]; i++) {
            float mu_i = (float)(start_value + (double)i*mu_step);
            buf[i] += amp * re_w_of_z(scale*((double)mu_i-vc_m), yterm);
          }
        }
      }
#else
      for (int mm = 0; mm < head->ndatapnts; mm++) {
        double vc_m = line_vc[mm];
        double amp = line_S[mm];
        double scale = line_scale[mm];
        double yterm = line_aL[mm];
        for (size_t i = line_start[mm]; i < line_end[mm]; i++) {
          float mu_i = (float)(start_value + (double)i*mu_step);
          thread_buf[i] += amp * re_w_of_z(scale*((double)mu_i-vc_m), yterm);
        }
      }
#endif

      // Merge per-thread slices into the output buffer
      for (int t = 0; t < nthreads; t++) {
        const double *buf = thread_buf + (size_t)t * mu_size;
        for (size_t i = 0; i < mu_size; i++) {
          y->data[i] += (float)buf[i];
        }
      }
    }

    result->spectra[comp].ndatapnts = mu_size;
    result->spectra[comp].wavenum = (halir_num*)calloc(mu_size, sizeof(halir_num));
    result->spectra[comp].data = (halir_num*)calloc(mu_size, sizeof(halir_num));
    if ((result->spectra[comp].wavenum == NULL) || (result->spectra[comp].data == NULL)) {
      fprintf(stderr, "Failed to allocate result arrays for composition index %zu\n", comp);
      goto comp_error;
    }

    for (size_t i = 0; i < mu_size; i++) {
      result->spectra[comp].wavenum[i] = (halir_num)(float)(start_value + (double)i*mu_step);
      result->spectra[comp].data[i] = (halir_num)y->data[i];
    }

    gsl_vector_float_free(y);
    gsl_vector_float_free(v0);
    gsl_vector_float_free(p_S);
    gsl_vector_float_free(a_B);
    gsl_vector_float_free(alphaL);
    gsl_vector_float_free(alphaD);
    gsl_vector_float_free(MM);
    gsl_vector_float_free(taB);
    free(line_vc);
    free(line_S);
    free(line_scale);
    free(line_aL);
    free(line_start);
    free(line_end);
    free(thread_buf);

    continue;

comp_error:
    gsl_vector_float_free(y);
    gsl_vector_float_free(v0);
    gsl_vector_float_free(p_S);
    gsl_vector_float_free(a_B);
    gsl_vector_float_free(alphaL);
    gsl_vector_float_free(alphaD);
    gsl_vector_float_free(MM);
    gsl_vector_float_free(taB);
    free(line_vc);
    free(line_S);
    free(line_scale);
    free(line_aL);
    free(line_start);
    free(line_end);
    free(thread_buf);
    goto calc_error;
  }

  return result;

calc_error:
  halir_result_free(result);
  return NULL;
}

int halir_test_calc(halir_simulation_setup *work)
{
  halir_result *result;

  result = halir_calculate_result(work);
  if (result == NULL) {
    return 1;
  }

  for (size_t comp = 0; comp < result->nspectra; comp++) {
    for (size_t i = 0; i < result->spectra[comp].ndatapnts; i++) {
      printf("%f %f\n", result->spectra[comp].wavenum[i], result->spectra[comp].data[i]);
    }
  }

  halir_result_free(result);
  return 0;
}
