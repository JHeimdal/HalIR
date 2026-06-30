#ifndef HALIR_H_
#define HALIR_H_

#include <ctype.h>
#include <stdio.h>
#include <stddef.h>
#include <gsl/gsl_vector.h>
#include <string.h>

#define UNIT_LEN 5
#define PATH_LEN 1024
#define COMMENT_LEN 4096

#define halir_num double

/**
 * Enumerates the available apodization methods for the HALIR algorithm.
 */
typedef enum
{
    /**
     * Boxcar apodization method.
     */
    HALIR_BOXCAR,

    /**
     * Triangle apodization method.
     */
    HALIR_TRIANGLE,

    /**
     * Hanning-Papoulis-Genzel (HAPP) apodization method.
     */
    HALIR_HAPPGENZEL,
} halir_apodization;

typedef enum
{
  HALIR_TRANSMISSION,
  HALIR_ABSORBANCE,
} halir_filetype;

/**
 * Enumerated type representing various units of measurement.
 */
typedef enum
{
    NONE,  /**< Unknown or unspecified unit */
    ATM,   /**< Atmospheres (unit of pressure) */
    MBAR,  /**< Millibars (unit of pressure) */
    BAR,   /**< Bars (unit of pressure) */
    PA,    /**< Pascals (unit of pressure) */
    HPA,   /**< Heighs of mercury (unit of pressure) */
    MMHG,  /**< Millimeters of mercury (unit of pressure) */
    PPM,   /**< Parts per million (unit of concentration) */
    PPT,   /**< Parts per trillion (unit of concentration) */
    PPB,   /**< Parts per billion (unit of concentration) */
    MM,    /**< Millimeters (unit of length) */
    CM,    /**< Centimeters (unit of length) */
    DM,    /**< Decimeters (unit of length) */
    M,     /**< Meters (unit of length) */
    KM,    /**< Kilometers (unit of length) */
    K,     /**< Kelvin (unit of temperature) */
    C,     /**< Celsius (unit of temperature) */
    F,     /**< Fahrenheit (unit of temperature) */
} halir_Units;

// Convert halir_Units to atm, cm, Kelvin
static double halir_Units_factors[] = {
  0.0, //None
  1.0, //ATM,
  9.86923266716013e-4, //MBAR,
  9.86923266716013e-1, //BAR,
  9.86923266716013e-6, //PA,
  9.86923266716013e-4, //HPA,
  1.31578947368421e-3, //MMHG,
  1.0e-6, //PPM,
  1.0e-12, //PPT,
  1.0e-9, //PPB,
  1.0e-1, //MM,
  1.0, //CM,
  1.0e1, //DM,
  1.0e2, //M,
  1.0e5, //KM,
  1.0, //CM,
  1.0, //CM,
  1.0, //CM,
};

static const char *halir_Units_to_str[] = {
  "NONE",
  "atm",
  "mbar",
  "bar",
  "pa",
  "hpa",
  "mmhg",
  "ppm",
  "ppt",
  "ppb",
  "mm",
  "cm",
  "dm",
  "m",
  "km",
  "K",
  "C",
  "F",
};

#define HALIR_H_NUM_UNITS 18

static inline double CtoK(const double *temp) {
  return *temp + 273.15;
}
static inline double FtoK(const double *temp) {
  return (*temp + 459.67) * 5./9.;
}
static inline double halir_Units_to_Hitran(const halir_Units *unit, const double *value, int *read_err) {
  switch (*unit) {
    case ATM:
      *read_err = 0;
      return halir_Units_factors[ATM] * (*value);
      break;
    case BAR:
      *read_err = 0;
      return halir_Units_factors[BAR] * (*value);
      break;
    case MBAR:
      *read_err = 0;
      return halir_Units_factors[MBAR] * (*value);
      break;
    case PA:
      *read_err = 0;
      return halir_Units_factors[PA] * (*value);
      break;
    case HPA:
      *read_err = 0;
      return halir_Units_factors[HPA] * (*value);
      break;
    case MMHG:
      *read_err = 0;
      return halir_Units_factors[MMHG] * (*value);
      break;
    case PPM:
      *read_err = 0;
      return halir_Units_factors[PPM] * (*value);
      break;
    case PPT:
      *read_err = 0;
      return halir_Units_factors[PPT] * (*value);
      break;
    case PPB:
      *read_err = 0;
      return halir_Units_factors[PPB] * (*value);
      break;
    case MM:
      *read_err = 0;
      return halir_Units_factors[MM] * (*value);
      break;
    case CM:
      *read_err = 0;
      return halir_Units_factors[CM] * (*value);
      break;
    case DM:
      *read_err = 0;
      return halir_Units_factors[DM] * (*value);
      break;
    case M:
      *read_err = 0;
      return halir_Units_factors[M] * (*value);
      break;
    case KM:
      *read_err = 0;
      return halir_Units_factors[KM] * (*value);
      break;
    case K:
      *read_err = 0;
      return (*value);
      break;
    case C:
      *read_err = 0;
      return CtoK(value);
      break;
    case F:
      *read_err = 0;
      return FtoK(value);
      break;
    default:
      *read_err = 3;
      return 0;
  }
}

/**
 * @brief Header metadata for a loaded HITRAN parameter file.
 */
typedef struct {
  char molec[6];      /**< Molecule name/identifier. */
  char isoname[12];   /**< Isotopologue label from the prmfile. */
  int nisotp;         /**< Number of isotopologue IDs in @ref molecs. */
  float roi_low;      /**< Lower spectral ROI bound (cm^-1). */
  float roi_high;     /**< Upper spectral ROI bound (cm^-1). */
  int ndatapnts;      /**< Number of line entries in @ref halir_compound::hitran_prms. */
  int *molecs;        /**< Dynamic array of isotopologue identifiers. */
} halir_HitranHead;

/**
 * @brief One spectral line entry from HITRAN parameter data.
 */
typedef struct {
  int    molec_num;        /**< HITRAN molecule number. */
  int    isotp_num;        /**< HITRAN isotopologue number. */
  halir_num trans_mu;      /**< Transition wavenumber (cm^-1). */
  halir_num line_I;        /**< Line intensity. */
  halir_num einstein_A;    /**< Einstein A coefficient. */
  halir_num air_B;         /**< Air-broadened half-width. */
  halir_num self_B;        /**< Self-broadened half-width. */
  halir_num low_state_en;  /**< Lower-state energy. */
  halir_num temp_air_B;    /**< Temperature dependence of air width. */
  halir_num pressure_S;    /**< Pressure shift term. */
  char   u_vib_quant[16];  /**< Upper vibrational quantum numbers. */
  char   l_vib_quant[16];  /**< Lower vibrational quantum numbers. */
  char   u_loc_quant[16];  /**< Upper local quantum numbers. */
  char   l_loc_quant[16];  /**< Lower local quantum numbers. */
  char   err_code[6];      /**< Error code string from source data. */
  char   ref_code[12];     /**< Reference code string from source data. */
  char   line_mix[2];      /**< Line-mixing flag. */
  halir_num u_stat_w;      /**< Upper-state statistical weight. */
  halir_num l_stat_w;      /**< Lower-state statistical weight. */
  halir_num abundance;     /**< Isotopologue abundance factor. */
  halir_num molecMass;     /**< Isotopologue molecular mass. */
} halir_HitranLine;

/**
 * @brief Composition entry describing one gas species in the workspace.
 */
typedef struct {
  char molec[6];                 /**< Molecule name/identifier. */
  char isotop[12];               /**< Isotopologue label. */
  halir_num vmr;                 /**< Volume mixing ratio, expected in [0,1]. */
  halir_num conc;                /**< Concentration value (normalized by API helpers). */
  halir_Units concU;             /**< Unit associated with @ref conc. */
  char prmfile[PATH_LEN];        /**< Optional prmfile path for line data loading. */
  halir_HitranHead hitran_head;  /**< Loaded prmfile header metadata. */
  halir_HitranLine *hitran_prms; /**< Dynamic array of loaded line parameters. */
} halir_compound;

/**
 * @brief Full simulation workspace for one HALIR run.
 */
typedef struct
{
  char pname[PATH_LEN];           /**< Project name. */
  char rootDir[PATH_LEN];         /**< Root directory for project assets. */
  char pcomments[COMMENT_LEN];    /**< Free-form project comments. */
  char **pfiles;                  /**< Dynamic list of project file paths. */
  size_t pfiles_length;           /**< Number of entries in @ref pfiles. */

  halir_num temp;                 /**< Temperature (canonical unit: K). */
  halir_num press;                /**< Pressure (canonical unit: ATM). */
  halir_num pathL;                /**< Path length (canonical unit: CM). */
  halir_num ROI[2];               /**< Spectral region of interest [low, high]. */
  halir_num res;                  /**< Spectral resolution. */
  halir_num fov;                  /**< Field-of-view parameter. */
  halir_Units tempU;              /**< Temperature unit marker (typically K after normalization). */
  halir_Units pressU;             /**< Pressure unit marker (typically ATM after normalization). */
  halir_Units pathLU;             /**< Path length unit marker (typically CM after normalization). */
  halir_apodization apod;         /**< Apodization mode. */
  halir_filetype ftype;           /**< Spectrum output type. */
  char bgfile[PATH_LEN];          /**< Optional background file path. */

  size_t composition_length;      /**< Number of composition entries. */
  halir_compound **composition;   /**< Dynamic array of composition pointers. */
} halir_simulation_setup;

/**
 * @brief Result container for a simulated spectrum.
 *
 * Ownership model:
 * - `wavenum` and `data` are owned by this struct instance.
 * - `composition` stores metadata copied from the workspace entry.
 * - Nested source buffers (`composition.hitran_prms`,
 *   `composition.hitran_head.molecs`) are not deep-copied from the source
 *   workspace by calculation APIs.
 */
typedef struct{
  size_t ndatapnts;          /**< Number of spectral samples. */
  halir_compound composition; /**< Composition metadata used for the result. */
  halir_num *data;           /**< Spectral values array. */
  halir_num *wavenum;        /**< Wavenumber values array. */
} halir_spectra;

/* Backward-compatible alias for historical typo. */
typedef halir_spectra hair_spectra;

/**
 * @brief Result container for a complete HalIR simulation run.
 *
 * Ownership model:
 * - `spectra` array is owned by this struct instance.
 * - Each spectra entry owns its `wavenum` and `data` buffers.
 * - `workspace` is a non-owning pointer for provenance only.
 */
typedef struct{
  halir_simulation_setup *workspace; /**< Pointer to the workspace used for the calculation. */
  size_t nspectra;           /**< Number of spectra in the result. */
  halir_spectra *spectra;    /**< Dynamic array of spectra. */
} halir_result;

/**
 * @defgroup halir_api HALIR Public API
 * @brief Public entry points for workspace construction, validation, and execution.
 *
 * Return codes are documented per function and summarized in
 * @ref md_docs_2error__codes "HALIR Error Code Matrix".
 */

/**
 * @defgroup halir_api_lifecycle Workspace Lifecycle
 * @ingroup halir_api
 * @brief Allocate, release, and validate workspace state.
 */

/**
 * @defgroup halir_api_setup Workspace Setup
 * @ingroup halir_api
 * @brief Configure project metadata and sample environment.
 */

/**
 * @defgroup halir_api_composition Composition Setup
 * @ingroup halir_api
 * @brief Add and configure gas composition entries.
 */

/**
 * @defgroup halir_api_io Data I/O and Execution
 * @ingroup halir_api
 * @brief Parse input, load prmfile data, and run calculation.
 */

/**
 * @defgroup halir_api_utils Utilities
 * @ingroup halir_api
 * @brief Supporting helper routines.
 */

/**
 * @brief Parse HALIR JSON input and build a workspace.
 * @ingroup halir_api_io
 *
 * @param inputFile Path to a JSON file or a JSON string buffer.
 * @return Pointer to a populated workspace on success, or NULL on failure.
 */
halir_simulation_setup* halir_parseJSONinput(const char* const inputFile);

/**
 * @brief Allocate a zero-initialized HALIR workspace.
 * @ingroup halir_api_lifecycle
 *
 * The returned workspace must be released with halir_simulation_setup_free().
 *
 * @return Pointer to a new workspace, or NULL if allocation fails.
 */
halir_simulation_setup* halir_simulation_setup_create(void);

/**
 * @brief Free a workspace and all nested dynamic allocations.
 * @ingroup halir_api_lifecycle
 *
 * This function is NULL-safe.
 *
 * @param work Workspace to release.
 */
void halir_simulation_setup_free(halir_simulation_setup *work);

/**
 * @brief Set project metadata fields in a workspace.
 * @ingroup halir_api_setup
 *
 * @param work Workspace instance.
 * @param pname Project name.
 * @param rootDir Root directory for project assets.
 * @param pcomments Free-form project comments.
 * @return 0 on success.
 * @return 1 when any argument is NULL.
 * @return 2 when bounded string copy fails.
 */
int halir_simulation_setup_set_project(halir_simulation_setup *work, const char *pname, const char *rootDir, const char *pcomments);

/**
 * @brief Append a project file path to the workspace file list.
 * @ingroup halir_api_setup
 *
 * @param work Workspace instance.
 * @param path File path to append.
 * @return 0 on success.
 * @return 1 when arguments are invalid.
 * @return 2 when memory allocation fails.
 */
int halir_simulation_setup_add_project_file(halir_simulation_setup *work, const char *path);

/**
 * @brief Set and normalize sample-environment parameters.
 * @ingroup halir_api_setup
 *
 * Temperature, pressure and path length are converted to canonical internal
 * units (K, ATM, CM) at the API boundary.
 *
 * @param work Workspace instance.
 * @param temp Temperature value.
 * @param tempU Temperature unit.
 * @param press Pressure value.
 * @param pressU Pressure unit.
 * @param pathL Optical path length.
 * @param pathLU Optical path length unit.
 * @param roi_low ROI lower bound in wavenumber.
 * @param roi_high ROI upper bound in wavenumber.
 * @param res Spectral resolution.
 * @param fov Field of view.
 * @param apod Apodization mode.
 * @param ftype Output spectrum type.
 * @param bgfile Optional background file path.
 * @return 0 on success.
 * @return 1 when input arguments fail validation.
 * @return 2 when conversion or string assignment fails.
 */
int halir_simulation_setup_set_sample_env(halir_simulation_setup *work,
                                   double temp, halir_Units tempU,
                                   double press, halir_Units pressU,
                                   double pathL, halir_Units pathLU,
                                   double roi_low, double roi_high,
                                   double res, double fov,
                                   halir_apodization apod,
                                   halir_filetype ftype,
                                   const char *bgfile);

/**
 * @brief Add a composition entry to the workspace.
 * @ingroup halir_api_composition
 *
 * @param work Workspace instance.
 * @param molec Molecule identifier.
 * @param isotop Isotopologue label.
 * @param out_index Output index of the new composition entry.
 * @return 0 on success.
 * @return 1 when arguments are invalid.
 * @return 2 when memory allocation or copy fails.
 */
int halir_simulation_setup_add_composition(halir_simulation_setup *work, const char *molec, const char *isotop, size_t *out_index);

/**
 * @brief Set volume mixing ratio for an existing composition entry.
 * @ingroup halir_api_composition
 *
 * @param work Workspace instance.
 * @param comp_index Composition index.
 * @param vmr Volume mixing ratio in [0, 1].
 * @return 0 on success.
 * @return 1 when workspace/index/vmr validation fails.
 */
int halir_compound_set_vmr(halir_simulation_setup *work, size_t comp_index, double vmr);

/**
 * @brief Set concentration for an existing composition entry.
 * @ingroup halir_api_composition
 *
 * Concentration is normalized using halir_Units_to_Hitran().
 *
 * @param work Workspace instance.
 * @param comp_index Composition index.
 * @param conc Concentration value.
 * @param concU Concentration unit.
 * @return 0 on success.
 * @return 1 when workspace/index/value validation fails.
 * @return 2 when unit conversion fails.
 */
int halir_compound_set_concentration(halir_simulation_setup *work, size_t comp_index, double conc, halir_Units concU);

/**
 * @brief Load HITRAN prmfile data for a composition entry.
 * @ingroup halir_api_io
 *
 * @param work Workspace instance.
 * @param comp_index Composition index.
 * @param prmfile_path Path to prmfile.
 * @return 0 on success.
 * @return 1 when arguments are invalid.
 * @return 2 when file access/read/format loading fails.
 */
int halir_compound_load_prmfile(halir_simulation_setup *work, size_t comp_index, const char *prmfile_path);

/**
 * @brief Validate that a workspace is complete and numerically consistent.
 * @ingroup halir_api_lifecycle
 *
 * Validation checks include required project metadata, sample environment
 * constraints, ROI ordering, and composition-level concentration bounds.
 *
 * @param work Workspace instance.
 * @return 0 when valid.
 * @return 1 when validation fails.
 */
int halir_simulation_setup_validate(halir_simulation_setup *work);

/**
 * @brief Find the nearest index in a float vector for a target value.
 * @ingroup halir_api_utils
 *
 * @param v GSL float vector.
 * @param val Target value.
 * @return Index of nearest value in the vector.
 */
size_t find_nearest_index(gsl_vector_float *v, float val);

/**
 * @brief Execute the HALIR spectral calculation.
 * @ingroup halir_api_io
 *
 * @param work Validated workspace.
 * @return 0 on success.
 * @return 1 on invalid workspace state or runtime setup failure.
 */
int halir_test_calc(halir_simulation_setup *work);

/**
 * @brief Allocate a standalone spectra container with data buffers.
 * @ingroup halir_api_lifecycle
 *
 * The returned object owns `data` and `wavenum` buffers and must be released
 * with halir_spectra_free().
 *
 * @param ndatapnts Number of data samples to allocate.
 * @return Pointer to a new spectra object, or NULL on failure.
 */
halir_spectra *halir_spectra_create(size_t ndatapnts);

/**
 * @brief Free a spectra object and its owned sample buffers.
 * @ingroup halir_api_lifecycle
 *
 * This function is NULL-safe.
 *
 * @param spectra Spectra object to release.
 */
void halir_spectra_free(halir_spectra *spectra);

/**
 * @brief Allocate a result container with space for `nspectra` entries.
 * @ingroup halir_api_lifecycle
 *
 * The returned object owns the spectra array and must be released with
 * halir_result_free(). The `workspace` pointer is stored as non-owning
 * metadata and is never freed by halir_result_free().
 *
 * @param workspace Workspace pointer associated with this result.
 * @param nspectra Number of spectra entries to allocate.
 * @return Pointer to a new result object, or NULL on failure.
 */
halir_result *halir_result_create(halir_simulation_setup *workspace, size_t nspectra);

/**
 * @brief Free a result object and all owned spectra buffers.
 * @ingroup halir_api_lifecycle
 *
 * This function is NULL-safe.
 *
 * @param result Result object to release.
 */
void halir_result_free(halir_result *result);

/**
 * @brief Execute HALIR calculation and return per-composition spectra.
 * @ingroup halir_api_io
 *
 * Returns one spectrum per composition. The caller owns the returned result and
 * must release it with halir_result_free().
 *
 * The result contains deep-copied composition scalar/string metadata and
 * freshly allocated `wavenum`/`data` buffers. It does not deep-copy large
 * source prmfile arrays from workspace composition entries.
 *
 * @param work Validated workspace.
 * @return Pointer to result object on success, or NULL on failure.
 */
halir_result *halir_calculate_result(halir_simulation_setup *work);

/**
 * @brief Write one spectra object to an SPC file.
 * @ingroup halir_api_io
 *
 * The written file uses the SPC new-format little-endian layout.
 *
 * @param spectra Spectra object to serialize.
 * @param file_path Output SPC file path.
 * @return 0 on success.
 * @return 1 when arguments are invalid.
 * @return 2 when file creation or write fails.
 */
int halir_spectra_write_spc(const halir_spectra *spectra, const char *file_path);

/**
 * @brief Write a result container as a multi-subfile SPC file.
 * @ingroup halir_api_io
 *
 * Each spectra entry is written as one SPC subfile.
 *
 * @param result Result object to serialize.
 * @param file_path Output SPC file path.
 * @return 0 on success.
 * @return 1 when arguments are invalid.
 * @return 2 when file creation or write fails.
 */
int halir_result_write_spc(const halir_result *result, const char *file_path);

/**
 * @brief Read the first spectra from an SPC file.
 * @ingroup halir_api_io
 *
 * For multi-subfile SPC files, the first subfile is returned.
 *
 * @param file_path Input SPC file path.
 * @return Newly allocated spectra object on success, or NULL on failure.
 */
halir_spectra *halir_spectra_read_spc(const char *file_path);

/**
 * @brief Read an SPC file into a result container.
 * @ingroup halir_api_io
 *
 * Single-subfile SPC files produce a one-spectra result. Multi-subfile SPC
 * files produce one spectra entry per subfile.
 *
 * The returned result has a NULL workspace pointer because SPC does not store
 * enough HALIR workspace metadata to reconstruct it.
 *
 * @param file_path Input SPC file path.
 * @return Newly allocated result object on success, or NULL on failure.
 */
halir_result *halir_result_read_spc(const char *file_path);

// input str is made lower case and compared halir_Units
static inline halir_Units halir_Unit_from_str(char *str) {
  size_t length = strlen(str);
  for (int i=0; i < length; i++)
    str[i] = tolower(str[i]);

  if (strcmp(str, "atm") == 0) { return ATM; }
  else if (strcmp(str, "mbar") == 0) { return MBAR; }
  else if (strcmp(str, "bar") == 0) { return BAR; }
  else if (strcmp(str, "pa") == 0) { return PA; }
  else if (strcmp(str, "hpa") == 0) { return HPA; }
  else if (strcmp(str, "mmhg") == 0) { return MMHG; }
  else if (strcmp(str, "ppm") == 0) { return PPM; }
  else if (strcmp(str, "ppt") == 0) { return PPT; }
  else if (strcmp(str, "ppb") == 0) { return PPB; }
  else if (strcmp(str, "mm") == 0) { return MM; }
  else if (strcmp(str, "cm") == 0) { return CM; }
  else if (strcmp(str, "dm") == 0) { return DM; }
  else if (strcmp(str, "m") == 0) { return M; }
  else if (strcmp(str, "km") == 0) { return KM; }
  else if (strcmp(str, "k") == 0) { return K; }
  else if (strcmp(str, "c") == 0) { return C; }
  else if (strcmp(str, "f") == 0) { return F; }
  return NONE;
}

static inline void halir_print_simulation_setup(halir_simulation_setup *work) {
  printf("Project\n");
  printf("  pname: %s\n", work->pname);
  printf("  rootDir: %s\n", work->rootDir);
  printf("  pcomments: %s\n", work->pcomments);
  printf("Sample Environment\n");
  printf("  temp: %3.2f K\n", work->temp);
  printf("  press: %3.6f atm\n", work->press);
  //printf("%s\n", halir_Units_to_str[work->pressU]);
  printf("  pathL: %3.6f cm\n", work->pathL);
  //printf("%s\n", halir_Units_to_str[work->pathLU]);
  printf("  ROI: [%f, %f]\n", work->ROI[0], work->ROI[1]);
  printf("  res: %3.6f\n", work->res);
  printf("  fov: %3.6f\n", work->fov);
  printf("  bgfile: %s\n", work->bgfile);
  printf("  apod: %d\n", work->apod);
  printf("  ftype: %d\n", work->ftype);
  printf("Composition\n");
  for (int cc = 0; cc < work->composition_length; cc++) {
    printf("Species %d\n", cc);
    printf("   molce: %s\n", work->composition[cc]->molec);
    printf("   isotop: %s\n", work->composition[cc]->isotop);
    printf("   vmr: %2.6e\n", work->composition[cc]->vmr);
    printf("   prmfile: %s\n", work->composition[cc]->prmfile);
    printf("   nlines: %d\n", work->composition[cc]->hitran_head.ndatapnts);
  }
}

#endif
