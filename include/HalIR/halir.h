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

typedef enum
{
  HALIR_BOXCAR,
  HALIR_TRIANGLE,
  HALIR_HAPPGENZEL,
} halir_apodization;

typedef enum
{
  HALIR_TRANSMISSION,
  HALIR_ABSORBANCE,
} halir_filetype;

typedef enum
{
  NONE,
  ATM,
  MBAR,
  BAR,
  PA,
  HPA,
  MMHG,
  PPM,
  PPT,
  PPB,
  MM,
  CM,
  DM,
  M,
  KM,
  K,
  C,
  F,
} halir_Units;

// Convert halir_Units to atm, cm, Kelvin
double halir_Units_factors[] = {
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

char *halir_Units_to_str[] = {
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

double CtoK(const double *temp) {
  return *temp + 273.15;
}
double FtoK(const double *temp) {
  return (*temp + 459.67) * 5./9.;
}
double halir_Units_to_Hitran(const halir_Units *unit, const double *value, int *read_err) {
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

typedef struct {
  char molec[6];
  char isoname[12];
  int nisotp;
  float roi_low;
  float roi_high;
  int ndatapnts;
  int *molecs;
} halir_HitranHead;

typedef struct {
  int    molec_num;	 // Molecule Number
  int    isotp_num;	 // Isotopologue Number
  halir_num trans_mu;	 // Transition wavenumber (cm^-1)
  halir_num line_I;	 // Line Intensity
  halir_num einstein_A;	 // Einstein A-coefficient
  halir_num air_B;		 // Air broadened with
  halir_num self_B;	 // Self broadened width
  halir_num low_state_en;	 // Lower state energy
  halir_num temp_air_B;	 // Temperature dependence of air width
  halir_num pressure_S;	 // Pressure Shift
  char   u_vib_quant[16];// Upper vibrational quanta
  char   l_vib_quant[16];// Lower vibrational quanta
  char   u_loc_quant[16];// Upper local quanta
  char   l_loc_quant[16];// Lower local quanta
  char   err_code[6];	 // Error codes
  char   ref_code[12];	 // Reference codes
  char   line_mix[2];	 // Flag for line mixing
  halir_num u_stat_w;	 // Upper statistical weight
  halir_num l_stat_w;	 // Lower statistical weight
  halir_num abundance; // abundace for Isotopologue
  halir_num molecMass; // Mass for this Isotopologue
} halir_HitranLine;

// Struct describing a compond (gas) in the system
typedef struct {
  char molec[6];
  char isotop[12];
  halir_num vmr;
  halir_num conc;
  halir_Units concU;
  char prmfile[PATH_LEN];
  halir_HitranHead hitran_head;
  halir_HitranLine *hitran_prms;
} halir_compound;

// Struct defining the system
typedef struct
{
  // project descriptors
  char pname[PATH_LEN];
  char rootDir[PATH_LEN];
  char pcomments[COMMENT_LEN];
  char **pfiles;
  size_t pfiles_length;

  // Environment of system
  halir_num temp;
  halir_num press;
  halir_num pathL;
  halir_num ROI[2];
  halir_num res;
  halir_num fov;
  halir_Units tempU;
  halir_Units pressU;
  halir_Units pathLU;
  halir_apodization apod;
  halir_filetype ftype;
  char bgfile[PATH_LEN];
  // Composition of system
  size_t composition_length;
  halir_compound **composition;
} halir_workspace;

// Result of simulated spectra
typedef struct{
  size_t ndatapnts;
  halir_compound *composition;
  halir_num *data;
} hair_spectra;

halir_workspace* halir_parseJSONinput(const char* const inputFile);
size_t find_nearest_index(gsl_vector_float *v, float val);
int halir_test_calc(halir_workspace *work);

// input str is made lower case and compared halir_Units
halir_Units halir_Unit_from_str(char *str) {
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

void halir_print_workspace(halir_workspace *work) {
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
