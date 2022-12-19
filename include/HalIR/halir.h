#ifndef HALIR_H_
#define HALIR_H_

#include <stdio.h>
#include <stddef.h>

#define UNIT_LEN 4
#define PATH_LEN 1024
#define COMMENT_LEN 4096

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
  ATM,
  MBAR,
  BAR,
  PA,
  HPA,
} halir_pressUnits;

typedef enum
{
  M,
  DM,
  CM,
  MM,
} halir_lengthUnits;

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
  double trans_mu;	 // Transition wavenumber (cm^-1)
  double line_I;	 // Line Intensity
  double einstein_A;	 // Einstein A-coefficient
  double air_B;		 // Air broadened with
  double self_B;	 // Self broadened width
  double low_state_en;	 // Lower state energy
  double temp_air_B;	 // Temperature dependence of air width
  double pressure_S;	 // Pressure Shift
  char   u_vib_quant[16];// Upper vibrational quanta
  char   l_vib_quant[16];// Lower vibrational quanta
  char   u_loc_quant[16];// Upper local quanta
  char   l_loc_quant[16];// Lower local quanta
  char   err_code[6];	 // Error codes
  char   ref_code[12];	 // Reference codes
  char   line_mix[2];	 // Flag for line mixing
  double u_stat_w;	 // Upper statistical weight
  double l_stat_w;	 // Lower statistical weight
  double abundance; // abundace for Isotopologue
  double molecMass; // Mass for this Isotopologue
} halir_HitranLine;

// Struct describing a compond (gas) in the system
typedef struct {
  char molec[6];
  char isotop[12];
  double vmr;
  double partialPress;
  char pressU[UNIT_LEN];
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
  double temp;
  double press;
  double pathL;
  double ROI[2];
  double res;
  double fov;
  char tempU[UNIT_LEN];
  char pressU[UNIT_LEN];
  char pathLU[UNIT_LEN];
  halir_apodization apod;
  halir_filetype ftype;
  char bgfile[PATH_LEN];
  // Composition of system
  size_t composition_length;
  halir_compound **composition;
} halir_workspace;

halir_workspace* halir_parseJSONinput(const char* const inputFile);
double CtoK(double temp) {
  return temp+273.15;
}

void halir_print_workspace(halir_workspace *work) {
  printf("Project\n");
  printf("  pname: %s\n", work->pname);
  printf("  rootDir: %s\n", work->rootDir);
  printf("  pcomments: %s\n", work->pcomments);
  printf("Sample Environment\n");
  printf("  temp: %3.2f ", work->temp);
  printf("%s\n", work->tempU);
  printf("  press: %3.6f ", work->press);
  printf("%s\n", work->pressU);
  printf("  pathL: %3.6f ", work->pathL);
  printf("%s\n", work->pathLU);
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
  }
}

#endif
