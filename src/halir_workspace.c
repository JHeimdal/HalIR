#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <HalIR/halir.h>
#include "halir_internal.h"

halir_simulation_setup*
halir_simulation_setup_create(void)
{
  return (halir_simulation_setup*)calloc(1, sizeof(halir_simulation_setup));
}

int
halir_simulation_setup_set_project(halir_simulation_setup *work, const char *pname, const char *rootDir, const char *pcomments)
{
  if ((work == NULL) || (pname == NULL) || (rootDir == NULL) || (pcomments == NULL)) {
    return 1;
  }

  if (!halir_copy_string_checked(work->pname, sizeof(work->pname), pname, "pname")) {
    return 2;
  }
  if (!halir_copy_string_checked(work->rootDir, sizeof(work->rootDir), rootDir, "rootDir")) {
    return 2;
  }
  if (!halir_copy_string_checked(work->pcomments, sizeof(work->pcomments), pcomments, "pcomments")) {
    return 2;
  }

  return 0;
}

int
halir_simulation_setup_add_project_file(halir_simulation_setup *work, const char *path)
{
  char **new_pfiles;
  size_t new_len;
  size_t path_len;

  if ((work == NULL) || (path == NULL)) {
    return 1;
  }

  new_len = work->pfiles_length + 1;
  new_pfiles = (char**)realloc(work->pfiles, new_len * sizeof(char*));
  if (new_pfiles == NULL) {
    return 2;
  }
  work->pfiles = new_pfiles;

  path_len = strlen(path);
  work->pfiles[new_len - 1] = (char*)malloc(path_len + 1);
  if (work->pfiles[new_len - 1] == NULL) {
    return 2;
  }
  memcpy(work->pfiles[new_len - 1], path, path_len + 1);
  work->pfiles_length = new_len;

  return 0;
}

int
halir_simulation_setup_set_sample_env(halir_simulation_setup *work,
                               double temp, halir_Units tempU,
                               double press, halir_Units pressU,
                               double pathL, halir_Units pathLU,
                               double roi_low, double roi_high,
                               double res, double fov,
                               halir_apodization apod,
                               halir_filetype ftype,
                               const char *bgfile)
{
  int read_err = 0;

  if ((work == NULL) || (bgfile == NULL)) {
    return 1;
  }
  if ((!isfinite(temp)) || (!isfinite(press)) || (!isfinite(pathL)) ||
      (!isfinite(roi_low)) || (!isfinite(roi_high)) ||
      (!isfinite(res)) || (!isfinite(fov))) {
    return 1;
  }
  if ((temp <= 0.0) || (press <= 0.0) || (pathL <= 0.0)) {
    return 1;
  }
  if ((roi_high <= roi_low) || (res <= 0.0) || (fov < 0.0)) {
    return 1;
  }

  work->temp = halir_Units_to_Hitran(&tempU, &temp, &read_err);
  if (read_err != 0) {
    return 2;
  }
  work->tempU = K;

  work->press = halir_Units_to_Hitran(&pressU, &press, &read_err);
  if (read_err != 0) {
    return 2;
  }
  work->pressU = ATM;

  work->pathL = halir_Units_to_Hitran(&pathLU, &pathL, &read_err);
  if (read_err != 0) {
    return 2;
  }
  work->pathLU = CM;

  work->ROI[0] = roi_low;
  work->ROI[1] = roi_high;
  work->res = res;
  work->fov = fov;
  work->apod = apod;
  work->ftype = ftype;

  if (!halir_copy_string_checked(work->bgfile, sizeof(work->bgfile), bgfile, "bgfile")) {
    return 2;
  }

  return 0;
}

int
halir_simulation_setup_add_composition(halir_simulation_setup *work, const char *molec, const char *isotop, size_t *out_index)
{
  halir_compound **new_composition;
  halir_compound *new_comp;
  size_t new_len;

  if ((work == NULL) || (molec == NULL) || (isotop == NULL) || (out_index == NULL)) {
    return 1;
  }

  new_len = work->composition_length + 1;
  new_composition = (halir_compound**)realloc(work->composition, new_len * sizeof(halir_compound*));
  if (new_composition == NULL) {
    return 2;
  }
  work->composition = new_composition;

  new_comp = (halir_compound*)calloc(1, sizeof(halir_compound));
  if (new_comp == NULL) {
    return 2;
  }

  if (!halir_copy_string_checked(new_comp->molec, sizeof(new_comp->molec), molec, "molec")) {
    free(new_comp);
    return 2;
  }
  if (!halir_copy_string_checked(new_comp->isotop, sizeof(new_comp->isotop), isotop, "isotop")) {
    free(new_comp);
    return 2;
  }

  work->composition[new_len - 1] = new_comp;
  work->composition_length = new_len;
  *out_index = new_len - 1;

  return 0;
}

int
halir_compound_set_vmr(halir_simulation_setup *work, size_t comp_index, double vmr)
{
  if ((work == NULL) || (comp_index >= work->composition_length)) {
    return 1;
  }
  if ((!isfinite(vmr)) || (vmr < 0.0) || (vmr > 1.0)) {
    return 1;
  }

  work->composition[comp_index]->vmr = vmr;
  return 0;
}

int
halir_compound_set_concentration(halir_simulation_setup *work, size_t comp_index, double conc, halir_Units concU)
{
  int read_err = 0;
  double normalized_conc;

  if ((work == NULL) || (comp_index >= work->composition_length)) {
    return 1;
  }
  if ((!isfinite(conc)) || (conc < 0.0)) {
    return 1;
  }

  normalized_conc = halir_Units_to_Hitran(&concU, &conc, &read_err);
  if (read_err != 0) {
    return 2;
  }

  work->composition[comp_index]->conc = normalized_conc;
  work->composition[comp_index]->concU = PPM;

  return 0;
}

void halir_simulation_setup_free(halir_simulation_setup *work)
{
  if (work == NULL) {
    return;
  }

  if (work->pfiles != NULL) {
    for (size_t i = 0; i < work->pfiles_length; i++) {
      free(work->pfiles[i]);
    }
    free(work->pfiles);
    work->pfiles = NULL;
  }

  if (work->composition != NULL) {
    for (size_t i = 0; i < work->composition_length; i++) {
      if (work->composition[i] != NULL) {
        free(work->composition[i]->hitran_head.molecs);
        free(work->composition[i]->hitran_prms);
        free(work->composition[i]);
      }
    }
    free(work->composition);
    work->composition = NULL;
  }

  free(work);
}

int
halir_simulation_setup_validate(halir_simulation_setup *work)
{
  if (work == NULL) {
    return 1;
  }
  if ((strlen(work->pname) == 0) || (strlen(work->rootDir) == 0)) {
    fprintf(stderr, "Validation: project pname or rootDir not set\n");
    return 1;
  }
  if ((!isfinite(work->temp)) || (!isfinite(work->press)) || (!isfinite(work->pathL))) {
    fprintf(stderr, "Validation: sample environment has non-finite values\n");
    return 1;
  }
  if ((work->temp <= 0.0) || (work->press <= 0.0) || (work->pathL <= 0.0)) {
    fprintf(stderr, "Validation: sample environment has non-positive values\n");
    return 1;
  }
  if ((!isfinite(work->ROI[0])) || (!isfinite(work->ROI[1])) ||
      (work->ROI[1] <= work->ROI[0])) {
    fprintf(stderr, "Validation: ROI invalid (must have high > low, both finite)\n");
    return 1;
  }
  if ((!isfinite(work->res)) || (work->res <= 0.0)) {
    fprintf(stderr, "Validation: resolution must be positive and finite\n");
    return 1;
  }
  if ((!isfinite(work->fov)) || (work->fov < 0.0)) {
    fprintf(stderr, "Validation: FOV must be non-negative and finite\n");
    return 1;
  }
  if (work->composition_length == 0) {
    fprintf(stderr, "Validation: no compositions defined\n");
    return 1;
  }
  for (size_t i = 0; i < work->composition_length; i++) {
    if (work->composition[i] == NULL) {
      fprintf(stderr, "Validation: composition %zu is NULL\n", i);
      return 1;
    }
    if ((strlen(work->composition[i]->molec) == 0) || (strlen(work->composition[i]->isotop) == 0)) {
      fprintf(stderr, "Validation: composition %zu missing molec or isotop\n", i);
      return 1;
    }
    if ((!isfinite(work->composition[i]->vmr)) || (work->composition[i]->vmr < 0.0) || (work->composition[i]->vmr > 1.0)) {
      fprintf(stderr, "Validation: composition %zu VMR out of range [0,1] or non-finite\n", i);
      return 1;
    }
  }
  return 0;
}
