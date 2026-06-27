#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_vector_float.h>
#include <gsl/gsl_math.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#include <HalIR/halir.h>
#include <HalIR/tips.h>
#include <HalIR/cJSON.h>
#include <cerf.h>

#ifdef HALIR_HAVE_OPENMP
#include <omp.h>
#endif

static int copy_string_checked(char *dst, size_t dst_size, const char *src, const char *field_name)
{
  size_t src_len;

  if ((dst == NULL) || (src == NULL) || (dst_size == 0)) {
    fprintf(stderr, "Invalid string copy args for field: %s\n", field_name);
    return 0;
  }

  src_len = strlen(src);
  if (src_len >= dst_size) {
    fprintf(stderr, "Field '%s' exceeds max length (%zu)\n", field_name, dst_size - 1);
    return 0;
  }

  memcpy(dst, src, src_len + 1);
  return 1;
}

halir_workspace*
halir_workspace_create(void)
{
  return (halir_workspace*)calloc(1, sizeof(halir_workspace));
}

int
halir_workspace_set_project(halir_workspace *work, const char *pname, const char *rootDir, const char *pcomments)
{
  if ((work == NULL) || (pname == NULL) || (rootDir == NULL) || (pcomments == NULL)) {
    return 1;
  }

  if (!copy_string_checked(work->pname, sizeof(work->pname), pname, "pname")) {
    return 2;
  }
  if (!copy_string_checked(work->rootDir, sizeof(work->rootDir), rootDir, "rootDir")) {
    return 2;
  }
  if (!copy_string_checked(work->pcomments, sizeof(work->pcomments), pcomments, "pcomments")) {
    return 2;
  }

  return 0;
}

int
halir_workspace_add_project_file(halir_workspace *work, const char *path)
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
halir_workspace_set_sample_env(halir_workspace *work,
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

  if (!copy_string_checked(work->bgfile, sizeof(work->bgfile), bgfile, "bgfile")) {
    return 2;
  }

  return 0;
}

int
halir_workspace_add_composition(halir_workspace *work, const char *molec, const char *isotop, size_t *out_index)
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

  if (!copy_string_checked(new_comp->molec, sizeof(new_comp->molec), molec, "molec")) {
    free(new_comp);
    return 2;
  }
  if (!copy_string_checked(new_comp->isotop, sizeof(new_comp->isotop), isotop, "isotop")) {
    free(new_comp);
    return 2;
  }

  work->composition[new_len - 1] = new_comp;
  work->composition_length = new_len;
  *out_index = new_len - 1;

  return 0;
}

int
halir_compound_set_vmr(halir_workspace *work, size_t comp_index, double vmr)
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
halir_compound_set_concentration(halir_workspace *work, size_t comp_index, double conc, halir_Units concU)
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

static int
load_prmfile(halir_compound *current_comp)
{
  FILE *fp;
  size_t readSize;

  if (access(current_comp->prmfile, F_OK|R_OK) != 0) {
    fprintf(stderr, "Cant access the prmfile: %s\n", current_comp->prmfile);
    return 0;
  }

  fp = fopen(current_comp->prmfile, "rb");
  if (fp == NULL) {
    fprintf(stderr, "Could not open prmfile: %s\n", current_comp->prmfile);
    return 0;
  }

  readSize = fread(&current_comp->hitran_head, sizeof(halir_HitranHead), 1, fp);
  if (readSize != 1) {
    fprintf(stderr, "Error while reading prmfile header: %s\n", current_comp->prmfile);
    fclose(fp);
    return 0;
  }

  if (current_comp->hitran_head.nisotp < 0 || current_comp->hitran_head.ndatapnts < 0) {
    fprintf(stderr, "Invalid prmfile metadata in: %s\n", current_comp->prmfile);
    fclose(fp);
    return 0;
  }

  current_comp->hitran_head.molecs = NULL;
  current_comp->hitran_prms = NULL;

  if (current_comp->hitran_head.nisotp > 0) {
    current_comp->hitran_head.molecs = (int*)malloc(sizeof(int) * (size_t)current_comp->hitran_head.nisotp);
    if (current_comp->hitran_head.molecs == NULL) {
      fprintf(stderr, "Memory allocation failed while reading prmfile: %s\n", current_comp->prmfile);
      fclose(fp);
      return 0;
    }

    readSize = fread(current_comp->hitran_head.molecs, sizeof(int), (size_t)current_comp->hitran_head.nisotp, fp);
    if (readSize != (size_t)current_comp->hitran_head.nisotp) {
      fprintf(stderr, "Error while reading prmfile isotopes: %s\n", current_comp->prmfile);
      fclose(fp);
      return 0;
    }
  }

  if (current_comp->hitran_head.ndatapnts > 0) {
    current_comp->hitran_prms = (halir_HitranLine*)malloc(sizeof(halir_HitranLine) * (size_t)current_comp->hitran_head.ndatapnts);
    if (current_comp->hitran_prms == NULL) {
      fprintf(stderr, "Memory allocation failed while reading prmfile lines: %s\n", current_comp->prmfile);
      fclose(fp);
      return 0;
    }

    readSize = fread(current_comp->hitran_prms, sizeof(halir_HitranLine), (size_t)current_comp->hitran_head.ndatapnts, fp);
    if (readSize != (size_t)current_comp->hitran_head.ndatapnts) {
      fprintf(stderr, "Error while reading prmfile lines: %s\n", current_comp->prmfile);
      fclose(fp);
      return 0;
    }
  }

  fclose(fp);
  return 1;
}

int
halir_compound_load_prmfile(halir_workspace *work, size_t comp_index, const char *prmfile_path)
{
  if ((work == NULL) || (prmfile_path == NULL) || (comp_index >= work->composition_length)) {
    return 1;
  }

  if (!copy_string_checked(work->composition[comp_index]->prmfile, sizeof(work->composition[comp_index]->prmfile), prmfile_path, "prmfile_path")) {
    return 1;
  }

  if (!load_prmfile(work->composition[comp_index])) {
    return 2;
  }

  return 0;
}

void halir_workspace_free(halir_workspace *work)
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

halir_spectra *
halir_spectra_create(size_t ndatapnts)
{
  halir_spectra *spectra;

  if (ndatapnts == 0) {
    return NULL;
  }

  spectra = (halir_spectra*)calloc(1, sizeof(halir_spectra));
  if (spectra == NULL) {
    return NULL;
  }

  spectra->wavenum = (halir_num*)calloc(ndatapnts, sizeof(halir_num));
  spectra->data = (halir_num*)calloc(ndatapnts, sizeof(halir_num));
  if ((spectra->wavenum == NULL) || (spectra->data == NULL)) {
    free(spectra->wavenum);
    free(spectra->data);
    free(spectra);
    return NULL;
  }

  spectra->ndatapnts = ndatapnts;
  return spectra;
}

void
halir_spectra_free(halir_spectra *spectra)
{
  if (spectra == NULL) {
    return;
  }

  free(spectra->wavenum);
  free(spectra->data);
  free(spectra->composition.hitran_head.molecs);
  free(spectra->composition.hitran_prms);
  free(spectra);
}

halir_result *
halir_result_create(halir_workspace *workspace, size_t nspectra)
{
  halir_result *result;

  if ((workspace == NULL) || (nspectra == 0)) {
    return NULL;
  }

  result = (halir_result*)calloc(1, sizeof(halir_result));
  if (result == NULL) {
    return NULL;
  }

  result->spectra = (halir_spectra*)calloc(nspectra, sizeof(halir_spectra));
  if (result->spectra == NULL) {
    free(result);
    return NULL;
  }

  result->workspace = workspace;
  result->nspectra = nspectra;
  return result;
}

void
halir_result_free(halir_result *result)
{
  size_t i;

  if (result == NULL) {
    return;
  }

  if (result->spectra != NULL) {
    for (i = 0; i < result->nspectra; i++) {
      free(result->spectra[i].wavenum);
      free(result->spectra[i].data);
      free(result->spectra[i].composition.hitran_head.molecs);
      free(result->spectra[i].composition.hitran_prms);
    }
    free(result->spectra);
  }

  free(result);
}

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

int
halir_workspace_validate(halir_workspace *work)
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

halir_workspace*
halir_parseJSONinput(const char* const inputFile)
{
  // Read the text file
  FILE *fp = NULL;
  size_t counter;
  size_t file_length;
  int arrayLength;
  int read_err = 0;
  const char *infile;
  char *infile_buffer = NULL;
  int infile_is_allocated = 0;

  const cJSON *input_field = NULL;
  const cJSON *project_field = NULL;
  const cJSON *sampleEnv_field = NULL;
  const cJSON *composition_field = NULL;
  const cJSON *field_value = NULL;

  if (access(inputFile, F_OK|R_OK) == 0) {
    fp = fopen(inputFile, "r");
    if (fp == NULL) {
      fprintf(stderr, "Could not open input file: %s\n", inputFile);
      read_err = 1;
      goto end;
    }
    if (fseek(fp, 0, SEEK_END) != 0) {
      fprintf(stderr, "Could not seek input file: %s\n", inputFile);
      read_err = 1;
      goto end;
    }
    file_length = ftell(fp);
    if ((long)file_length < 0) {
      fprintf(stderr, "Could not determine input file length: %s\n", inputFile);
      read_err = 1;
      goto end;
    }
    rewind(fp);
    infile_buffer = (char*)malloc((file_length+1)*sizeof(char));
    if (infile_buffer == NULL) {
      fprintf(stderr, "Memory allocation failed for input file buffer\n");
      read_err = 1;
      goto end;
    }
    if (fread(infile_buffer, sizeof(char), file_length, fp) != file_length) {
      fprintf(stderr, "Could not read full input file: %s\n", inputFile);
      read_err = 1;
      goto end;
    }
    fclose(fp);
    fp = NULL;
    infile_buffer[file_length] = '\0';
    infile = infile_buffer;
    infile_is_allocated = 1;
  } else {
    infile = inputFile;
  }

  halir_workspace *ret_workspace = (halir_workspace*)calloc(1, sizeof(halir_workspace));
  if (ret_workspace == NULL) {
    fprintf(stderr, "Memory allocation failed for workspace\n");
    read_err = 1;
    goto end;
  }

  cJSON *root = cJSON_Parse(infile);
  if (root == NULL) {
    const char *err_ptr = cJSON_GetErrorPtr();
    if (err_ptr != NULL) {
      fprintf(stderr, "Error before: %s\n", err_ptr);
      read_err = 1;
    }
    read_err = 2;
    goto end;
  }
  input_field = cJSON_GetObjectItem(root, "input");
  if (!cJSON_IsObject(input_field)) {
    fprintf(stderr, "Missing input object or wrong type\n");
    read_err = 3;
    goto end;
  }
  project_field = cJSON_GetObjectItem(input_field, "project");
  if (!cJSON_IsObject(project_field)) {
    fprintf(stderr, "Missing project object or wrong type\n");
    read_err = 3;
    goto end;
  }
  field_value = cJSON_GetObjectItem(project_field, "pname");
  if (cJSON_IsString(field_value)) {
    if (!copy_string_checked(ret_workspace->pname, sizeof(ret_workspace->pname), field_value->valuestring, "pname")) {
      read_err = 3;
      goto end;
    }
    /*printf("pname: %s\n", ret_workspace->pname);*/
  } else {
    fprintf(stderr, "Missing temp keyword or wrong type\n");
    read_err = 3;
    goto end;
  }
  field_value = cJSON_GetObjectItem(project_field, "rootDir");
  if (cJSON_IsString(field_value)) {
    if (!copy_string_checked(ret_workspace->rootDir, sizeof(ret_workspace->rootDir), field_value->valuestring, "rootDir")) {
      read_err = 3;
      goto end;
    }
    /*printf("rootDir: %s\n", ret_workspace->rootDir);*/
  } else {
    fprintf(stderr, "Missing rootDir keyword or wrong type\n");
    read_err = 3;
    goto end;
  }
  field_value = cJSON_GetObjectItem(project_field, "pcomments");
  if (cJSON_IsString(field_value)) {
    if (!copy_string_checked(ret_workspace->pcomments, sizeof(ret_workspace->pcomments), field_value->valuestring, "pcomments")) {
      read_err = 3;
      goto end;
    }
    /*printf("pcomments: %s\n", ret_workspace->pcomments);*/
  } else {
    fprintf(stderr, "Missing pcomments keyword or wrong type\n");
    read_err = 3;
    goto end;
  }
  field_value = cJSON_GetObjectItem(project_field, "pfiles");
  counter = 0;
  cJSON *pfile = NULL;
  if (cJSON_IsArray(field_value)) {
    arrayLength = cJSON_GetArraySize(field_value);
    ret_workspace->pfiles = (char**)calloc((size_t)arrayLength, sizeof(char*));
    if (ret_workspace->pfiles == NULL) {
      fprintf(stderr, "Could not allocate pfiles array\n");
      read_err = 3;
      goto end;
    }
    ret_workspace->pfiles_length = arrayLength;
    cJSON_ArrayForEach(pfile, field_value) {
      if (cJSON_IsString(pfile)) {
        size_t pfile_len = strlen(pfile->valuestring);
        ret_workspace->pfiles[counter] = (char*)malloc(pfile_len + 1);
        if (ret_workspace->pfiles[counter] == NULL) {
          fprintf(stderr, "Could not allocate pfiles entry\n");
          read_err = 3;
          goto end;
        }
        memcpy(ret_workspace->pfiles[counter], pfile->valuestring, pfile_len + 1);
        /*printf("pfiles: %s\n", ret_workspace->pfiles[counter]);*/
      } else if (cJSON_IsNull(pfile)) {
        ret_workspace->pfiles[counter] = (char*)malloc(1);
        if (ret_workspace->pfiles[counter] == NULL) {
          fprintf(stderr, "Could not allocate pfiles null entry\n");
          read_err = 3;
          goto end;
        }
        ret_workspace->pfiles[counter][0] = '\0';
      } else {
        fprintf(stderr, "Invalid pfiles entry type\n");
        read_err = 3;
        goto end;
      }
      counter++;
    }
  } else {
    fprintf(stderr, "Missing pfiles keyword or wrong type\n");
    read_err = 3;
    goto end;
  }

  // Read the sampleEnv field
  sampleEnv_field = cJSON_GetObjectItem(input_field, "sampleEnv");
  if (!cJSON_IsObject(sampleEnv_field)) {
    fprintf(stderr, "Missing sampleEnv object or wrong type\n");
    read_err = 3;
    goto end;
  }
  field_value = cJSON_GetObjectItem(sampleEnv_field, "tempU");
  if (cJSON_IsString(field_value)) {
    ret_workspace->tempU = halir_Unit_from_str(field_value->valuestring);
  } else {
    fprintf(stderr, "Missing tempU keyword or wrong type\n");
    read_err = 3;
    goto end;
  }
  field_value = cJSON_GetObjectItem(sampleEnv_field, "temp");
  if (cJSON_IsNumber(field_value)) {
    ret_workspace->temp = halir_Units_to_Hitran(&ret_workspace->tempU, &field_value->valuedouble, &read_err);
    // Converted to from tempU -> kelvin (halir_Units) update workspace tempU
    ret_workspace->tempU = K;
    if (read_err != 0) {
      fprintf(stderr, "Temperature unit not supported\n");
      goto end;
    }
  } else {
    fprintf(stderr, "Missing temp keyword or wrong type\n");
    read_err = 3;
    goto end;
  }
  field_value = cJSON_GetObjectItem(sampleEnv_field, "pressU");
  if (cJSON_IsString(field_value)) {
    ret_workspace->pressU = halir_Unit_from_str(field_value->valuestring);
    /*printf("%s\n", ret_workspace->pressU);*/
  } else {
    fprintf(stderr, "Missing pressU keyword or wrong type\n");
    read_err = 3;
    goto end;
  }
  field_value = cJSON_GetObjectItem(sampleEnv_field, "press");
  if (cJSON_IsNumber(field_value)) {
    /*printf("%s: fac: %e\n", halir_Units_to_str[ret_workspace->pressU], halir_Units_factors[ret_workspace->pressU]);*/
    ret_workspace->press = halir_Units_to_Hitran(&ret_workspace->pressU, &field_value->valuedouble, &read_err);
    // Converted to from pressU -> atm (halir_Units) update workspace pressU
    ret_workspace->pressU = ATM;
    if (read_err != 0) {
      fprintf(stderr, "Pressure unit not supported\n");
      goto end;
    }
  } else {
    fprintf(stderr, "Missing press keyword or wrong type\n");
    read_err = 3;
    goto end;
  }
  field_value = cJSON_GetObjectItem(sampleEnv_field, "pathLU");
  if (cJSON_IsString(field_value)) {
    ret_workspace->pathLU = halir_Unit_from_str(field_value->valuestring);
    /*printf("%s\n", ret_workspace->pathLU);*/
  } else {
    fprintf(stderr, "Missing pathLU keyword or wrong type\n");
    read_err = 3;
    goto end;
  }
  field_value = cJSON_GetObjectItem(sampleEnv_field, "pathL");
  if (cJSON_IsNumber(field_value)) {
    /*printf("%s: fac: %e\n", halir_Units_to_str[ret_workspace->pathLU], halir_Units_factors[ret_workspace->pathLU]);*/
    ret_workspace->pathL = halir_Units_to_Hitran(&ret_workspace->pathLU, &field_value->valuedouble, &read_err);
    if (read_err != 0) {
      fprintf(stderr, "Path length unit not supported\n");
      goto end;
    }
  } else {
    fprintf(stderr, "Missing pathL keyword or wrong type\n");
    read_err = 3;
    goto end;
  }
  field_value = cJSON_GetObjectItem(sampleEnv_field, "bgfile");
  if (cJSON_IsString(field_value)) {
    if (!copy_string_checked(ret_workspace->bgfile, sizeof(ret_workspace->bgfile), field_value->valuestring, "bgfile")) {
      read_err = 3;
      goto end;
    }
    /*printf("%s\n", ret_workspace->pathLU);*/
  } else {
    fprintf(stderr, "Missing bgfile keyword or wrong type\n");
    read_err = 3;
    goto end;
  }
  field_value = cJSON_GetObjectItem(sampleEnv_field, "ROI");
  counter = 0;
  cJSON *value = NULL;
  if (cJSON_IsArray(field_value)) {
    arrayLength = cJSON_GetArraySize(field_value);
    if (arrayLength != 2) {
      fprintf(stderr, "ROI must contain exactly 2 values\n");
      read_err = 3;
      goto end;
    }
    cJSON_ArrayForEach(value, field_value) {
      if (!cJSON_IsNumber(value)) {
        fprintf(stderr, "ROI entries must be numeric\n");
        read_err = 3;
        goto end;
      }
      ret_workspace->ROI[counter] = value->valuedouble;
      counter++;
    }
  } else {
    fprintf(stderr, "Missing ROI keyword or wrong type\n");
    read_err = 3;
    goto end;
  }
  field_value = cJSON_GetObjectItem(sampleEnv_field, "res");
  if (cJSON_IsNumber(field_value)) {
    ret_workspace->res = field_value->valuedouble;
    /*printf("pathL: %3.1f ", ret_workspace->pathL);*/
  } else {
    fprintf(stderr, "Missing res keyword or wrong type\n");
    read_err = 3;
    goto end;
  }
  field_value = cJSON_GetObjectItem(sampleEnv_field, "fov");
  if (cJSON_IsNumber(field_value)) {
    ret_workspace->fov = field_value->valuedouble;
    /*printf("pathL: %3.1f ", ret_workspace->pathL);*/
  } else {
    fprintf(stderr, "Missing fov keyword or wrong type\n");
    read_err = 3;
    goto end;
  }
  field_value = cJSON_GetObjectItem(sampleEnv_field, "apod");
  if (cJSON_IsString(field_value)) {
    char apod[96];
    if (!copy_string_checked(apod, sizeof(apod), field_value->valuestring, "apod")) {
      read_err = 3;
      goto end;
    }
    for (int i=0; i < strlen(apod); i++)
      apod[i] = tolower(apod[i]);
    if (strcmp(apod, "boxcar") == 0)
      ret_workspace->apod = HALIR_BOXCAR;
    else if (strcmp(apod, "triangle") == 0)
      ret_workspace->apod = HALIR_TRIANGLE;
    else if (strcmp(apod, "happgenzel") == 0)
      ret_workspace->apod = HALIR_HAPPGENZEL;
    else {
      fprintf(stderr, "No \"%s\" apodisation function implemented\n", apod);
      read_err = 3;
      goto end;
    }
    /*printf("pathL: %3.1f ", ret_workspace->pathL);*/
  } else {
    fprintf(stderr, "Missing apod keyword or wrong type\n");
    read_err = 3;
    goto end;
  }
  field_value = cJSON_GetObjectItem(sampleEnv_field, "ftype");
  if (cJSON_IsString(field_value)) {
    char ftype[96];
    if (!copy_string_checked(ftype, sizeof(ftype), field_value->valuestring, "ftype")) {
      read_err = 3;
      goto end;
    }
    for (int i=0; i < strlen(ftype); i++)
      ftype[i] = tolower(ftype[i]);
    if (strcmp(ftype, "transmission") == 0)
      ret_workspace->ftype = HALIR_TRANSMISSION;
    else if (strcmp(ftype, "absorbance") == 0)
      ret_workspace->ftype = HALIR_ABSORBANCE;
    else {
      fprintf(stderr, "No \"%s\" file type implemented\n", ftype);
      read_err = 3;
      goto end;
    }
    /*printf("pathL: %3.1f ", ret_workspace->pathL);*/
  } else {
    fprintf(stderr, "Missing ftype keyword or wrong type\n");
    read_err = 3;
    goto end;
  }
  // Read the compositions
  composition_field = cJSON_GetObjectItem(input_field, "composition");
  if (cJSON_IsArray(composition_field)) {
    arrayLength = cJSON_GetArraySize(composition_field);
    ret_workspace->composition_length = arrayLength;
    ret_workspace->composition = (halir_compound**)calloc((size_t)arrayLength, sizeof(halir_compound*));
    if (ret_workspace->composition == NULL) {
      fprintf(stderr, "Could not allocate composition array\n");
      read_err = 4;
      goto end;
    }
    halir_compound *current_comp;
    cJSON *comp;
    counter = 0;
    cJSON_ArrayForEach(comp, composition_field) {
      ret_workspace->composition[counter] = (halir_compound*)malloc(sizeof(halir_compound));
      if (ret_workspace->composition[counter] == NULL) {
        fprintf(stderr, "Could not allocate composition entry\n");
        read_err = 4;
        goto end;
      }
      memset(ret_workspace->composition[counter], 0, sizeof(halir_compound));
      current_comp = ret_workspace->composition[counter];
      current_comp->hitran_prms = NULL;
      value = cJSON_GetObjectItem(comp, "molec");
      if (cJSON_IsString(value)) {
        if (!copy_string_checked(current_comp->molec, sizeof(current_comp->molec), value->valuestring, "composition.molec")) {
          read_err = 4;
          goto end;
        }
        /*printf("molec: %s\n", current_comp->molec);*/
      } else {
        fprintf(stderr, "Missing molec keyword or wrong type\n");
        read_err = 4;
        goto end;
      }
      value = cJSON_GetObjectItem(comp, "isotop");
      if (cJSON_IsString(value)) {
        if (!copy_string_checked(current_comp->isotop, sizeof(current_comp->isotop), value->valuestring, "composition.isotop")) {
          read_err = 4;
          goto end;
        }
        /*printf("isotop: %s\n", current_comp->isotop);*/
      } else {
        fprintf(stderr, "Missing isotop keyword or wrong type\n");
        read_err = 4;
        goto end;
      }
      // Here are are three different cases of keywords
      // 1. There can be only concU and conc
      // 2, There can be only vmr
      // 3, Both 1 and 2
      cJSON *concU = cJSON_GetObjectItem(comp, "concU");
      cJSON *conc = cJSON_GetObjectItem(comp, "conc");
      cJSON *vmr = cJSON_GetObjectItem(comp, "vmr");
      if ((concU != NULL) && (conc != NULL) && (vmr == NULL)) {
        if (cJSON_IsString(concU)) {
          current_comp->concU = halir_Unit_from_str(concU->valuestring);
        } else {
          fprintf(stderr, "Missing pressU keyword or wrong type\n");
          read_err = 4;
          goto end;
        }
        if (cJSON_IsNumber(conc)) {
          switch (current_comp->concU) {
            case PPM:
            case PPT:
            case PPB:
              current_comp->conc = halir_Units_to_Hitran(&current_comp->concU, &conc->valuedouble, &read_err);
              current_comp->vmr = current_comp->conc;
              break;
            default:
              current_comp->conc = halir_Units_to_Hitran(&current_comp->concU, &conc->valuedouble, &read_err);
              current_comp->vmr = current_comp->conc/ret_workspace->press;
          }
          if (read_err != 0) {
            fprintf(stderr, "Conc unit not supported\n");
            goto end;
          }
        } else {
          fprintf(stderr, "Missing conc keyword or wrong type\n");
          read_err = 4;
          goto end;
        }
      } else if ((concU == NULL) && (conc == NULL) && (vmr != NULL)) {
        if (cJSON_IsNumber(vmr)) {
          current_comp->vmr = vmr->valuedouble;
          /*printf("vmr: %e\n", current_comp->vmr);*/
        } else {
          fprintf(stderr, "Missing vmr keyword or wrong type\n");
          read_err = 4;
          goto end;
        }
      } else if ((concU != NULL) && (conc != NULL) && (vmr != NULL)) {
        if (cJSON_IsString(concU)) {
          current_comp->concU = halir_Unit_from_str(concU->valuestring);
          /*printf("pressU: %s\n", current_comp->pressU);*/
        } else {
          fprintf(stderr, "Missing pressU keyword or wrong type\n");
          read_err = 4;
          goto end;
        }
        if (cJSON_IsNumber(conc)) {
          current_comp->conc = halir_Units_factors[current_comp->concU] * conc->valuedouble;
          /*printf("pressU: %s\n", current_comp->pressU);*/
        } else {
          fprintf(stderr, "Missing pressU keyword or wrong type\n");
          read_err = 4;
          goto end;
        }
        if (cJSON_IsNumber(vmr)) {
          current_comp->vmr = vmr->valuedouble;
          /*printf("vmr: %e\n", current_comp->vmr);*/
        } else {
          fprintf(stderr, "Missing vmr keyword or wrong type\n");
          read_err = 4;
          goto end;
        }
      } else {
        fprintf(stderr, "Missing concentration keyword either vmr or a conc, concU pair needed\n");
        read_err = 4;
        goto end;
      }
      value = cJSON_GetObjectItem(comp, "prmfile");
      if (cJSON_IsString(value)) {
        if (!copy_string_checked(current_comp->prmfile, sizeof(current_comp->prmfile), value->valuestring, "composition.prmfile")) {
          read_err = 4;
          goto end;
        }
        /*printf("prmfile: %s\n", current_comp->prmfile);*/
      } else {
        fprintf(stderr, "Missing prmfile or wrong type\n");
        read_err = 4;
        goto end;
      }
      // Is there a file!
      if (strcmp(current_comp->prmfile, "") != 0) {
        if (!load_prmfile(current_comp)) {
          read_err = 5;
          goto end;
        }

        /*printf("From hitran file\n");*/
        /*printf("molec: %s\n", current_comp->hitran_head.molec);*/
        /*printf("ROI: %f %f\n", current_comp->hitran_head.roi_low,current_comp->hitran_head.roi_high);*/
        /*printf("ndatapnts: %d\n", current_comp->hitran_head.ndatapnts);*/
      }
      /*for (int i=0; i < current_comp->hitran_head.ndatapnts; i++) {*/
               /*//1 2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21*/
        /*printf("%d%d %f %e %f %f %f %f %f %f %s %s %s %s %s %s %s %f %f %f %f\n",*/
            /*current_comp->hitran_prms[i].molec_num,*/
            /*current_comp->hitran_prms[i].isotp_num,*/
            /*current_comp->hitran_prms[i].trans_mu,*/
            /*current_comp->hitran_prms[i].line_I,*/
            /*current_comp->hitran_prms[i].einstein_A,*/
            /*current_comp->hitran_prms[i].air_B,*/
            /*current_comp->hitran_prms[i].self_B,*/
            /*current_comp->hitran_prms[i].low_state_en,*/
            /*current_comp->hitran_prms[i].temp_air_B,*/
            /*current_comp->hitran_prms[i].pressure_S,*/
            /*current_comp->hitran_prms[i].u_vib_quant,*/
            /*current_comp->hitran_prms[i].l_vib_quant,*/
            /*current_comp->hitran_prms[i].u_loc_quant,*/
            /*current_comp->hitran_prms[i].l_loc_quant,*/
            /*current_comp->hitran_prms[i].err_code,*/
            /*current_comp->hitran_prms[i].ref_code,*/
            /*current_comp->hitran_prms[i].line_mix,*/
            /*current_comp->hitran_prms[i].u_stat_w,*/
            /*current_comp->hitran_prms[i].l_stat_w,*/
            /*current_comp->hitran_prms[i].abundance,*/
            /*current_comp->hitran_prms[i].molecMass);*/
      /*}*/
      counter++;
    }
  } else {
    fprintf(stderr, "Missing composition keyword or wrong type\n");
    read_err = 3;
    goto end;
  }

end:
  if (fp != NULL) {
    fclose(fp);
  }
  cJSON_Delete(root);
  if (infile_is_allocated) {
    free(infile_buffer);
  }
  if (read_err) {
    fprintf(stderr, "Something went wrong while reading input\n");
    halir_workspace_free(ret_workspace);
    return NULL;
  } else {
    return ret_workspace;
  }
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

halir_result *halir_calculate_result(halir_workspace *work)
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
  const double c_si = 299792458.; // m/s — SI speed of light; used in MM (molecular speed term)
  const double c_erg = 29979245800.; // cm/s — CGS speed of light; used via hc_erg
  const double atmmass_si=1.6605e-27; // kg/amu — atomic mass unit; used in MM
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

int halir_test_calc(halir_workspace *work)
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
