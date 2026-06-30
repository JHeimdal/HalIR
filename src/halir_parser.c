#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <ctype.h>

#include <HalIR/halir.h>
#include <HalIR/cJSON.h>
#include <HalIR/spc.h>
#include "halir_internal.h"

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
halir_compound_load_prmfile(halir_simulation_setup *work, size_t comp_index, const char *prmfile_path)
{
  if ((work == NULL) || (prmfile_path == NULL) || (comp_index >= work->composition_length)) {
    return 1;
  }

  if (!halir_copy_string_checked(work->composition[comp_index]->prmfile, sizeof(work->composition[comp_index]->prmfile), prmfile_path, "prmfile_path")) {
    return 1;
  }

  if (!load_prmfile(work->composition[comp_index])) {
    return 2;
  }

  return 0;
}

halir_simulation_setup*
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

  halir_simulation_setup *ret_workspace = (halir_simulation_setup*)calloc(1, sizeof(halir_simulation_setup));
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
    if (!halir_copy_string_checked(ret_workspace->pname, sizeof(ret_workspace->pname), field_value->valuestring, "pname")) {
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
    if (!halir_copy_string_checked(ret_workspace->rootDir, sizeof(ret_workspace->rootDir), field_value->valuestring, "rootDir")) {
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
    if (!halir_copy_string_checked(ret_workspace->pcomments, sizeof(ret_workspace->pcomments), field_value->valuestring, "pcomments")) {
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
    if (!halir_copy_string_checked(ret_workspace->bgfile, sizeof(ret_workspace->bgfile), field_value->valuestring, "bgfile")) {
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
    if (!halir_copy_string_checked(apod, sizeof(apod), field_value->valuestring, "apod")) {
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
    if (!halir_copy_string_checked(ftype, sizeof(ftype), field_value->valuestring, "ftype")) {
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
        if (!halir_copy_string_checked(current_comp->molec, sizeof(current_comp->molec), value->valuestring, "composition.molec")) {
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
        if (!halir_copy_string_checked(current_comp->isotop, sizeof(current_comp->isotop), value->valuestring, "composition.isotop")) {
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
        if (!halir_copy_string_checked(current_comp->prmfile, sizeof(current_comp->prmfile), value->valuestring, "composition.prmfile")) {
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
    halir_simulation_setup_free(ret_workspace);
    return NULL;
  } else {
    return ret_workspace;
  }
}
