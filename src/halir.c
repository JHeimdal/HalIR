#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <ctype.h>

#include "HalIR/halir.h"
#include "HalIR/cJSON.h"

halir_workspace*
halir_parseJSONinput(const char* const inputFile)
{
  // Read the text file
  FILE *fp;
  size_t counter;
  size_t readSize;
  int arrayLength;
  int read_err = 0;
  const cJSON *input_field = NULL;
  const cJSON *project_field = NULL;
  const cJSON *sampleEnv_field = NULL;
  const cJSON *composition_field = NULL;
  const cJSON *field_value = NULL;

  if (access(inputFile, F_OK|R_OK) == 0)
    fp = fopen(inputFile, "r");
  else {
    fprintf(stderr, "No such file or no access: %s\n", inputFile);
    return NULL;
  }
  fseek(fp, 0, SEEK_END);
  size_t file_length = ftell(fp);
  rewind(fp);

  char infile[file_length+1];
  fread(infile, sizeof(char), file_length, fp);
  fclose(fp);
  infile[file_length] = '\0';

  halir_workspace *ret_workspace = (halir_workspace*)malloc(sizeof(halir_workspace));

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
  project_field = cJSON_GetObjectItem(input_field, "project");
  field_value = cJSON_GetObjectItem(project_field, "pname");
  if (cJSON_IsString(field_value)) {
    strcpy(ret_workspace->pname, field_value->valuestring);
    /*printf("pname: %s\n", ret_workspace->pname);*/
  } else {
    fprintf(stderr, "Missing temp keyword or wrong type\n");
    read_err = 3;
    goto end;
  }
  field_value = cJSON_GetObjectItem(project_field, "rootDir");
  if (cJSON_IsString(field_value)) {
    strcpy(ret_workspace->rootDir, field_value->valuestring);
    /*printf("rootDir: %s\n", ret_workspace->rootDir);*/
  } else {
    fprintf(stderr, "Missing rootDir keyword or wrong type\n");
    read_err = 3;
    goto end;
  }
  field_value = cJSON_GetObjectItem(project_field, "pcomments");
  if (cJSON_IsString(field_value)) {
    strcpy(ret_workspace->pcomments, field_value->valuestring);
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
    ret_workspace->pfiles = malloc(arrayLength*sizeof(char*));
    ret_workspace->pfiles_length = arrayLength;
    cJSON_ArrayForEach(pfile, field_value) {
      if (cJSON_IsString(pfile)) {
        ret_workspace->pfiles[counter] = (char*)malloc(sizeof(pfile->valuestring));
        strcpy(ret_workspace->pfiles[counter], pfile->valuestring);
        /*printf("pfiles: %s\n", ret_workspace->pfiles[counter]);*/
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
  field_value = cJSON_GetObjectItem(sampleEnv_field, "temp");
  if (cJSON_IsNumber(field_value)) {
    ret_workspace->temp = field_value->valuedouble;
    /*printf("temp: %3.1f ", ret_workspace->temp);*/
  } else {
    fprintf(stderr, "Missing temp keyword or wrong type\n");
    read_err = 3;
    goto end;
  }
  field_value = cJSON_GetObjectItem(sampleEnv_field, "tempU");
  if (cJSON_IsString(field_value)) {
    strncpy(ret_workspace->tempU, field_value->valuestring, UNIT_LEN-1);
    ret_workspace->tempU[UNIT_LEN-1] = '\0';
    /*printf("%s\n", ret_workspace->tempU);*/
  } else {
    fprintf(stderr, "Missing tempU keyword or wrong type\n");
    read_err = 3;
    goto end;
  }
  field_value = cJSON_GetObjectItem(sampleEnv_field, "press");
  if (cJSON_IsNumber(field_value)) {
    ret_workspace->press = field_value->valuedouble;
    /*printf("press: %3.1f ", ret_workspace->press);*/
  } else {
    fprintf(stderr, "Missing press keyword or wrong type\n");
    read_err = 3;
    goto end;
  }
  field_value = cJSON_GetObjectItem(sampleEnv_field, "pressU");
  if (cJSON_IsString(field_value)) {
    strncpy(ret_workspace->pressU, field_value->valuestring, UNIT_LEN-1);
    ret_workspace->pressU[UNIT_LEN-1] = '\0';
    /*printf("%s\n", ret_workspace->pressU);*/
  } else {
    fprintf(stderr, "Missing pressU keyword or wrong type\n");
    read_err = 3;
    goto end;
  }
  field_value = cJSON_GetObjectItem(sampleEnv_field, "pathL");
  if (cJSON_IsNumber(field_value)) {
    ret_workspace->pathL = field_value->valuedouble;
    /*printf("pathL: %3.1f ", ret_workspace->pathL);*/
  } else {
    fprintf(stderr, "Missing pathL keyword or wrong type\n");
    read_err = 3;
    goto end;
  }
  field_value = cJSON_GetObjectItem(sampleEnv_field, "pathLU");
  if (cJSON_IsString(field_value)) {
    strncpy(ret_workspace->pathLU, field_value->valuestring, UNIT_LEN-1);
    ret_workspace->pathLU[UNIT_LEN-1] = '\0';
    /*printf("%s\n", ret_workspace->pathLU);*/
  } else {
    fprintf(stderr, "Missing pathLU keyword or wrong type\n");
    read_err = 3;
    goto end;
  }
  field_value = cJSON_GetObjectItem(sampleEnv_field, "bgfile");
  if (cJSON_IsString(field_value)) {
    strcpy(ret_workspace->bgfile, field_value->valuestring);
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
    cJSON_ArrayForEach(value, field_value) {
      if (cJSON_IsNumber(value)) {
        ret_workspace->ROI[counter] = value->valuedouble;
      }
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
    strcpy(apod, field_value->valuestring);
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
    strcpy(ftype, field_value->valuestring);
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
    ret_workspace->composition = malloc(arrayLength*sizeof(halir_compound*));
    halir_compound *current_comp;
    cJSON *comp;
    counter = 0;
    cJSON_ArrayForEach(comp, composition_field) {
      ret_workspace->composition[counter] = (halir_compound*)malloc(sizeof(halir_compound));
      current_comp = ret_workspace->composition[counter];
      current_comp->hitran_prms = NULL;
      value = cJSON_GetObjectItem(comp, "molec");
      if (cJSON_IsString(value)) {
        strcpy(current_comp->molec, value->valuestring);
        /*printf("molec: %s\n", current_comp->molec);*/
      } else {
        fprintf(stderr, "Missing molec keyword or wrong type\n");
        read_err = 4;
        goto end;
      }
      value = cJSON_GetObjectItem(comp, "isotop");
      if (cJSON_IsString(value)) {
        strcpy(current_comp->isotop, value->valuestring);
        /*printf("isotop: %s\n", current_comp->isotop);*/
      } else {
        fprintf(stderr, "Missing isotop keyword or wrong type\n");
        read_err = 4;
        goto end;
      }
      value = cJSON_GetObjectItem(comp, "pressU");
      if (cJSON_IsString(value)) {
        strcpy(current_comp->pressU, value->valuestring);
        /*printf("pressU: %s\n", current_comp->pressU);*/
      } else {
        fprintf(stderr, "Missing pressU keyword or wrong type\n");
        read_err = 4;
        goto end;
      }
      value = cJSON_GetObjectItem(comp, "vmr");
      if (cJSON_IsNumber(value)) {
        current_comp->vmr = value->valuedouble;
        /*printf("vmr: %e\n", current_comp->vmr);*/
      } else {
        fprintf(stderr, "Missing vmr keyword or wrong type\n");
        read_err = 4;
        goto end;
      }
      value = cJSON_GetObjectItem(comp, "prmfile");
      if (cJSON_IsString(value)) {
        strcpy(current_comp->prmfile, value->valuestring);
        /*printf("prmfile: %s\n", current_comp->prmfile);*/
      } else {
        fprintf(stderr, "Missing prmfile or wrong type\n");
        read_err = 4;
        goto end;
      }
      // Is there a file!
      if (strcmp(current_comp->prmfile, "") != 0) {
        FILE *fp;
        if (access(current_comp->prmfile, F_OK|R_OK) == 0)
          fp = fopen(current_comp->prmfile, "rb");
        else {
          fprintf(stderr, "Cant access the prmfile: %s\n", current_comp->prmfile);
          read_err = 4;
          goto end;
        }
        readSize = fread(&current_comp->hitran_head, sizeof(halir_HitranHead), 1, fp);
        if (readSize != 1) {
          fprintf(stderr, "Error while reading prmfile: %s\n", current_comp->prmfile);
          read_err = 5;
          goto end;
        }
        current_comp->hitran_head.molecs = (int*)malloc(sizeof(int)*current_comp->hitran_head.nisotp);
        readSize = fread(current_comp->hitran_head.molecs, sizeof(int), current_comp->hitran_head.nisotp, fp);
        if (readSize != current_comp->hitran_head.nisotp) {
          fprintf(stderr, "Error while reading prmfile: %s\n", current_comp->prmfile);
          read_err = 5;
          goto end;
        }
        current_comp->hitran_prms = (halir_HitranLine*)malloc(sizeof(halir_HitranLine)*current_comp->hitran_head.ndatapnts);
        readSize = fread(current_comp->hitran_prms, sizeof(halir_HitranLine), current_comp->hitran_head.ndatapnts, fp);
        if (readSize != current_comp->hitran_head.ndatapnts) {
          fprintf(stderr, "Error while reading prmfile: %s\n", current_comp->prmfile);
          read_err = 5;
          goto end;
        }

        /*printf("From hitran file\n");*/
        /*printf("molec: %s\n", current_comp->hitran_head.molec);*/
        /*printf("ROI: %f %f\n", current_comp->hitran_head.roi_low,current_comp->hitran_head.roi_high);*/
        /*printf("ndatapnts: %d\n", current_comp->hitran_head.ndatapnts);*/
      }
      /*for (int i=0; i < current_comp->hitran_head.ndatapnts; i++) {*/
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
  }

end:
  cJSON_Delete(root);
  if (read_err) {
    fprintf(stderr, "Something went wrong while reading input\n");
    return NULL;
  } else
    return ret_workspace;
}
