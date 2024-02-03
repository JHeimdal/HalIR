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

#include "HalIR/halir.h"
#include "HalIR/tips.h"
#include "HalIR/cJSON.h"
#include <cerf.h>

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
          printf("concU: %s\n", halir_Units_to_str[current_comp->concU]);
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
        /*printf("reading prmfile: %s with %d number of lines\n", current_comp->prmfile, current_comp->hitran_head.ndatapnts);*/
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

int halir_test_calc(halir_workspace *work)
{
  double q296, qT;
  double T = work->temp;
  double sqT = sqrtf(T);
  double P = work->press;
  double vc, S, q, tfac, mu_step, mu_off1;
  double sig_v = 0.1;
  float alphaD_max, alphaD_min, alphaL_max,  alphaL_min;

  const double kb_si = 1.3806503e-23; // J/K
  const double kb_erg = 1.3806503e-16; // erg/K
  const double invkb_erg = 1./kb_erg; // K/erg
  const double invkb_si = 1./kb_si; // K/J
  const double c_si = 299792458.; // m/s
  const double c_erg = 29979245800.; // cm/s
  const double atmmass_si=1.6605e-27; // g
  const double atmmass_erg=1.6605e-24; // g
  const double hc_erg = 6.62607015e-27*c_erg;
  const double hc_si = 6.62607015e-34*c_si;
  const double c1_erg = hc_erg/kb_erg;
  const double NL = 2.686780524e19;
  const double pi = 3.14159265358979;
  const double ln2 = log(2);
  const double sqrt_pi = sqrt(pi);
  const double sqrt_ln2 = sqrt(ln2);

  for (int comp=0; comp < work->composition_length; comp++) {
    halir_HitranHead *head = &work->composition[comp]->hitran_head;
    halir_HitranLine *prms = work->composition[comp]->hitran_prms;

    q = work->composition[comp]->vmr;
    tfac = sqrtf(2*ln2*kb_si*T);
    // Allocate vectors used in the calculations
    gsl_vector_float *v0 = gsl_vector_float_alloc(head->ndatapnts);
    /*gsl_vector_float *t_mu = gsl_vector_float_alloc(head->ndatapnts);*/
    gsl_vector_float *p_S  = gsl_vector_float_alloc(head->ndatapnts);
    gsl_vector_float *a_B  = gsl_vector_float_alloc(head->ndatapnts);
    /*gsl_vector_float *s_B  = gsl_vector_float_alloc(head->ndatapnts);*/
    gsl_vector_float *MM   = gsl_vector_float_alloc(head->ndatapnts);
    gsl_vector_float *taB  = gsl_vector_float_alloc(head->ndatapnts);
    gsl_vector_float *alphaD  = gsl_vector_float_alloc(head->ndatapnts);
    gsl_vector_float *alphaL  = gsl_vector_float_alloc(head->ndatapnts);

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
    mu_off1 = ceilf(50*GSL_MAX(alphaD_max, alphaL_max));

    /*printf("mu_step: %f mu_off1: %f\n", mu_step, mu_off1);*/
    size_t mu_size = (size_t)ceilf(((work->ROI[1]-work->ROI[0])+2*mu_off1)/mu_step);
    mu_step = ((work->ROI[1]-work->ROI[0])+2*mu_off1)/(mu_size-1);
    gsl_vector_float *mu = gsl_vector_float_alloc(mu_size);
    float start_value = work->ROI[0]-mu_off1;
    for (int i=0; i < mu_size; i++) {
      gsl_vector_float_set(mu, i, start_value+i*mu_step);
    }
    size_t off1 = (size_t)ceilf(mu_off1/mu_step);
    /*printf("mu_step: %f mu_size: %ld\n", mu_step, mu_size);*/
    /*printf("mu[0]: %f mu[-1]: %f\n", gsl_vector_float_get(mu, 0), gsl_vector_float_get(mu, mu_size-1));*/

    gsl_vector_float *y = gsl_vector_float_alloc(mu_size);
    // Print some information (possible in later callback?)
    /*printf("alphaL_min: %f alphaL_max %f\n", alphaL_min, alphaL_max);*/
    /*printf("alphaD_min: %f alphaD_max %f\n", alphaD_min, alphaD_max);*/
    for (int mm=0; mm < head->ndatapnts; mm++) {
          vc = gsl_vector_float_get(v0, mm);
          size_t idx = find_nearest_index(mu, vc);
          q296 = tips_2020(prms[mm].molec_num, prms[mm].isotp_num, 296);
          qT = tips_2020(prms[mm].molec_num, prms[mm].isotp_num, T);
          S = prms[mm].line_I * q296/qT * exp(c1_erg*prms[mm].low_state_en/T)/exp(c1_erg*prms[mm].low_state_en/296)*((1-exp(-c1_erg*vc/T))/(1-exp(-c1_erg*vc/296)));
          float a_D = gsl_vector_float_get(alphaD, mm);
          float a_L = gsl_vector_float_get(alphaL, mm);

          for (int i=(idx-off1); i < (idx+off1); i++) {
            // y->data[i] += S/a_D * q * NL * work->pathL * 298/T * voigt(sqrt_ln2*(mu->data[i]-vc)/a_D, 1/a_D, sqrt_ln2*a_L/a_D);
            y->data[i] += sqrt_ln2*S/(sqrt_pi*a_D) * q * work->pathL/100 * P*101325/1e4 * invkb_si/T * re_w_of_z(sqrt_ln2*(mu->data[i]-vc)/a_D, sqrt_ln2*a_L/a_D);
          }
    }
    for (int i = 0; i < mu_size; i++) {
      printf("%f %f\n", mu->data[i], y->data[i]);
    }
    gsl_vector_float_free(v0);
    gsl_vector_float_free(p_S);
    gsl_vector_float_free(a_B);
    gsl_vector_float_free(alphaL);
    gsl_vector_float_free(alphaD);
    gsl_vector_float_free(MM);
    gsl_vector_float_free(taB);
  }
  return 0;
}
