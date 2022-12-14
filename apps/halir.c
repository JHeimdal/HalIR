#include "HalIR/halir.h"
#include "apps/cJSON.h"
#include <stdio.h>
#include <unistd.h>

int parseJSONinput(const char* const metadata)
{
  size_t counter;
  const cJSON *input_field = NULL;
  const cJSON *project_field = NULL;
  const cJSON *sampleEnv_field = NULL;
  const cJSON *field_value = NULL;

  cJSON *root = cJSON_Parse(metadata);
  if (root == NULL) {
    const char *err_ptr = cJSON_GetErrorPtr();
    if (err_ptr != NULL)
      fprintf(stderr, "Error before: %s\n", err_ptr);
    goto end;
  }
  input_field = cJSON_GetObjectItem(root, "input");
  project_field = cJSON_GetObjectItem(input_field, "project");
  field_value = cJSON_GetObjectItem(project_field, "pname");
  if (cJSON_IsString(field_value)) {
    printf("pname: %s\n", field_value->valuestring);
  }
  field_value = cJSON_GetObjectItem(project_field, "rootDir");
  if (cJSON_IsString(field_value)) {
    printf("rootDir: %s\n", field_value->valuestring);
  }
  field_value = cJSON_GetObjectItem(project_field, "pcomments");
  if (cJSON_IsString(field_value)) {
    printf("pcomments: %s\n", field_value->valuestring);
  }
  field_value = cJSON_GetObjectItem(project_field, "pfiles");
  counter = 0;
  cJSON *pfile = NULL;
  if (cJSON_IsArray(field_value)) {
    cJSON_ArrayForEach(pfile, field_value) {
      if (cJSON_IsString(pfile)) {
        printf("pfiles: %s\n", pfile->valuestring);
      }
    }
  }

  // Read the sampleEnv field
  sampleEnv_field = cJSON_GetObjectItem(input_field, "sampleEnv");
  field_value = cJSON_GetObjectItem(sampleEnv_field, "temp");
  if (cJSON_IsNumber(field_value)) {
    printf("temp: %3.1f ", field_value->valuedouble);
  }
  field_value = cJSON_GetObjectItem(sampleEnv_field, "tempU");
  if (cJSON_IsString(field_value)) {
    printf("%s\n", field_value->valuestring);
  }
  field_value = cJSON_GetObjectItem(sampleEnv_field, "press");
  if (cJSON_IsNumber(field_value)) {
    printf("press: %3.1f ", field_value->valuedouble);
  }
  field_value = cJSON_GetObjectItem(sampleEnv_field, "pressU");
  if (cJSON_IsString(field_value)) {
    printf("%s\n", field_value->valuestring);
  }
  field_value = cJSON_GetObjectItem(sampleEnv_field, "pathL");
  if (cJSON_IsNumber(field_value)) {
    printf("pathL: %3.1f ", field_value->valuedouble);
  }
  field_value = cJSON_GetObjectItem(sampleEnv_field, "pathLU");
  if (cJSON_IsString(field_value)) {
    printf("%s\n", field_value->valuestring);
  }
  field_value = cJSON_GetObjectItem(sampleEnv_field, "ROI");
  counter = 0;
  cJSON *value = NULL;
  if (cJSON_IsArray(field_value)) {
    printf("ROI: [");
    cJSON_ArrayForEach(value, field_value) {
      if (cJSON_IsNumber(value)) {
        printf("%f ", value->valuedouble);
      }
    }
    printf("]\n");
  }
end:
  cJSON_Delete(root);
  return 0;
}

int main()
{
  FILE *fp;
  char inputFile[] = "input.json";
  if (access(inputFile, F_OK|R_OK) == 0)
    fp = fopen(inputFile, "r");
  else
    return 9;
  fseek(fp, 0, SEEK_END);
  size_t file_length = ftell(fp);
  rewind(fp);

  char infile[file_length+1];
  fread(infile, sizeof(char), file_length, fp);
  fclose(fp);
  infile[file_length] = '\0';
  /*printf("%s\n", infile);*/
  parseJSONinput(infile);
  return 0;
}
