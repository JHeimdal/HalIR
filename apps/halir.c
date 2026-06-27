#include "HalIR/halir.h"
#include <stdio.h>
//#include <string.h>
#include <unistd.h>
#include <stdlib.h>


int main(int argc, char **argv)
{
  char *inputFile;
  halir_result *result;
  if (argc > 1)
    inputFile = argv[1];
  else {
    printf("Usage %s input.json\n", argv[0]);
    exit(99);
  }
  // char inputFile[] = "input.json";

  halir_workspace *work = halir_parseJSONinput(inputFile);
  if (work == NULL)
    return 99;
  halir_print_workspace(work);

  result = halir_calculate_result(work);
  if (result == NULL) {
    halir_workspace_free(work);
    return 1;
  }

  for (size_t comp = 0; comp < result->nspectra; comp++) {
    for (size_t i = 0; i < result->spectra[comp].ndatapnts; i++) {
      printf("%f %f\n", result->spectra[comp].wavenum[i], result->spectra[comp].data[i]);
    }
  }

  halir_result_free(result);
  halir_workspace_free(work);
  return 0;
}
