#include "HalIR/halir.h"
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>


int main(int argc, char **argv)
{
  char *inputFile;
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
  /*halir_test_calc(work);*/
  return 0;
}
