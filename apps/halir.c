#include "HalIR/halir.h"
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>


int main(int argc, char **argv)
{
  char *inputFile;
  if (argc > 1)
      inputFile = argv[1];
  // char inputFile[] = "input.json";

  halir_workspace *work = halir_parseJSONinput(inputFile);
  if (work == NULL)
    return 99;
  halir_print_workspace(work);
  return 0;
}
