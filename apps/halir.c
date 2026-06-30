#include "HalIR/halir.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

static const char *HALIR_EXAMPLE_INPUT =
  "{\n"
  "  \"input\": {\n"
  "    \"project\": {\n"
  "      \"pname\": \"ExampleProject\",\n"
  "      \"rootDir\": \"./output\",\n"
  "      \"hapi_db\": \"\",\n"
  "      \"pcomments\": \"Example HalIR input\",\n"
  "      \"pfiles\": []\n"
  "    },\n"
  "    \"sampleEnv\": {\n"
  "      \"temp\": 303.0,\n"
  "      \"tempU\": \"K\",\n"
  "      \"press\": 1.0,\n"
  "      \"pressU\": \"atm\",\n"
  "      \"pathL\": 300.0,\n"
  "      \"pathLU\": \"m\",\n"
  "      \"ROI\": [2172.0, 2174.0],\n"
  "      \"res\": 0.1,\n"
  "      \"apod\": \"Boxcar\",\n"
  "      \"fov\": 0.1,\n"
  "      \"ftype\": \"Transmission\",\n"
  "      \"bgfile\": \"\"\n"
  "    },\n"
  "    \"composition\": [\n"
  "      {\n"
  "        \"molec\": \"CO\",\n"
  "        \"isotop\": \"Natural\",\n"
  "        \"vmr\": 5.921539600296e-05,\n"
  "        \"prmfile\": \"\"\n"
  "      }\n"
  "    ]\n"
  "  }\n"
  "}\n";

static void print_usage(const char *program_name)
{
  fprintf(stderr,
          "Usage: %s [--example] [--spc <output.spc>] <input.json>\n"
          "  --example, -e   Print an example input file to stdout\n"
          "  --spc, -s       Write calculated result to an SPC file\n"
          "  --help, -h      Show this help message\n",
          program_name);
}

int main(int argc, char **argv)
{
  const char *inputFile = NULL;
  const char *spcOutputFile = NULL;
  halir_result *result;

  for (int argi = 1; argi < argc; argi++)
  {
    if ((strcmp(argv[argi], "--help") == 0) || (strcmp(argv[argi], "-h") == 0))
    {
      print_usage(argv[0]);
      return 0;
    }

    if ((strcmp(argv[argi], "--example") == 0) || (strcmp(argv[argi], "-e") == 0))
    {
      fputs(HALIR_EXAMPLE_INPUT, stdout);
      return 0;
    }

    if ((strcmp(argv[argi], "--spc") == 0) || (strcmp(argv[argi], "-s") == 0))
    {
      argi++;
      if (argi >= argc)
      {
        fprintf(stderr, "Missing output path for --spc option.\n");
        print_usage(argv[0]);
        return 99;
      }
      spcOutputFile = argv[argi];
      continue;
    }

    if (argv[argi][0] == '-')
    {
      fprintf(stderr, "Unknown option: %s\n", argv[argi]);
      print_usage(argv[0]);
      return 99;
    }

    if (inputFile != NULL)
    {
      fprintf(stderr, "Only one input file may be provided.\n");
      print_usage(argv[0]);
      return 99;
    }

    inputFile = argv[argi];
  }

  if (inputFile == NULL)
  {
    print_usage(argv[0]);
    return 99;
  }

  halir_simulation_setup *work = halir_parseJSONinput(inputFile);
  if (work == NULL)
    return 99;
  halir_print_simulation_setup(work);

  result = halir_calculate_result(work);
  if (result == NULL)
  {
    halir_simulation_setup_free(work);
    return 1;
  }

  if (spcOutputFile != NULL)
  {
    if (halir_result_write_spc(result, spcOutputFile) != 0)
    {
      fprintf(stderr, "Failed to write SPC output: %s\n", spcOutputFile);
      halir_result_free(result);
      halir_simulation_setup_free(work);
      return 1;
    }

    halir_result_free(result);
    halir_simulation_setup_free(work);
    return 0;
  }

  for (size_t comp = 0; comp < result->nspectra; comp++)
  {
    for (size_t i = 0; i < result->spectra[comp].ndatapnts; i++)
    {
      printf("%f %f\n", result->spectra[comp].wavenum[i], result->spectra[comp].data[i]);
    }
  }

  halir_result_free(result);
  halir_simulation_setup_free(work);
  return 0;
}
