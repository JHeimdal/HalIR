#include <stdlib.h>

#include <HalIR/halir.h>

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
halir_result_create(halir_simulation_setup *workspace, size_t nspectra)
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
