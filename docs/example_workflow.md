# HALIR API Workflow Example

This page shows a minimal direct-API flow for configuring a workspace,
loading spectral line parameters, validating inputs, running the calculation,
and releasing resources.

## Steps

1. Create workspace.
2. Set project metadata and sample environment.
3. Add composition and set concentration (or VMR).
4. Load HITRAN prmfile for the composition.
5. Validate workspace.
6. Run calculation.
7. Free workspace.

## Example

```c
#include "HalIR/halir.h"

int run_example(void)
{
  halir_simulation_setup *work = NULL;
  size_t comp_idx = 0;
  int rc = 0;

  work = halir_simulation_setup_create();
  if (work == NULL) {
    return 1;
  }

  rc = halir_simulation_setup_set_project(work, "ExampleProject", "/tmp", "API example");
  if (rc != 0) {
    halir_simulation_setup_free(work);
    return rc;
  }

  rc = halir_simulation_setup_set_sample_env(
      work,
      296.0, K,
      1.0, ATM,
      100.0, CM,
      2000.0, 2245.0,
      0.1, 0.1,
      HALIR_BOXCAR,
      HALIR_TRANSMISSION,
      "");
  if (rc != 0) {
    halir_simulation_setup_free(work);
    return rc;
  }

  rc = halir_simulation_setup_add_composition(work, "CO", "Natural", &comp_idx);
  if (rc != 0) {
    halir_simulation_setup_free(work);
    return rc;
  }

  rc = halir_compound_set_vmr(work, comp_idx, 5.921539600296e-05);
  if (rc != 0) {
    halir_simulation_setup_free(work);
    return rc;
  }

  rc = halir_compound_load_prmfile(work, comp_idx, "tests/CO_test/CO.hpar");
  if (rc != 0) {
    halir_simulation_setup_free(work);
    return rc;
  }

  rc = halir_simulation_setup_validate(work);
  if (rc != 0) {
    halir_simulation_setup_free(work);
    return rc;
  }

  rc = halir_test_calc(work);
  halir_simulation_setup_free(work);
  return rc;
}
```

## Notes

- `halir_simulation_setup_set_sample_env` normalizes core environment units to K, ATM, and CM.
- `halir_simulation_setup_validate` should be called before `halir_test_calc`.
- Always call `halir_simulation_setup_free` on every success and failure path.
