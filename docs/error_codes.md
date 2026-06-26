# HALIR Error Code Matrix

This page summarizes the current return-code conventions used by the public API.
Use it as a quick reference while integrating HALIR.

## General conventions

- `0`: Success.
- Non-zero: Failure.
- `1`: Input/state validation failure in most APIs.
- `2`: Conversion, allocation, or file/data-loading failure in APIs that perform those operations.

## Function matrix

| Function | `0` | `1` | `2` |
|---|---|---|---|
| `halir_workspace_set_project` | Success | Invalid argument (NULL) | Bounded string copy failed |
| `halir_workspace_add_project_file` | Success | Invalid argument | Allocation failed |
| `halir_workspace_set_sample_env` | Success | Validation failure | Unit conversion or string assignment failure |
| `halir_workspace_add_composition` | Success | Invalid argument | Allocation or bounded copy failure |
| `halir_compound_set_vmr` | Success | Invalid workspace/index/range | N/A |
| `halir_compound_set_concentration` | Success | Invalid workspace/index/value | Unit conversion failure |
| `halir_compound_load_prmfile` | Success | Invalid argument | File access/read/format loading failure |
| `halir_workspace_validate` | Valid workspace | Validation failure | N/A |
| `halir_test_calc` | Success | Invalid state/runtime setup failure | N/A |

## Notes

- Codes are intentionally compact and function-local; meaning is determined by the API being called.
- For user-facing diagnostics, prefer checking the specific API and emitting a context-rich message.
- Keep cleanup deterministic: call `halir_workspace_free` on all failure paths after workspace creation.
