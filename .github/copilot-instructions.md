# Copilot Instructions for HalIR

## Project Overview
- HalIR is a mixed-language scientific codebase for IR spectrum analysis and modeling.
- Primary implementation languages are C, C++, and Fortran.
- Core library target: HalIR (shared library) defined in src/CMakeLists.txt.
- Main executable target: halir in apps/CMakeLists.txt.
- Core parser and numerical kernel currently live in src/halir.c.

## Build and Test Workflow
- Configure and build with CMake:
  - cmake -S . -B Build
  - cmake --build Build
- Run tests with CTest:
  - ctest --test-dir Build --output-on-failure
- Always build and run tests after non-trivial edits.

## Coding Priorities
- Preserve scientific correctness and numerical stability over stylistic refactors.
- Prefer small, focused changes in existing files over large rewrites.
- Keep C source valid C (do not introduce C++ syntax in C files).
- Keep public API changes minimal and explicit.

## Memory and Resource Safety
- Use halir_workspace_free for cleanup of objects returned by halir_parseJSONinput.
- Ensure every opened file handle is closed on all success and failure paths.
- For parser code, keep deterministic cleanup behavior (single cleanup path is preferred).
- Avoid unchecked string copies into fixed-size buffers.

## Parser and Input Handling
- JSON parsing behavior is centralized in halir_parseJSONinput.
- Validate field types before use.
- Reject malformed input with clear stderr messages and NULL return.
- Preserve unit-conversion semantics used by halir_Units_to_Hitran.
- Keep ROI validation strict (numeric and expected shape).

## HITRAN prmfile Handling
- Preserve binary layout assumptions for:
  - halir_HitranHead
  - halir_HitranLine
- Treat prmfile input as untrusted:
  - validate metadata
  - check fread sizes for each section
  - fail cleanly on truncated/corrupt files
- Add/maintain regression tests for prmfile failure modes.

## Numerical Kernel Guardrails
- For spectral-grid and line-shape loops, guard against:
  - non-finite values
  - non-positive step sizes
  - out-of-bounds index ranges
- Initialize accumulation buffers before summation.
- Keep loops bounds-safe near edge indices.

## Tests and Regression Expectations
- If behavior changes in parsing or prmfile I/O, add/update tests in:
  - tests/test_input.cpp
  - tests/test_hpar.cpp
- Keep tests deterministic and self-contained.
- Prefer constructing temporary corrupted fixtures in tests rather than modifying committed binary fixtures.

## CMake and Dependency Notes
- Required dependencies include GSL and cerf.
- Fortran support is enabled in src/CMakeLists.txt and is part of library build.
- Do not remove or bypass Fortran build settings without a clear replacement plan.

## Style and Change Discipline
- Follow existing project style and naming patterns.
- Do not reformat unrelated code.
- Keep comments short and only where logic is non-obvious.
- For bugfixes, include the smallest test that reproduces and prevents regression.
