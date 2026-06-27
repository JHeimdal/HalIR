#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

PRESET="${1:-debug}"
BUILD_DIR="${2:-Debug}"
FIXTURE_PATH="$SCRIPT_DIR/fixtures/calc_result_metrics.txt"
TMP_PATH="${FIXTURE_PATH}.tmp"

cd "$REPO_ROOT"

cmake --preset "$PRESET"
cmake --build "$BUILD_DIR" --target test_calc

DUMP_OUTPUT="$($REPO_ROOT/$BUILD_DIR/tests/test_calc --dump)"

PARSED="$(printf '%s\n' "$DUMP_OUTPUT" | awk '
  BEGIN { count = 0 }
  /^spec=/ {
    delete kv
    for (i = 1; i <= NF; i++) {
      split($i, pair, "=")
      kv[pair[1]] = pair[2]
    }
    if (!("n" in kv) || !("w0" in kv) || !("wN" in kv) ||
        !("sum" in kv) || !("max" in kv) || !("mid" in kv)) {
      print "Failed to parse dump line: " $0 > "/dev/stderr"
      exit 2
    }
    printf "%s %s %s %s %s %s\n", kv["n"], kv["w0"], kv["wN"], kv["sum"], kv["max"], kv["mid"]
    count++
  }
  END {
    if (count == 0) {
      print "No spec lines found in test_calc --dump output" > "/dev/stderr"
      exit 3
    }
  }
')"

printf '%s\n' "$PARSED" > "$TMP_PATH"
mv "$TMP_PATH" "$FIXTURE_PATH"

echo "Updated fixture: $FIXTURE_PATH"
