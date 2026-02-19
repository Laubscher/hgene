#!/usr/bin/env bash
set -euo pipefail

echo "Running hgene tests ..."

TEST_DIR="tests"

if [[ ! -d "$TEST_DIR" ]]; then
  echo "No tests directory found."
  exit 1
fi

fail=0

for t in "$TEST_DIR"/*.sh; do
  echo "--------------------------------"
  echo "Running $t"
  if bash "$t"; then
    echo "PASS: $t"
  else
    echo "FAIL: $t"
    fail=1
  fi
done

echo "--------------------------------"

if [[ $fail -eq 0 ]]; then
  echo "All tests passed."
else
  echo "Some tests failed."
  exit 1
fi

