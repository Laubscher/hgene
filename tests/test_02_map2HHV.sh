#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]:-$0}" )" &> /dev/null && pwd )"
TEST_NAME="test_02_hsv1_XX"
FASTQ="${TEST_NAME}.trimmed.fastq"
SAM="${TEST_NAME}.sam"

cleanup() { rm -f "$FASTQ" "$SAM"; }
trap cleanup EXIT

cp "$SCRIPT_DIR/data/$FASTQ" .

set +e
bash "$SCRIPT_DIR/../src/map2HHV.sh" $FASTQ HHV1 2
rc=$?
set -e

if [[ $rc -ne 0 ]]; then
  echo "OK test_02_map2HHV: script a échoué (attendu pour chimérique)"
  exit 0
fi

[[ -f "$SAM" ]] || { echo "FAIL: SAM not generated"; exit 1; }

# KO si un read passe
mapped="$(samtools view -S -c -F 4 "$SAM")"
[[ "$mapped" -eq 0 ]] || { echo "FAIL: $mapped read(s) ont passé"; exit 1; }

echo "OK test_02_map2HHV: 0 read a passé"

