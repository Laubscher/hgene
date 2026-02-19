#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]:-$0}" )" &> /dev/null && pwd )"
TEST_NAME="test_02_hsv1_XX"
FASTQ="${TEST_NAME}.trimmed.fastq"
SAM="${TEST_NAME}.sam"

cleanup() { rm -f -- "$FASTQ" "$SAM"; }
trap cleanup EXIT

[[ -f "$SCRIPT_DIR/data/$FASTQ" ]] || { echo "FAIL: missing test FASTQ: $SCRIPT_DIR/data/$FASTQ"; exit 1; }
cp -- "$SCRIPT_DIR/data/$FASTQ" .

set +e
bash "$SCRIPT_DIR/../src/hgene_map.sh" "$FASTQ" HHV1 2
rc=$?
set -e

[[ -s "$SAM" ]] || { echo "FAIL: SAM not generated (or empty)"; exit 1; }

mapped="$(samtools view -c -F 4 "$SAM")"
[[ "$mapped" -eq 0 ]] || { echo "FAIL: $mapped mapped read(s) passed"; exit 1; }

echo "OK $TEST_NAME: 0 mapped reads"
