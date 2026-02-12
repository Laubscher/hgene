#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]:-$0}" )" &> /dev/null && pwd )"
TEST_NAME="test_01_hsv1_XX"
FASTQ="${TEST_NAME}.trimmed.fastq"
SAM="${TEST_NAME}.sam"

cleanup() { rm -f "$FASTQ" "$SAM"; }
trap cleanup EXIT

cp "$SCRIPT_DIR/data/$FASTQ" .

bash "$SCRIPT_DIR/../src/hgene_map.sh" "$FASTQ" HHV1 2

[[ -f "$SAM" ]] || { echo "FAIL: SAM not generated"; exit 1; }

# Compte les reads mappés (flag 0x4 = unmapped)
mapped="$(samtools view -S -c -F 4 "$SAM")"
[[ "$mapped" -gt 0 ]] || { echo "FAIL: SAM contains no mapped alignments"; exit 1; }

echo "OK test_01_map2HHV (SAM non-empty, mapped alignments present: $mapped)"

