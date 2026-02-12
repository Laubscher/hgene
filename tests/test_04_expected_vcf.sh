#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]:-$0}" )" &> /dev/null && pwd )"
TEST_NAME="test_04_hsv1_XX"
FASTQ="${TEST_NAME}.fastq"


cleanup() { rm -f "$FASTQ"; rm -r "${TEST_NAME}_output" ; }
trap cleanup EXIT

cp "$SCRIPT_DIR/data/$FASTQ" .

bash "$SCRIPT_DIR/../src/hgene" -v HHV1 $TEST_NAME 2

VCF="${TEST_NAME}_output/${TEST_NAME}.vcf.gz"

[[ -f "$VCF" ]] || { echo "FAIL: VCF not generated"; exit 1; }

UL23=$(gunzip -c $VCF | tail -2 | head -1 | cut -f5 -d ";" | tr -d '\n')

UL30=$(gunzip -c $VCF | tail -1 | cut -f5 -d ";" | tr -d '\n')

[[ "$(echo "$UL23")" != "BCSQ=missense|UL23-HHV1|UL23-HHV1|protein_coding|+|10A>10V|29C>T
" ]] || { echo "FAIL: VCF contains not expected BCSQ for UL23"; exit 1; }


[[ "$UL30" != "BCSQ=missense|UL30-HHV1|UL30-HHV1|protein_coding|+|10S>10F|29C>T

" ]] || { echo "FAIL: VCF contains not expected BCSQ for UL30"; exit 1; }

echo "OK test_04_expected"

rm test_04_hsv1_XX.log
