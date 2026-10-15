#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]:-$0}" )" &> /dev/null && pwd )"
TEST_NAME="test_03_hsv1_XX"
FASTQ="${TEST_NAME}.fastq"


cleanup() {
    rm -f "$FASTQ" "${TEST_NAME}.log"
    rm -rf "${TEST_NAME}_output"
}
trap cleanup EXIT

cp "$SCRIPT_DIR/data/$FASTQ" .

bash "$SCRIPT_DIR/../src/hgene" -v HHV1 "$TEST_NAME" 2

VCF="${TEST_NAME}_output/${TEST_NAME}.vcf.gz"

[[ -f "$VCF" ]] || { echo "FAIL: VCF not generated"; exit 1; }

VCF_CONTENT=$(gunzip -c -- "$VCF")
EXPECTED_UL23='BCSQ=missense|UL23-HHV1|UL23-HHV1|protein_coding|+|10A>10V|29C>T'
EXPECTED_UL30='BCSQ=missense|UL30-HHV1|UL30-HHV1|protein_coding|+|10S>10F|29C>T'

[[ "$VCF_CONTENT" == *"$EXPECTED_UL23"* ]] || {
    echo "FAIL: expected BCSQ annotation for UL23 is missing"
    exit 1
}

[[ "$VCF_CONTENT" == *"$EXPECTED_UL30"* ]] || {
    echo "FAIL: expected BCSQ annotation for UL30 is missing"
    exit 1
}

echo "OK test_03_expected"
