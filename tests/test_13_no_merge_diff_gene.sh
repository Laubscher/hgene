#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]:-$0}" )" &> /dev/null && pwd )"
SRC="$SCRIPT_DIR/../src"
DATA="$SCRIPT_DIR/data"
DB="$SCRIPT_DIR/db"
FASTA="$DB/HHV1.fasta"
GFF="$DB/HHV1.gff"
VCF_IN="$DATA/test_13_hsv1_XX.bcf.vcf"
BAM="$DATA/test_13_hsv1_XX.bam"
BAI="$DATA/test_13_hsv1_XX.bam.bai"
OUT="test_13.out.vcf"
LOG="test_13.err.log"
cleanup() { rm -f "$OUT" "$LOG"; }
trap cleanup EXIT
[[ -s "$FASTA" ]] || { echo "FAIL: missing FASTA"; exit 1; }
[[ -s "$GFF" ]] || { echo "FAIL: missing GFF"; exit 1; }
[[ -s "$VCF_IN" ]] || { echo "FAIL: missing VCF"; exit 1; }
[[ -s "$BAM" ]] || { echo "FAIL: missing BAM"; exit 1; }
[[ -s "$BAI" ]] || samtools index "$BAM"
python3 "$SRC/hgene_codon_haplotype_merge.py" "$VCF_IN" "$FASTA" "$GFF" "$BAM" > "$OUT" 2> "$LOG"
! grep -q "MERGED_POS=" "$OUT" || { echo "FAIL: unexpected merge"; cat "$OUT"; exit 1; }
echo "OK test_13_no_merge_diff_gene"
