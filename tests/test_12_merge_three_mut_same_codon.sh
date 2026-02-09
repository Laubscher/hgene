#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]:-$0}" )" &> /dev/null && pwd )"
SRC="$SCRIPT_DIR/../src/"
DATA="$SCRIPT_DIR/data"
DB="$SCRIPT_DIR/db"
OUT="test_12.out.vcf"

cleanup() { rm -f "$OUT"; }
trap cleanup EXIT

python3 "$SRC/merge_codon_mutations.py" \
  "$DATA/test12_three_mut_same_codon.vcf" \
  "$DB/TEST.fasta" \
  "$DB/TEST.gff" 0.5 > "$OUT"

count="$(awk '!/^#/ {c++} END{print c+0}' "$OUT")"
[[ "$count" -eq 1 ]] || { echo "FAIL: expected 1 variant line, got $count"; cat "$OUT"; exit 1; }

# GGC -> TAA, POS=7
grep -Fq $'TEST-A\t7\t.\tGGC\tTAA\t' "$OUT" || { echo "FAIL: merged REF/ALT not as expected (POS=7 REF=GGC ALT=TAA)"; cat "$OUT"; exit 1; }

grep -Fq 'MERGED_POS=7,8,9' "$OUT" || { echo "FAIL: MERGED_POS missing or wrong"; cat "$OUT"; exit 1; }

# BCSQ doit passer en stop_gained et AA 3G>3*
grep -Fq 'BCSQ=stop_gained|geneA|geneA|protein_coding|+|3G>3*|7GGC>TAA' "$OUT" || {
  echo "FAIL: BCSQ not updated as expected for stop_gained"; cat "$OUT"; exit 1;
}

echo "OK test_12_merge_three_mut_same_codon"
