#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]:-$0}" )" &> /dev/null && pwd )"
SRC="$SCRIPT_DIR/../src/"
DATA="$SCRIPT_DIR/data"
DB="$SCRIPT_DIR/db"
OUT="test_14.out.vcf"

cleanup() { rm -f "$OUT"; }
trap cleanup EXIT

python3 "$SRC/merge_codon_mutations.py" \
  "$DATA/test14_adjacent_diff_codons.vcf" \
  "$DB/TEST.fasta" \
  "$DB/TEST.gff" 0.5 > "$OUT"

count="$(awk '!/^#/ {c++} END{print c+0}' "$OUT")"
[[ "$count" -eq 2 ]] || { echo "FAIL: expected 2 variant lines (no merge), got $count"; cat "$OUT"; exit 1; }

# Aucun MERGED_POS attendu
grep -Fq 'MERGED_POS=' "$OUT" && { echo "FAIL: found MERGED_POS but should not merge"; cat "$OUT"; exit 1; }

echo "OK test_14_no_merge_adjacent_codons"
