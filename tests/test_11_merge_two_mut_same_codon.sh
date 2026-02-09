#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]:-$0}" )" &> /dev/null && pwd )"
SRC="$SCRIPT_DIR/../src/"
DATA="$SCRIPT_DIR/data"
DB="$SCRIPT_DIR/db"
OUT="test_11.out.vcf"

cleanup() { rm -f "$OUT"; }
trap cleanup EXIT

python3 "$SRC/merge_codon_mutations.py" \
  "$DATA/test11_two_mut_same_codon.vcf" \
  "$DB/TEST.fasta" \
  "$DB/TEST.gff" 0.5 > "$OUT"

# 1) Doit y avoir exactement 1 variant (2 -> 1)
count="$(awk '!/^#/ {c++} END{print c+0}' "$OUT")"
[[ "$count" -eq 1 ]] || { echo "FAIL: expected 1 variant line, got $count"; cat "$OUT"; exit 1; }

# 2) Doit contenir le merge et le REF/ALT codon-level attendus
grep -Fq $'TEST-A\t4\t.\tAAT\tGGT\t' "$OUT" || { echo "FAIL: merged REF/ALT not as expected (POS=4 REF=AAT ALT=GGT)"; cat "$OUT"; exit 1; }

# 3) INFO MERGED_POS
grep -Fq 'MERGED_POS=4,5' "$OUT" || { echo "FAIL: MERGED_POS missing or wrong"; cat "$OUT"; exit 1; }

# 4) BCSQ doit être mis à jour (AA + nuc)
grep -Fq 'BCSQ=missense|geneA|geneA|protein_coding|+|2N>2G|4AAT>GGT' "$OUT" || {
  echo "FAIL: BCSQ not updated as expected"; cat "$OUT"; exit 1;
}

echo "OK test_11_merge_two_mut_same_codon"
