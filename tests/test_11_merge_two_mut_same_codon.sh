#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]:-$0}" )" &> /dev/null && pwd )"

SRC="$SCRIPT_DIR/../src"
DATA="$SCRIPT_DIR/data"
DB="$SCRIPT_DIR/db"

FASTA="$DB/HHV1.fasta"
GFF="$DB/HHV1.gff"

# Inputs (dans tests/data)
VCF_IN="$DATA/test_11_hsv1_XX.bcf.vcf"
BAM="$DATA/test_11_hsv1_XX.bam"
BAI="$DATA/test_11_hsv1_XX.bam.bai"

OUT="test_11.out.vcf"
LOG="test_11.err.log"

cleanup() { rm -f "$OUT" "$LOG"; }
trap cleanup EXIT

# --- Sanity checks
[[ -s "$VCF_IN" ]] || { echo "FAIL: missing input VCF: $VCF_IN"; exit 1; }
[[ -s "$BAM" ]]    || { echo "FAIL: missing input BAM: $BAM"; exit 1; }

# Ensure BAM index exists
if [[ ! -s "$BAI" ]]; then
  if command -v samtools >/dev/null 2>&1; then
    samtools index "$BAM"
  else
    echo "FAIL: BAM index missing ($BAI) and samtools not available to create it"
    exit 1
  fi
fi

# --- Extract the 2 SNV positions from the input VCF (first two variant lines)
# We assume test_11 VCF contains exactly the 2 SNVs to merge (same codon).
read -r CHROM POS1 POS2 < <(
  zcat -f "$VCF_IN" \
  | awk 'BEGIN{c=0} !/^#/ {c++; if(c==1){chr=$1;p1=$2} else if(c==2){p2=$2; print chr, p1, p2; exit}}'
)

[[ -n "${CHROM:-}" && -n "${POS1:-}" && -n "${POS2:-}" ]] || {
  echo "FAIL: could not read first two SNVs from $VCF_IN"
  exit 1
}

# sorted MERGED_POS
if (( POS1 < POS2 )); then
  MERGED_POS_EXPECTED="${POS1},${POS2}"
else
  MERGED_POS_EXPECTED="${POS2},${POS1}"
fi

# --- Run hgene codon haplotype merge
# Args: <in.vcf(.gz)> <ref.fasta> <ref.gff> <bam> [hap_af_min] [min_inf_reads] [min_mapq] [min_baseq]
python3 "$SRC/hgene_codon_haplotype_merge.py" \
  "$VCF_IN" \
  "$FASTA" \
  "$GFF" \
  "$BAM" \
  0.10 10 40 20 \
  > "$OUT" 2> "$LOG"

# ============== Assertions ==============

# 1) Output must contain at least 1 variant line
out_count="$(awk '!/^#/ {c++} END{print c+0}' "$OUT")"
[[ "$out_count" -ge 1 ]] || {
  echo "FAIL: expected >=1 output variant, got $out_count"
  echo "---- OUT ----"; cat "$OUT"
  echo "---- LOG ----"; cat "$LOG"
  exit 1
}

# 2) Must contain MERGED_POS=pos1,pos2 (the codon merge happened)
grep -Fq "MERGED_POS=${MERGED_POS_EXPECTED}" "$OUT" || {
  echo "FAIL: MERGED_POS=${MERGED_POS_EXPECTED} not found in output"
  echo "---- OUT ----"; cat "$OUT"
  echo "---- LOG ----"; cat "$LOG"
  exit 1
}

# 3) Must contain haplotype info fields (at least HAP + HAP_DP + HAP_CT)
grep -Eq 'HAP=[RAO]{2,3}' "$OUT" || {
  echo "FAIL: HAP field missing (expected HAP=.. in INFO)"
  echo "---- OUT ----"; cat "$OUT"
  exit 1
}
grep -Fq 'HAP_DP=' "$OUT" || { echo "FAIL: HAP_DP missing"; echo "---- OUT ----"; cat "$OUT"; exit 1; }
grep -Fq 'HAP_CT=' "$OUT" || { echo "FAIL: HAP_CT missing"; echo "---- OUT ----"; cat "$OUT"; exit 1; }

# 4) The original SNVs should have been replaced (no remaining single-base records at POS1 or POS2)
# (We check for a record with POS=POS1/2 AND REF/ALT length 1; if merge happened, POS becomes codon anchor and REF/ALT are length 3)
if awk -v p="$POS1" '(!/^#/ && $2==p && length($4)==1 && length($5)==1){exit 0} END{exit 1}' "$OUT"; then
  echo "FAIL: original SNV at POS=$POS1 still present (expected it to be merged away)"
  echo "---- OUT ----"; cat "$OUT"
  exit 1
fi
if awk -v p="$POS2" '(!/^#/ && $2==p && length($4)==1 && length($5)==1){exit 0} END{exit 1}' "$OUT"; then
  echo "FAIL: original SNV at POS=$POS2 still present (expected it to be merged away)"
  echo "---- OUT ----"; cat "$OUT"
  exit 1
fi

# 5) Codon-level annotations should exist
grep -Fq 'CODON_REF=' "$OUT" || { echo "FAIL: CODON_REF missing"; echo "---- OUT ----"; cat "$OUT"; exit 1; }
grep -Fq 'CODON_ALT=' "$OUT" || { echo "FAIL: CODON_ALT missing"; echo "---- OUT ----"; cat "$OUT"; exit 1; }
grep -Fq 'AA_REF=' "$OUT"    || { echo "FAIL: AA_REF missing"; echo "---- OUT ----"; cat "$OUT"; exit 1; }
grep -Fq 'AA_ALT=' "$OUT"    || { echo "FAIL: AA_ALT missing"; echo "---- OUT ----"; cat "$OUT"; exit 1; }
grep -Fq 'DP_CODON=' "$OUT"  || { echo "FAIL: DP_CODON missing"; echo "---- OUT ----"; cat "$OUT"; exit 1; }

echo "OK test_11_merge_two_mut_same_codon"

