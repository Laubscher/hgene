#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]:-$0}" )" &> /dev/null && pwd )"
SRC="$SCRIPT_DIR/../src"
DATA="$SCRIPT_DIR/data"
DB="$SCRIPT_DIR/db"

FASTA="$DB/HHV1.fasta"
GFF="$DB/HHV1.gff"

VCF_IN="$DATA/test_15_hsv1_XX.bcf.vcf"
BAM="$DATA/test_15_hsv1_XX.bam"
BAI="$DATA/test_15_hsv1_XX.bam.bai"

OUT="test_15.out.vcf"
LOG="test_15.err.log"

cleanup() { rm -f "$OUT" "$LOG"; }
trap cleanup EXIT

[[ -s "$FASTA" ]]  || { echo "FAIL: missing FASTA: $FASTA"; exit 1; }
[[ -s "$GFF"   ]]  || { echo "FAIL: missing GFF: $GFF"; exit 1; }
[[ -s "$VCF_IN" ]] || { echo "FAIL: missing VCF: $VCF_IN"; exit 1; }
[[ -s "$BAM"   ]]  || { echo "FAIL: missing BAM: $BAM"; exit 1; }

if [[ ! -s "$BAI" ]]; then
  command -v samtools >/dev/null 2>&1 || { echo "FAIL: missing BAI and samtools not available"; exit 1; }
  samtools index "$BAM"
fi

python3 "$SRC/hgene_codon_haplotype_merge.py"   "$VCF_IN" "$FASTA" "$GFF" "$BAM"   0.10 10 40 20   > "$OUT" 2> "$LOG"

get_info() {
  awk -v info="$1" -v key="$2" 'BEGIN{
    n=split(info,a,";");
    for(i=1;i<=n;i++){
      if(a[i] ~ ("^"key"=")){
        sub("^"key"=","",a[i]); print a[i]; exit
      }
    }
  }'
}

declare -A CT
declare -A DP
CODON_REF=""

while IFS=$'\t' read -r CHR POS ID REF ALT QUAL FILT INFO; do
  ca="$(get_info "$INFO" "CODON_ALT")"
  cr="$(get_info "$INFO" "CODON_REF")"
  hc="$(get_info "$INFO" "HAP_CT")"
  hd="$(get_info "$INFO" "HAP_DP")"

  [[ -n "$cr" && -z "$CODON_REF" ]] && CODON_REF="$cr"
  [[ -n "$ca" && -n "$hc" ]] && CT["$ca"]="$hc"
  [[ -n "$ca" && -n "$hd" ]] && DP["$ca"]="$hd"
done < <(awk '!/^#/ {print}' "$OUT")

[[ -n "$CODON_REF" ]] || { echo "FAIL: missing CODON_REF in output"; cat "$OUT"; exit 1; }

if [[ -n "${CT[$CODON_REF]:-}" ]]; then
  echo "FAIL: reference codon was emitted"
  cat "$OUT"
  exit 1
fi

n_lines="$(awk '!/^#/{c++} END{print c+0}' "$OUT")"
[[ "$n_lines" -eq 3 ]] || { echo "FAIL: expected 3 output variants, got $n_lines"; cat "$OUT"; exit 1; }

[[ "${CT[GAG]:-}" == "20" ]] || { echo "FAIL: expected GAG HAP_CT=20"; cat "$OUT"; exit 1; }
[[ "${CT[GCG]:-}" == "30" ]] || { echo "FAIL: expected GCG HAP_CT=30"; cat "$OUT"; exit 1; }
[[ "${CT[TAG]:-}" == "10" ]] || { echo "FAIL: expected TAG HAP_CT=10"; cat "$OUT"; exit 1; }

echo "OK test_15_mixed_linked_unlinked_populations"
