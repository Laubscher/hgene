#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

usage() {
  cat >&2 <<'EOF'
Usage:
  hgene_virotype_report.sh <virus> <prefix>

EOF
  exit 2
}

timestamp() { date "+%Y-%m-%d %H:%M:%S"; }
log() { local level="$1"; shift; echo "[$(timestamp)] [$level] $*"; }
info() { log "INFO" "$*"; }
step() { log "STEP" "$*"; }
error() { log "ERROR" "$*"; }

virus="${1:-}"; prefix="${2:-}"
[[ -n "$virus" && -n "$prefix" ]] || usage

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]:-$0}")" >/dev/null 2>&1 && pwd)"

VT_ROOT="$(cd -- "${SCRIPT_DIR}/../" && pwd)"

# Map virus key -> virotyper DB_ID
case "$virus" in
  HHV1) DB_ID="hsv1" ;;
  HHV2) DB_ID="hsv2" ;;
  *) echo "ERROR: virotyper report enabled only for HHV1/HHV2 (got: $virus)" >&2; exit 1 ;;
esac

OUTDIR="${prefix}_output"
REPORT_DIR="${OUTDIR}/report"
mkdir -p "$REPORT_DIR"

VCF_SRC="${OUTDIR}/${prefix}.vcf.gz"
VCF_CSI="${VCF_SRC}.csi"
BAM_SRC="${OUTDIR}/BAM_${prefix}/${prefix}.bam"
BAM_BAI="${BAM_SRC}.bai"

[[ -s "$VCF_SRC" ]] || { echo "ERROR: missing VCF: $VCF_SRC" >&2; exit 1; }
[[ -s "$VCF_CSI" ]] || { echo "ERROR: missing VCF index: $VCF_CSI" >&2; exit 1; }
[[ -s "$BAM_SRC" ]] || { echo "ERROR: missing BAM: $BAM_SRC" >&2; exit 1; }
[[ -s "$BAM_BAI" ]] || { echo "ERROR: missing BAM index: $BAM_BAI" >&2; exit 1; }

[[ -d "$VT_ROOT" ]] || { echo "ERROR: virotyper directory not found: $VT_ROOT" >&2; exit 1; }
[[ -f "${VT_ROOT}/notebooks/bam_report.Rmd" ]] || { echo "ERROR: missing Rmd: ${VT_ROOT}/notebooks/bam_report.Rmd" >&2; exit 1; }
[[ -f "${VT_ROOT}/notebooks/vcf_report.Rmd" ]] || { echo "ERROR: missing Rmd: ${VT_ROOT}/notebooks/vcf_report.Rmd" >&2; exit 1; }
[[ -d "${VT_ROOT}/data/db/${DB_ID}" ]] || { echo "ERROR: missing DB dir: ${VT_ROOT}/data/db/${DB_ID}" >&2; exit 1; }

pushd "$VT_ROOT" >/dev/null

# Copy inputs locally
cp -f "$VCF_SRC" .
cp -f "$VCF_CSI" .
cp -f "$BAM_SRC" .
cp -f "$BAM_BAI" .

VCF_LOCAL="${VT_ROOT}/$(basename "$VCF_SRC")"
BAM_LOCAL="${VT_ROOT}/$(basename "$BAM_SRC")"

# Render BAM report
Rscript -e "rmarkdown::render(
  'notebooks/bam_report.Rmd',
  params=list(input_bam_file='${BAM_LOCAL}'),
  knit_root_dir='${VT_ROOT}',
  output_dir='${REPORT_DIR}',
  output_file='${prefix}_bam'
)"

# Render VCF report (+ docx path via params)
Rscript -e "rmarkdown::render(
  'notebooks/vcf_report.Rmd',
  params=list(
    input_vcf_file='${VCF_LOCAL}',
    input_db_dir='${VT_ROOT}/data/db/${DB_ID}',
    input_fasta='${VT_ROOT}/db/${virus}.fasta',
    output_docx_report='${prefix}.vcf.gz.${DB_ID}.docx'
  ),
  knit_root_dir='${VT_ROOT}',
  output_dir='${REPORT_DIR}',
  output_file='${prefix}'
)"


if [[ -f "${VT_ROOT}/${prefix}.vcf.gz.${DB_ID}.docx" ]]; then
  mv -f "${VT_ROOT}/${prefix}.vcf.gz.${DB_ID}.docx" "${REPORT_DIR}/"
fi

# Cleanup local copies
rm -f "$(basename "$VCF_SRC")" "$(basename "$VCF_CSI")" "$(basename "$BAM_SRC")" "$(basename "$BAM_BAI")"

popd >/dev/null

info "Virotyper reports written to ${REPORT_DIR}/"

