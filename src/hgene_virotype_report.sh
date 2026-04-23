#!/usr/bin/env bash
# Author: Florian Laubscher
# Part of: hgene pipeline

set -euo pipefail
IFS=$'\n\t'

usage() {
  cat >&2 <<'EOF'
Usage:
  hgene_virotype_report.sh <virus> <prefix> [template.docx] [resistance_db_root]

EOF
  exit 2
}

timestamp() { date "+%Y-%m-%d %H:%M:%S"; }
log() { local level="$1"; shift; echo "[$(timestamp)] [$level] $*"; }
info() { log "INFO" "$*"; }
step() { log "STEP" "$*"; }
error() { log "ERROR" "$*"; }

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]:-$0}")" >/dev/null 2>&1 && pwd)"

VT_ROOT="$(cd -- "${SCRIPT_DIR}/../" && pwd)"

virus="${1:-}"; prefix="${2:-}"
[[ -n "$virus" && -n "$prefix" ]] || usage

# Map virus key -> virotyper DB_ID
case "$virus" in
  HHV1) DB_ID="hsv1" ;;
  HHV2) DB_ID="hsv2" ;;
  HHV5) DB_ID="cmv" ;;
  *) echo "ERROR: virotyper report enabled only for HHV1/HHV2/HHV5 (got: $virus)" >&2; exit 1 ;;
esac

USER_TEMPLATE="${3:-}"
RESISTANCE_DB_ROOT="${4:-}"

DEFAULT_DB_DIR="${VT_ROOT}/data/db/${DB_ID}"
DB_DIR="${DEFAULT_DB_DIR}"

if [[ -n "${RESISTANCE_DB_ROOT:-}" ]]; then
  case "$virus" in
    HHV1)
      for candidate in \
        "${RESISTANCE_DB_ROOT}/HHV1" \
        "${RESISTANCE_DB_ROOT}/hsv1"
      do
        if [[ -s "${candidate}/resistances.chk.xlsx" ]]; then
          DB_DIR="${candidate}"
          info "Using external resistance DB directory: ${DB_DIR}"
          break
        fi
      done
      ;;
    HHV2)
      for candidate in \
        "${RESISTANCE_DB_ROOT}/HHV2" \
        "${RESISTANCE_DB_ROOT}/hsv2"
      do
        if [[ -s "${candidate}/resistances.chk.xlsx" ]]; then
          DB_DIR="${candidate}"
          info "Using external resistance DB directory: ${DB_DIR}"
          break
        fi
      done
      ;;
    HHV5)
      for candidate in \
        "${RESISTANCE_DB_ROOT}/HHV5" \
        "${RESISTANCE_DB_ROOT}/cmv"
      do
        if [[ -s "${candidate}/resistances.chk.xlsx" ]]; then
          DB_DIR="${candidate}"
          info "Using external resistance DB directory: ${DB_DIR}"
          break
        fi
      done
      ;;
  esac
fi

if [[ ! -s "${DB_DIR}/resistances.chk.xlsx" ]]; then
  echo "ERROR: missing resistance DB file: ${DB_DIR}/resistances.chk.xlsx" >&2
  exit 1
fi

# template par défaut
DB_TEMPLATE="${DEFAULT_DB_DIR}/template.docx"


if [[ -n "${USER_TEMPLATE:-}" ]]; then
  if [[ -s "$USER_TEMPLATE" ]]; then
    TEMPLATE_DOCX="$USER_TEMPLATE"
    info "Using user template: $TEMPLATE_DOCX"
  else
    warn "User template not found: $USER_TEMPLATE"
    warn "Falling back to default template: $DB_TEMPLATE"
    TEMPLATE_DOCX=${DB_TEMPLATE}
  fi
else
  info "Using default template: $DB_TEMPLATE"
  TEMPLATE_DOCX=${DB_TEMPLATE}
fi


RUN_DIR="$(pwd)"
OUTDIR="${RUN_DIR}/${prefix}_output"
REPORT_DIR="${OUTDIR}/REPORT_${prefix}"
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


# Render BAM report
cp -f "${VT_ROOT}/notebooks/bam_report.Rmd" "${REPORT_DIR}/bam_report.Rmd"

Rscript -e "rmarkdown::render(
  '${REPORT_DIR}/bam_report.Rmd',
  params=list(input_bam_file='${BAM_SRC}',output_coverage_tsv='${REPORT_DIR}/coverage_summary.tsv'),
  knit_root_dir='${REPORT_DIR}',
  output_dir='${REPORT_DIR}',
  output_file='${prefix}_bam'
)"

# Render VCF report (+ docx path via params)
cp -f "${VT_ROOT}/notebooks/vcf_report.Rmd" "${REPORT_DIR}/vcf_report.Rmd"
Rscript -e "rmarkdown::render(
  '${REPORT_DIR}/vcf_report.Rmd',
  params=list(
    input_vcf_file='${VCF_SRC}',
    input_db_dir='${DB_DIR}',
    input_fasta='${VT_ROOT}/db/${virus}.fasta',
    template_docx='${TEMPLATE_DOCX}',
    output_docx_report='${prefix}.vcf.gz.${DB_ID}.docx',
    input_coverage_tsv='${REPORT_DIR}/coverage_summary.tsv'
  ),
  knit_root_dir='${REPORT_DIR}',
  output_dir='${REPORT_DIR}',
  output_file='${prefix}'
)"

# Cleanup local copies
rm -f "${REPORT_DIR}/bam_report.Rmd"
rm -f "${REPORT_DIR}/vcf_report.Rmd"

if [[ -f "${VT_ROOT}/${prefix}.vcf.gz.${DB_ID}.docx" ]]; then
  mv -f "${VT_ROOT}/${prefix}.vcf.gz.${DB_ID}.docx" "${REPORT_DIR}/"
fi

info "Virotyper reports written to ${REPORT_DIR}/"
