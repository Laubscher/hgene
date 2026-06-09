#!/usr/bin/env bash
# Author: Florian Laubscher
# Part of: hgene pipeline

set -euo pipefail
IFS=$'\n\t'

get_hg_version() {
  local here verfile
  here="$(cd -- "$(dirname -- "${BASH_SOURCE[0]:-$0}")" &>/dev/null && pwd)"

  verfile="$here/VERSION"
  if [[ -s "$verfile" ]]; then
    tr -d '\r\n' < "$verfile"
    return 0
  fi

  # fallback: git describe only if available (dev mode)
  if command -v git >/dev/null 2>&1; then
    if git -C "$here" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
      git -C "$here" describe --tags --always --dirty 2>/dev/null && return 0
    fi
  fi

  echo "unknown"
}

HG_VERSION="$(get_hg_version)"
export HG_VERSION

usage() {
  cat >&2 <<'USAGE'
Usage:
  hgene_main.sh <prefix> <CPU> <virus> [resistance_db_dir] [mode] [input1] [input2]

Modes:
  ont                  ONT mode. Expects <prefix>.fastq.
  illumina-comparison  Paired-end Illumina comparison mode. Expects input1=R1 and input2=R2.
USAGE
  exit 2
}

timestamp() { date "+%Y-%m-%d %H:%M:%S"; }
log() { local level="$1"; shift; echo "[$(timestamp)] [$level] $*"; }
info() { log "INFO" "$*"; }
step() { log "STEP" "$*"; }
warn() { log "WARN" "$*"; }
error() { log "ERROR" "$*"; }
die() { error "$*"; exit 1; }
need_cmd() { command -v "$1" >/dev/null 2>&1 || die "Missing required command: $1"; }

# Args:
#   $1 prefix
#   $2 CPU
#   $3 virus
#   $4 resistance DB directory
#   $5 mode
#   $6 input1
#   $7 input2

[[ $# -ge 3 && $# -le 7 ]] || usage

prefix="$1"
CPU="$2"
virus="$3"
HG_RESISTANCE_DB_DIR="${4:-}"
mode="${5:-ont}"
input1="${6:-}"
input2="${7:-}"

[[ "$mode" == "ont" || "$mode" == "illumina-comparison" ]] || die "Invalid mode: $mode. Expected: ont or illumina-comparison"

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]:-$0}")" >/dev/null 2>&1 && pwd)"

info "HG_VERSION=$HG_VERSION"
info "hgene_main started"
info "prefix=$prefix virus=$virus cpu=$CPU mode=$mode"

cleanup() {
  rm -f "${prefix}.sam" "${prefix}.trimmed.fastq" 2>/dev/null || true
}
trap cleanup EXIT

need_cmd bash

# -------------------- preprocessing + mapping --------------------

if [[ "$mode" == "ont" ]]; then
  need_cmd porechop

  ONT_FASTQ="${prefix}.fastq"
  [[ -s "$ONT_FASTQ" ]] || die "Input FASTQ not found: $ONT_FASTQ"

  step "Adapter trimming (porechop)"
  bash "$SCRIPT_DIR/hgene_split_and_porechop.sh" "$ONT_FASTQ" "${prefix}.trimmed.fastq" "$CPU"

  step "Mapping to reference (hgene_map.sh)"
  bash "${SCRIPT_DIR}/hgene_map.sh" "${prefix}.trimmed.fastq" "$virus" "$CPU" "$prefix"

elif [[ "$mode" == "illumina-comparison" ]]; then
  [[ -n "$input1" && -n "$input2" ]] || die "Illumina-comparison mode expects R1 and R2 FASTQ files"
  [[ -s "$input1" ]] || die "R1 FASTQ not found: $input1"
  [[ -s "$input2" ]] || die "R2 FASTQ not found: $input2"

  warn "Illumina-comparison mode: skipping porechop"
  warn "Illumina-comparison mode is intended for technical comparison only and is not validated"

  step "Mapping paired-end Illumina reads to reference (hgene_map_illumina.sh)"
  bash "${SCRIPT_DIR}/hgene_map_illumina.sh" "$input1" "$input2" "$virus" "$CPU" "$prefix"
fi

# -------------------- variant calling --------------------

step "Variant calling (hgene_variant_call.sh)"
bash "${SCRIPT_DIR}/hgene_variant_call.sh" "$prefix" "$virus" "$CPU" "$mode"


# --- Auto user template based on virus ---

HG_TEMPLATE_DOCX=""

TEMPLATE_ROOT="${HG_TEMPLATE_ROOT:-$HOME/template}"
case "$virus" in
  HHV1) CANDIDATE="${TEMPLATE_ROOT}/hsv1/template.docx" ;;
  HHV2) CANDIDATE="${TEMPLATE_ROOT}/hsv2/template.docx" ;;
  HHV5) CANDIDATE="${TEMPLATE_ROOT}/cmv/template.docx" ;;
  *) CANDIDATE="" ;;
esac

if [[ -n "$CANDIDATE" && -s "$CANDIDATE" ]]; then
  HG_TEMPLATE_DOCX="$CANDIDATE"
  info "User template detected: $HG_TEMPLATE_DOCX"
else
  if [[ -n "$CANDIDATE" ]]; then
    warn "User template not found: $CANDIDATE"
    warn "Falling back to default DB template"
  fi
fi

if [[ "${virus}" == "HHV1" || "${virus}" == "HHV2" || "${virus}" == "HHV5" ]]; then
  step "step Virotyper report (hgene_virotype_report.sh)"
  bash "${SCRIPT_DIR}/hgene_virotype_report.sh" "${virus}" "${prefix}" "${HG_TEMPLATE_DOCX:-}" "${HG_RESISTANCE_DB_DIR:-}"
fi

info "hgene_main finished"