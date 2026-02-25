#!/usr/bin/env bash
# Author: Florian Laubscher
# Part of: hgene pipeline

set -euo pipefail


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

if declare -F info >/dev/null 2>&1; then
  info "HG_VERSION=$HG_VERSION"
fi
# ------------------------------------------------------------------

IFS=$'\n\t'

# Args:
#   $1 prefix
#   $2 CPU
#   $3 virus

usage() {
  cat >&2 <<'USAGE'
Usage: hgene_main.sh <prefix> <CPU> <virus>

Expects:
  <prefix>.fastq
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

cleanup() {
    rm -f "${prefix}.sam" "${prefix}.trimmed.fastq" 2>/dev/null || true
}
trap cleanup EXIT

[[ $# -eq 3 ]] || usage
prefix="$1"
CPU="$2"
virus="$3"

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]:-$0}")" >/dev/null 2>&1 && pwd)"

need_cmd porechop
need_cmd bash

[[ -s "${prefix}.fastq" ]] || die "Input FASTQ not found: ${prefix}.fastq"

step "Adapter trimming (porechop)"
# keep porechop quiet; errors still surface
porechop -t "$CPU" --discard_middle -i "${prefix}.fastq" -o "${prefix}.trimmed.fastq" >/dev/null # trimming

step "Mapping to reference (hgene_map.sh)"
bash "${SCRIPT_DIR}/hgene_map.sh" "${prefix}.trimmed.fastq" "$virus" "$CPU" "$prefix"

step "Variant calling (hgene_variant_call.sh)"
bash "${SCRIPT_DIR}/hgene_variant_call.sh" "$prefix" "$virus" "$CPU"


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
  bash "${SCRIPT_DIR}/hgene_virotype_report.sh" "${virus}" "${prefix}" "${HG_TEMPLATE_DOCX:-}"
fi
