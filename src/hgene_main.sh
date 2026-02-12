#!/usr/bin/env bash
set -euo pipefail
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