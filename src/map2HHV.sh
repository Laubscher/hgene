#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# Args:
#   $1 prefix
#   $2 virus
#   $3 CPU (minimap2 threads)

usage() {
  cat >&2 <<'USAGE'
Usage: map2HHV.sh <prefix> <virus> <CPU>

Expects:
  <prefix>_tr.fastq
Produces:
  <prefix>.sam
USAGE
  exit 2
}

timestamp() { date "+%Y-%m-%d %H:%M:%S"; }
log() { local level="$1"; shift; echo "[$(timestamp)] [$level] $*"; }
step() { log "STEP" "$*"; }
error() { log "ERROR" "$*"; }

die() { error "$*"; exit 1; }
need_cmd() { command -v "$1" >/dev/null 2>&1 || die "Missing required command: $1"; }

[[ $# -eq 3 ]] || usage
prefix="$1"
virus="$2"
CPU="$3"

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]:-$0}")" >/dev/null 2>&1 && pwd)"
REF="${SCRIPT_DIR}/../db/${virus}.fasta"

need_cmd minimap2

[[ -s "${prefix}_tr.fastq" ]] || die "Input FASTQ not found: ${prefix}_tr.fastq"
[[ -s "$REF" ]] || die "Reference FASTA not found: $REF"

step "minimap2 -ax map-ont (threads=$CPU)"
minimap2 -t "$CPU" -ax map-ont "$REF" "${prefix}_tr.fastq" > "${prefix}.sam"