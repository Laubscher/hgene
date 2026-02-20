#!/usr/bin/env bash
# Author: Florian Laubscher
# Part of: hgene pipeline

set -euo pipefail
IFS=$'\n\t'

# Args:
#   $1 input_fastq
#   $2 virus
#   $3 CPU (minimap2 threads)
#   $4 prefix

usage() {
  cat >&2 <<'USAGE'
Usage: hgene_map.sh <input_fastq> <virus> <CPU> [prefix]

Expects:
  <prefix>.trimmed.fastq
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


[[ $# -ge 3 && $# -le 4 ]] || usage
input_fastq="$1"
virus="$2"
CPU="$3"

if [[ $# -ge 4 ]]; then
    prefix="$4"
else
    # derive from filename (cut at first dot)
    base=$(basename "$input_fastq")
    prefix="${base%%.*}"
fi


SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]:-$0}")" >/dev/null 2>&1 && pwd)"
REF="${SCRIPT_DIR}/../db/${virus}.fasta"

need_cmd minimap2

[[ -s "${input_fastq}" ]] || die "Input FASTQ not found: ${input_fastq}"
[[ -s "$REF" ]] || die "Reference FASTA not found: $REF"

step "minimap2 -ax map-ont (threads=$CPU)"
minimap2 -t "$CPU" -ax map-ont -secondary=no -O 6,6 -E 2,2 "$REF" "${input_fastq}" > "${prefix}.sam"