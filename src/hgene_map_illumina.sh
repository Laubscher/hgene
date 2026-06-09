#!/usr/bin/env bash
# Author: Florian Laubscher
# Part of: hgene pipeline

set -euo pipefail
IFS=$'\n\t'

timestamp() { date "+%Y-%m-%d %H:%M:%S"; }
log() { local level="$1"; shift; echo "[$(timestamp)] [$level] $*"; }
info() { log "INFO" "$*"; }
step() { log "STEP" "$*"; }
warn() { log "WARN" "$*"; }
error() { log "ERROR" "$*"; }
die() { error "$*"; exit 1; }
need_cmd() { command -v "$1" >/dev/null 2>&1 || die "Missing required command: $1"; }

usage() {
  cat >&2 <<'USAGE'
Usage:
  hgene_map_illumina.sh <R1.fastq[.gz]> <R2.fastq[.gz]> <virus> <CPU> <prefix>

Arguments:
  <R1>      Illumina read 1 FASTQ file
  <R2>      Illumina read 2 FASTQ file
  <virus>   Virus reference key, e.g. HHV1, HHV2, HHV5
  <CPU>     Number of threads
  <prefix>  Output prefix

Produces:
  <prefix>.sam
USAGE
  exit 2
}

[[ $# -eq 5 ]] || usage

R1="$1"
R2="$2"
virus="$3"
CPU="$4"
prefix="$5"

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]:-$0}")" >/dev/null 2>&1 && pwd)"
DB_DIR="${SCRIPT_DIR}/../db"
FASTA="${DB_DIR}/${virus}.fasta"

[[ -s "$R1" ]] || die "R1 FASTQ not found: $R1"
[[ -s "$R2" ]] || die "R2 FASTQ not found: $R2"
[[ -s "$FASTA" ]] || die "Virus reference not found: $FASTA"

need_cmd minimap2

info "Illumina mapping started"
info "R1=$R1"
info "R2=$R2"
info "virus=$virus"
info "reference=$FASTA"
info "prefix=$prefix"
info "cpu=$CPU"

step "Mapping Illumina paired-end reads with minimap2 -ax sr -> ${prefix}.sam"

minimap2 -t "$CPU" -ax sr "$FASTA" "$R1" "$R2" > "${prefix}.sam"

[[ -s "${prefix}.sam" ]] || die "SAM output is empty or missing: ${prefix}.sam"

info "Illumina mapping finished"
info "SAM: ${prefix}.sam"