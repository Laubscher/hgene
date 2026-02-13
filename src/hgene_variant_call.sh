#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

usage() {
  cat >&2 <<'EOF'
Usage: hgene_variant_call.sh <prefix> <virus> <CPU>

This is a pipeline step. It expects upstream files to exist:
  - <prefix>.sam

It produces:
  - <prefix>_output/ (BAM + VCF + final VCF.GZ + index)
EOF
  exit 2
}

die() { echo "ERROR: $*" >&2; exit 1; }
need_cmd() { command -v "$1" >/dev/null 2>&1 || die "Missing required command: $1"; }

timestamp() { date "+%Y-%m-%d %H:%M:%S"; }
log() { local level="$1"; shift; echo "[$(timestamp)] [$level] $*"; }
info() { log "INFO" "$*"; }
step() { log "STEP" "$*"; }
warn() { log "WARN" "$*"; }

STEP_START=0
start_timer() { STEP_START=$(date +%s); }
end_timer() {
  local end dur
  end=$(date +%s)
  dur=$((end - STEP_START))
  log "TIME" "Step duration: ${dur}s"
}

rmf() { rm -f -- "$@" 2>/dev/null || true; }
index_vcfgz() { bcftools index --csi "$1"; }

# -------------------- args --------------------
[[ $# -eq 3 ]] || usage
prefix="$1"
virus="$2"
CPU="$3"

# -------------------- deps --------------------
need_cmd samtools
need_cmd lofreq
need_cmd bcftools
need_cmd bgzip
need_cmd gunzip
need_cmd python3
need_cmd awk

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]:-$0}")" >/dev/null 2>&1 && pwd)"
DB_DIR="${SCRIPT_DIR}/../db"
FASTA="${DB_DIR}/${virus}.fasta"
GFF="${DB_DIR}/${virus}.gff"

# --- Provenance: reference accession via .info + FASTA SHA tracking ---
# Looks for: <DB_DIR>/<virus>.info  (same stem as <virus>.fasta)
# Recommended keys:
#   ACCESSION=KX035107.1
#   FASTA_SHA256=<expected sha>
# Optional:
#   GFF_SHA256=<expected sha>

REF_INFO="${HG_REF_INFO:-${DB_DIR}/${virus}.info}"

REF_FASTA_SHA256_ACTUAL="$(sha256sum "$FASTA" | awk '{print $1}')"
REF_GFF_SHA256_ACTUAL="$(sha256sum "$GFF" | awk '{print $1}')"

REF_FASTA_SHA256_EXPECTED=""
REF_GFF_SHA256_EXPECTED=""
REF_INFO_STATUS="no_info"
REF_ACCESSION="sha256:${REF_FASTA_SHA256_ACTUAL}"

if [[ -s "$REF_INFO" ]]; then
  # shellcheck disable=SC1090
  source "$REF_INFO" || true
  REF_INFO_STATUS="loaded"
  REF_FASTA_SHA256_EXPECTED="${FASTA_SHA256:-}"
  REF_GFF_SHA256_EXPECTED="${GFF_SHA256:-}"

  if [[ -n "${REF_FASTA_SHA256_EXPECTED}" && "$REF_FASTA_SHA256_ACTUAL" == "$REF_FASTA_SHA256_EXPECTED" ]]; then
    REF_INFO_STATUS="ok"
    if [[ -n "${ACCESSION:-}" ]]; then
      REF_ACCESSION="$ACCESSION"
    else
      # info exists and sha matches, but accession missing -> keep sha identifier
      REF_ACCESSION="sha256:${REF_FASTA_SHA256_ACTUAL}"
    fi
  else
    REF_INFO_STATUS="fasta_sha_mismatch"
    REF_ACCESSION="sha256:${REF_FASTA_SHA256_ACTUAL}"
  fi
fi

info "HG_VERSION=${HG_VERSION:-unknown}"
info "REF_INFO=$REF_INFO"
info "REF_INFO_STATUS=$REF_INFO_STATUS"
info "REF_ACCESSION=$REF_ACCESSION"
info "REF_FASTA_SHA256=$REF_FASTA_SHA256_ACTUAL"
[[ -n "$REF_FASTA_SHA256_EXPECTED" ]] && info "REF_FASTA_SHA256_EXPECTED=$REF_FASTA_SHA256_EXPECTED"
info "REF_GFF_SHA256=$REF_GFF_SHA256_ACTUAL"
[[ -n "$REF_GFF_SHA256_EXPECTED" ]] && info "REF_GFF_SHA256_EXPECTED=$REF_GFF_SHA256_EXPECTED"
# ------------------------------------------------------------------


[[ -s "${prefix}.sam" ]] || die "Input SAM not found: ${prefix}.sam"
[[ -s "$FASTA" ]] || die "Reference FASTA not found: $FASTA"
[[ -s "$GFF"   ]] || die "Reference GFF not found:   $GFF"

export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}:/usr/local/lib"

# -------------------- workdir --------------------
TMPDIR="$(mktemp -d "${prefix}.tmp.XXXXXX")"
trap 'rm -rf "$TMPDIR" 2>/dev/null || true' EXIT

# Optional: remove previous adapter-trim fastq if present
rmf "${prefix}_tr.fastq"

# -------------------- BAM preprocessing --------------------
step "samtools sort/index (raw)"
start_timer
samtools sort "${prefix}.sam" -o "${TMPDIR}/${prefix}.sorted1.bam"
samtools index "${TMPDIR}/${prefix}.sorted1.bam"
end_timer

step "samtools view MAPQ>=40"
start_timer
samtools view -bq 40 "${TMPDIR}/${prefix}.sorted1.bam" > "${TMPDIR}/${prefix}.filtered2.bam"
end_timer

step "samtools view length filter"
start_timer

if [[ "$virus" == "HHV5" ]]; then
  samtools view -e 'rlen>99'  -O BAM -o "${TMPDIR}/${prefix}.filtered.bam" "${TMPDIR}/${prefix}.filtered2.bam"
else
  samtools view -e 'rlen>999' -O BAM -o "${TMPDIR}/${prefix}.filtered.bam" "${TMPDIR}/${prefix}.filtered2.bam"
fi
end_timer

step "samtools sort/index (filtered)"
start_timer
samtools sort "${TMPDIR}/${prefix}.filtered.bam" -o "${TMPDIR}/${prefix}.sorted.bam"
samtools index "${TMPDIR}/${prefix}.sorted.bam"
end_timer

# remove intermediate BAMs produced inside TMPDIR
rmf "${TMPDIR}/${prefix}.sorted1.bam" "${TMPDIR}/${prefix}.sorted1.bam.bai"
rmf "${TMPDIR}/${prefix}.filtered2.bam" "${TMPDIR}/${prefix}.filtered.bam"

# -------------------- LoFreq indelqual -> final BAM --------------------
step "LoFreq indelqual -> ${prefix}.bam"
start_timer
lofreq indelqual -u 16 -f "$FASTA" "${TMPDIR}/${prefix}.sorted.bam" > "${prefix}.bam"
rmf "${TMPDIR}/${prefix}.sorted.bam" "${TMPDIR}/${prefix}.sorted.bam.bai"
samtools index "${prefix}.bam"
end_timer

# -------------------- LoFreq calling --------------------
step "LoFreq call-parallel -> ${prefix}.lofreq_raw.vcf.gz"
start_timer
lofreq call-parallel --pp-threads "$CPU" \
  -f "$FASTA" \
  --no-default-filter -A -B -a 1 -b 1 \
  --force-overwrite --call-indels \
  "${prefix}.bam" \
  -o "${prefix}.lofreq_raw.vcf.gz"
end_timer

step "gunzip LoFreq VCF -> ${prefix}.lofreq_raw.vcf"
start_timer
gunzip -f "${prefix}.lofreq_raw.vcf.gz"
rmf "${prefix}.lofreq_raw.vcf.gz.tbi"
end_timer

# -------------------- LoFreq filter --------------------
step "LoFreq filter (cov>=20, af>=0.1) -> ${prefix}.lofreq_filtered.vcf"
start_timer
lofreq filter --no-defaults --cov-min 20 --af-min 0.1 \
  -i "${prefix}.lofreq_raw.vcf" -o "${prefix}.lofreq_filtered.vcf"
end_timer

# -------------------- Custom filtering hgene_filter_custom.py --------------------
step "hgene_filter_custom.py -> ${prefix}.filtered.vcf"
start_timer
python3 "${SCRIPT_DIR}/hgene_filter_custom.py" "${prefix}.lofreq_filtered.vcf" > "${prefix}.filtered.vcf"
end_timer

step "bgzip + bcftools index ${prefix}.filtered.vcf"
start_timer
bgzip -f "${prefix}.filtered.vcf"
index_vcfgz "${prefix}.filtered.vcf.gz"
end_timer

# -------------------- Normalize + csq --------------------
step "bcftools norm -m- | bcftools csq -> ${prefix}.bcf.vcf.gz"
start_timer
bcftools norm -Ou -m- -f "$FASTA" "${prefix}.filtered.vcf.gz" \
| bcftools csq -Ob --local-csq --force -o "${prefix}.bcf.vcf.gz" \
    --fasta-ref="$FASTA" \
    --gff="$GFF"
end_timer

# -------------------- Codon haplotypes --------------------
step "codon haplotypes (hgene_codon_haplotype_merge.py) -> ${prefix}.vcf"
start_timer
python3 "${SCRIPT_DIR}/hgene_codon_haplotype_merge.py" \
  "${prefix}.bcf.vcf.gz" \
  "$FASTA" \
  "$GFF" \
  "${prefix}.bam" \
  0.1 \
  10 \
  40 \
  20 \
  > "${prefix}.vcf"
end_timer

# --- Add provenance lines to final VCF header (plain VCF before bgzip) ---
if command -v bcftools >/dev/null 2>&1; then
  tmp_vcf="${prefix}.vcf.hgene.tmp"
  bcftools annotate -h <(cat <<EOF
##hgene_version=${HG_VERSION:-unknown}
##hgene_ref_info=${REF_INFO}
##hgene_ref_info_status=${REF_INFO_STATUS}
##hgene_ref_accession=${REF_ACCESSION}
##hgene_ref_fasta_sha256=${REF_FASTA_SHA256_ACTUAL}
##hgene_ref_gff_sha256=${REF_GFF_SHA256_ACTUAL}
EOF
  if [[ -n "${REF_FASTA_SHA256_EXPECTED:-}" ]]; then echo "##hgene_ref_fasta_sha256_expected=${REF_FASTA_SHA256_EXPECTED}"; fi
  if [[ -n "${REF_GFF_SHA256_EXPECTED:-}" ]]; then echo "##hgene_ref_gff_sha256_expected=${REF_GFF_SHA256_EXPECTED}"; fi
  cat <<EOF2
##hgene_ref_fasta=${FASTA}
##hgene_ref_gff=${GFF}
EOF2
  ) -Ov -o "$tmp_vcf" "${prefix}.vcf"
  mv -f "$tmp_vcf" "${prefix}.vcf"
else
  warn "bcftools not found: skipping hgene provenance header injection into ${prefix}.vcf"
fi
# -------------------------------------------------------------------


step "bgzip + bcftools index final -> ${prefix}.vcf.gz"
start_timer
bgzip -f "${prefix}.vcf"
index_vcfgz "${prefix}.vcf.gz"
end_timer

# -------------------- Packaging --------------------
step "packaging outputs -> ${prefix}_output/"
start_timer

mkdir -p "BAM_${prefix}" "VCF_${prefix}" "${prefix}_output"

mv -f "${prefix}.bam" "BAM_${prefix}/"
mv -f "${prefix}.bam.bai" "BAM_${prefix}/"

mv -f "${prefix}.lofreq_filtered.vcf" "VCF_${prefix}/"
mv -f "${prefix}.lofreq_raw.vcf" "VCF_${prefix}/"
mv -f "${prefix}.filtered.vcf.gz" "VCF_${prefix}/"
mv -f "${prefix}.filtered.vcf.gz.csi" "VCF_${prefix}/"
mv -f "${prefix}.bcf.vcf.gz" "VCF_${prefix}/"

mv -f "VCF_${prefix}" "${prefix}_output/"
mv -f "BAM_${prefix}" "${prefix}_output/"

mv -f "${prefix}.vcf.gz" "${prefix}_output/"
mv -f "${prefix}.vcf.gz.csi" "${prefix}_output/"

end_timer