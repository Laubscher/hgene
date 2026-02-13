#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat >&2 <<'EOF'
Usage: make_ref_info.sh <ref.fasta> <ref.gff> <accession> [out.info]

Writes a <ref>.info file with:
  ACCESSION=...
  FASTA_SHA256=...
  GFF_SHA256=...

If out.info is not provided:
  <ref.fasta without .fasta>.info
EOF
  exit 2
}

[[ $# -ge 3 ]] || usage

FASTA="$1"
GFF="$2"
ACC="$3"
OUT="${4:-}"

[[ -s "$FASTA" ]] || { echo "ERROR: FASTA not found: $FASTA" >&2; exit 1; }
[[ -s "$GFF"   ]] || { echo "ERROR: GFF not found: $GFF" >&2; exit 1; }

if [[ -z "$OUT" ]]; then
  if [[ "$FASTA" == *.fasta ]]; then
    OUT="${FASTA%.fasta}.info"
  else
    OUT="${FASTA}.info"
  fi
fi

FASTA_SHA="$(sha256sum "$FASTA" | awk '{print $1}')"
GFF_SHA="$(sha256sum "$GFF" | awk '{print $1}')"

cat > "$OUT" <<EOF
ACCESSION=$ACC
FASTA_SHA256=$FASTA_SHA
GFF_SHA256=$GFF_SHA
EOF

echo "Wrote: $OUT"
