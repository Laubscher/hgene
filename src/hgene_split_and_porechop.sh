#!/usr/bin/env bash
# Author: Florian Laubscher
# Part of: hgene pipeline

set -euo pipefail

# Parallel Porechop wrapper for hgene
# Usage:
#   hgene_split_and_porechop.sh <input.fastq> <output.fastq> <cpu>

input_fastq="${1:?ERROR: input FASTQ required}"
output_fastq="${2:?ERROR: output FASTQ required}"
cpu="${3:-1}"

if [ ! -f "$input_fastq" ]; then
    echo "ERROR: input FASTQ not found: $input_fastq" >&2
    exit 1
fi

if ! command -v porechop >/dev/null 2>&1; then
    echo "ERROR: porechop not found in PATH" >&2
    exit 1
fi

if ! [[ "$cpu" =~ ^[0-9]+$ ]]; then
    echo "ERROR: CPU must be an integer, got: $cpu" >&2
    exit 1
fi

if [ "$cpu" -lt 1 ]; then
    cpu=1
fi

tmpdir="${TMPDIR:-.}/hgene_porechop_parallel_$$"
mkdir -p "$tmpdir/chunks" "$tmpdir/trimmed"

cleanup() {
    rm -rf "$tmpdir"
}
trap cleanup EXIT

if [ "$cpu" -le 1 ]; then
    echo "[hgene] Running porechop single-thread" >&2
    porechop -i "$input_fastq" -o "$output_fastq" --discard_middle >/dev/null # trimming
    exit 0
fi

echo "[hgene] Running porechop in parallel with $cpu jobs" >&2

# FASTQ must have a line count divisible by 4.
total_lines=$(wc -l < "$input_fastq")
if [ $((total_lines % 4)) -ne 0 ]; then
    echo "ERROR: input FASTQ line count is not divisible by 4: $input_fastq" >&2
    exit 1
fi

# Split while preserving complete FASTQ records.
# Keep a minimum of one read per chunk.
lines_per_chunk=$(( (total_lines / cpu / 4) * 4 ))
if [ "$lines_per_chunk" -lt 4 ]; then
    lines_per_chunk=4
fi

split -d -l "$lines_per_chunk" --additional-suffix=.fastq "$input_fastq" "$tmpdir/chunks/chunk_"

for chunk in "$tmpdir"/chunks/*.fastq; do
    chunk_name=$(basename "$chunk" .fastq)
    (
        porechop -i "$chunk" -o "$tmpdir/trimmed/${chunk_name}_tr.fastq" --discard_middle >/dev/null
    ) &

    # Limit concurrent jobs to $cpu.
    while [ "$(jobs -rp | wc -l)" -ge "$cpu" ]; do
        sleep 0.2
    done
done

wait

cat "$tmpdir"/trimmed/chunk_*_tr.fastq > "$output_fastq"

echo "[hgene] Porechop output written to: $output_fastq" >&2
