#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]:-$0}" )" &> /dev/null && pwd )"

cleanup() { rm test_03_in.vcf test_03_out.vcf ; }
trap cleanup EXIT

# VCF minimal avec une insertion en homopolymère à 35% (AF=0.35, HRUN=8)
cat > "test_03_in.vcf" <<'VCF'
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=SB,Number=1,Type=Integer,Description="Strand bias">
##INFO=<ID=DP4,Number=4,Type=Integer,Description="Counts for ref-fwd, ref-rev, alt-fwd, alt-rev">
##INFO=<ID=HRUN,Number=1,Type=Integer,Description="Homopolymer run length">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
UL23-HHV1	103	.	CA	C	0	.	DP=7996;AF=0.35;SB=0;DP4=3000,2180,1400,1416;INDEL;HRUN=8
VCF

python3 "$SCRIPT_DIR/../src/readVCF.py" "test_03_in.vcf" > "test_03_out.vcf"

# La ligne variant ne doit plus être là
if grep -q $'UL23-HHV1\t103\t.\tCA\tC' "test_03_out.vcf"; then
  echo "FAIL: homopolymer insertion (AF=0.35, HRUN=8) not removed"
  exit 1
fi

# Et le header doit rester
grep -q "^#CHROM" "test_03_out.vcf" || { echo "FAIL: VCF header missing"; exit 1; }

echo "OK test_03_hpoly_insertion_removed"

