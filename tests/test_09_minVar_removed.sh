#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]:-$0}" )" &> /dev/null && pwd )"

cleanup() { rm test_09_in.vcf test_09_out.vcf ; }
trap cleanup EXIT

cat > "test_09_in.vcf" <<'VCF'
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=SB,Number=1,Type=Integer,Description="Strand bias">
##INFO=<ID=DP4,Number=4,Type=Integer,Description="Counts for ref-fwd, ref-rev, alt-fwd, alt-rev">
##INFO=<ID=HRUN,Number=1,Type=Integer,Description="Homopolymer run length">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
UL23-HHV1	103	.	C	A	0	.	DP=99;AF=0.15;SB=0;DP4=0,74,15,0
VCF

python3 "$SCRIPT_DIR/../src/readVCF.py" "test_09_in.vcf" > "test_09_out.vcf"

# La ligne variant ne doit plus être là
if grep -q $'UL23-HHV1\t103\t.\tC\tA' "test_09_out.vcf"; then
  echo "FAIL: minority variant with low depth (DP=99) not removed"
  exit 1
fi

# Et le header doit rester
grep -q "^#CHROM" "test_09_out.vcf" || { echo "FAIL: VCF header missing"; exit 1; }

echo "OK test_09_minVar_removed"

