#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]:-$0}" )" &> /dev/null && pwd )"

cleanup() { rm test_10_in.vcf test_10_out.vcf ; }
trap cleanup EXIT

cat > "test_10_in.vcf" <<'VCF'
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=SB,Number=1,Type=Integer,Description="Strand bias">
##INFO=<ID=DP4,Number=4,Type=Integer,Description="Counts for ref-fwd, ref-rev, alt-fwd, alt-rev">
##INFO=<ID=HRUN,Number=1,Type=Integer,Description="Homopolymer run length">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
UL23-HHV1	103	.	C	A	0	.	DP=100;AF=0.15;SB=0;DP4=75,0,15,0
VCF

python3 "$SCRIPT_DIR/../src/hgene_filter_custom.py" "test_10_in.vcf" > "test_10_out.vcf"

grep -q "^#CHROM" "test_10_out.vcf" || { echo "FAIL: VCF header missing"; exit 1; }

grep -q $'UL23-HHV1\t103\t.\tC\tA' "test_10_out.vcf" || {   echo "FAIL: minority variant (DP=100) not kept"; exit 1; }

echo "OK test_10_minVar_kept"

