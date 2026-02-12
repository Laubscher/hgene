#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]:-$0}" )" &> /dev/null && pwd )"

cleanup() { rm test_07_in.vcf test_07_out.vcf ; }
trap cleanup EXIT

cat > "test_07_in.vcf" <<'VCF'
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=SB,Number=1,Type=Integer,Description="Strand bias">
##INFO=<ID=DP4,Number=4,Type=Integer,Description="Counts for ref-fwd, ref-rev, alt-fwd, alt-rev">
##INFO=<ID=HRUN,Number=1,Type=Integer,Description="Homopolymer run length">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
UL23-HHV1	103	.	CA	C	0	.	DP=7996;AF=0.35;SB=0;DP4=3000,2180,1400,1416;INDEL
VCF

python3 "$SCRIPT_DIR/../src/hgene_filter_custom.py" "test_07_in.vcf" > "test_07_out.vcf"

grep -q "^#CHROM" "test_07_out.vcf" || { echo "FAIL: VCF header missing"; exit 1; }

grep -q $'UL23-HHV1\t103\t.\tCA\tC' "test_07_out.vcf" || {   echo "FAIL: deletion (AF=0.35, no HRUN field) not kept"; exit 1; }

echo "OK test_07_noHpoly_deletion_kept"

