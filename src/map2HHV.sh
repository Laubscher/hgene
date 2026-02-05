#!/bin/bash
#$1 prefix file name
#$2 Herpes ; HHV1 HHV2 VZV(HHV3) CMV(HHV5)

SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]:-$0}"; )" &> /dev/null && pwd 2> /dev/null; )";

minimap2 -ax map-ont $SCRIPT_DIR/../db/$2.fasta $1_tr.fastq > $1.sam
