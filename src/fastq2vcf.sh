#!/bin/bash
#$1 = prefix file
CPU=$2
virus=$3
echo "Starting.."
log="log.txt"
SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]:-$0}"; )" &> /dev/null && pwd 2> /dev/null; )";

echo "logfile; Sample:$1" > $log

echo "Trimming.."
echo "Trimming.." >> $log

porechop -t $CPU --discard_middle -i $1.fastq -o $1_tr.fastq 1> /dev/null # trimming

echo "Map to reference.." >> $log
echo "Map to reference.."
echo $virus >> $log

bash $SCRIPT_DIR/map2HHV.sh $1 $virus $CPU >> $log

echo "Generate consensus.."

bash $SCRIPT_DIR/2vcf.sh $1 $virus $CPU
