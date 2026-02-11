#!/bin/bash

#$1 prefix $2 virus $3 CPU

virus=$2
CPU=$3

export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib"

SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]:-$0}"; )" &> /dev/null && pwd 2> /dev/null; )";

rm $1_tr.fastq

samtools sort $1.sam -o $1.sorted1.bam
samtools index $1.sorted1.bam
samtools view -bq 40 $1.sorted1.bam > $1.filtered2.bam

if [ "$virus" = "HHV5" ]; then
  samtools view -e 'rlen>99' -O BAM -o $1.filtered.bam $1.filtered2.bam
else
  samtools view -e 'rlen>999' -O BAM -o $1.filtered.bam $1.filtered2.bam
fi

samtools sort $1.filtered.bam -o $1.sorted.bam
samtools index $1.sorted.bam

rm $1.sam $1.sorted1.bam $1.filtered2.bam $1.filtered.bam $1.sorted1.bam.bai

lofreq indelqual -u 16 -f $SCRIPT_DIR/../db/$2.fasta $1.sorted.bam > $1.bam # Add indel quality

rm $1.sorted.bam $1.sorted.bam.bai

samtools index $1.bam

lofreq call-parallel --pp-threads $CPU -f $SCRIPT_DIR/../db/$2.fasta --no-default-filter -A -B -a 1 -b 1 --force-overwrite --call-indels $1.bam -o $1.lofreq.vcf.gz

gunzip $1.lofreq.vcf.gz

lofreq filter --no-defaults --cov-min 20 --af-min 0.1 -i $1.lofreq.vcf -o $1_fl.vcf #--sb-alpha 0.01 --sb-incl-indels

# Contre la référence
python3 $SCRIPT_DIR/readVCF.py $1_fl.vcf > $1.3.vcf # on applique un threshold different au deletion (car homopolymer biais -> différent selon la longueur de la region HRUN= dans le .vcf)

bgzip -f $1.3.vcf

bcftools index $1.3.vcf.gz

#norm mettre deux lignes pour deux variants même positions

#csq -> conséquences AA
bcftools norm -Ou -m- -f $SCRIPT_DIR/../db/$2.fasta $1.3.vcf.gz | bcftools csq -Ob --local-csq --force -o $1.4.vcf.gz \
              --fasta-ref=$SCRIPT_DIR/../db/$2.fasta \
              --gff=$SCRIPT_DIR/../db/$2.gff

#Fusion des SNV majoritaires dans un même codon
python3 $SCRIPT_DIR/merge_codon_mutations.py  $1.4.vcf.gz $SCRIPT_DIR/../db/$2.fasta $SCRIPT_DIR/../db/$2.gff 1 $1.bam 0.9 0.1 > $1.vcf

bgzip -f $1.vcf
bcftools index $1.vcf.gz



rm $1.3.vcf.gz $1.3.vcf.gz.csi

mkdir BAM_$1
mv $1.bam BAM_$1/
mv $1.bam.bai BAM_$1/
mkdir VCF_$1
mv $1_fl.vcf VCF_$1/
mv $1.lofreq.vcf VCF_$1/
mv $1.lofreq.vcf.gz.tbi VCF_$1/
mv $1.4.vcf.gz VCF_$1/
mkdir $1_output
mv VCF_$1 $1_output
mv BAM_$1 $1_output
mv $1.vcf.gz $1_output
mv $1.vcf.gz.csi $1_output
mv $1.log $1_output
