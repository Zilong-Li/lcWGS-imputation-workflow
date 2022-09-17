#!/usr/bin/env sh

set -e

if [ "$#" -ne 4 ]; then
    echo "Illegal number of parameters"
    echo "Usage: $0 prefix in-vcf out-dir samples"
    exit
fi

PREFIX=$1
IN=$2
OUTDIR=$3
SAMPLES=$4
THREADS=4
OUT1=$OUTDIR/$PREFIX.bcf ## vcf for glimpse
OUT2=$OUTDIR/$PREFIX     ## prefix of hap legend sample

# bcftools norm -m -any $IN -Ou --threads 4 |  bcftools view -v snps -m2 -M2 -s $SAMPLES --threads 4 | bcftools norm -d snps -Ob -o $OUT1 --threads 4 && bcftools index -f $OUT1
## glimpse accept multiallelics? but it seems quilt says no. so use SNP only
bcftools view -v snps -m2 -M2 -s $SAMPLES --threads $THREADS $IN | bcftools norm - -d snps -Ob -o $OUT1 --threads $THREADS && bcftools index -f $OUT1
bcftools convert --haplegendsample $OUT2 $OUT1

OUT3=$OUTDIR/$PREFIX.sites.vcf.gz
bcftools view -G $OUT1 -Oz -o $OUT3 && tabix -f $OUT3

echo "$@ done"
