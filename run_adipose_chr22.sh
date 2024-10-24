#!/bin/bash


####### CHR 22 #######

TISSUE="adipose"
CHROM=22
START=10519389
END=50807702

i=$START

while [ $i -le $END ]
do
 j=$(( $i + 999999))
 bsub -m "judicator" -R "rusage[mem=4GB]" -o /project/voight_circ-eQTLs/jwoerner/data/out/all/log/${TISSUE}_chr${CHROM}.out Rscript /project/voight_circ-eQTLs/jwoerner/circQTL-tools/circQTL.R $CHROM $i $j $TISSUE /project/voight_circ-eQTLs/jwoerner/data/covariates/${TISSUE}_covariates.txt /project/voight_circ-eQTLs/jwoerner/data/time/${TISSUE}_time.csv /project/voight_circ-eQTLs/jwoerner/data/expression/${TISSUE}_expression.csv /project/voight_circ-eQTLs/jwoerner/data/snps/snps_chr${CHROM}.txt 500000 /project/voight_circ-eQTLs/jwoerner/data/out/all/results/
 i=$(( $i + 1000000 ))
done

####################
