
#!/bin/bash
#$ -S /bin/bash
#$ -q abaqus.q
#$ -l qname=abaqus.q
#$ -cwd
#$ -V
#$ -l mf=192G
#$ -j y
#$ -o /home/hpc2862/repos/impute/logs/$JOB_NAME.txt

#!/bin/bash

#all other options come from config

CHR=$1
START=$2
STOP=$3
OPTION=$4


cd ~/bootcamp/PsychChip/impute

impute2 \
-known_haps_g PAWS_snps_filtered_recoded.chr${CHR}.phased.with.ref.haps  \
-h ALL.integrated_phase1_SHAPEIT_16-06-14.nomono/ALL.chr${CHR}.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.haplotypes.gz \
-l ALL.integrated_phase1_SHAPEIT_16-06-14.nomono/ALL.chr${CHR}.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.legend.gz \
-m ALL.integrated_phase1_SHAPEIT_16-06-14.nomono/genetic_map_chr${CHR}_combined_b37.txt \
-int ${START} ${STOP} \
-Ne 20000 \
-buffer 250 \
-o PAWS_psychchip_chr${CHR}.flipped.phased.imputed.${START}.${STOP} ${OPTION}


