#Darina's pipeline

#1. QC on the genotype data (removing SNPs/IDs with a low callrate (below 98%), removal of related IDs (one from each pair), checking for sex consistency and mean heterozygosity (usually I exclude all IDs with a heterozygosity > 4sd from the mean) , using MDS to check for possible population outliers, removing SNPs with low MAF (below 1%), removing of AT/CG SNPs, removing of SNPs with HWE p-value < e-05).

#2. Phasing the qc-ed data using shapeit2 (this includes updating the positions to the specific reference sample and checking for strand inconsistencies)

#3. Imputing the shapeit2 data using impute2 (using a chunk size of 5 MB)

#4. QC of the imputed data (info score of at least 0.8, HWE p-value not lower than e-05).

#I usually then convert the imputed data into best guessed genotypes (using a probabiliy threshold of 90%) and then rerun a check on SNP and ID callrates on these data.



## built using https://github.com/Chris1221/impute as a guide

### need to recode AB calls to AGCT for matching 1000 genomes
Recode_for_phasing.R

#You must split the dataset by chromosomes prior to phasing since SHAPEIT phases only one chromosome at a time. To do so, you can use the following Plink command for example:

for chr in $(seq 1 22); do
     plink --noweb --file PAWS_snps_filtered_recoded \
           --chr $chr \
           --recode \
           --out PAWS_snps_filtered_recoded.chr$chr ;
done


## Make bed files out of ped and map

for chr in $(seq 1 22); do
     plink --noweb --file PAWS_snps_filtered_recoded.chr$chr --make-bed --out PAWS_snps_filtered_recoded.chr$chr ;
done



## need to remove the duplicates with an exm name
#This creates a file named exclude_snp_list.txt
Rscript remove_dups.R

## will need to re split after removing duplicates so delete these files
rm PAWS_snps_filtered_recoded*chr*


## exclude duplicates in chromosome split
for chr in $(seq 1 20); do
     plink --noweb --file PAWS_snps_filtered_recoded \
     --exclude exclude_snp_list.chr$chr.txt \
           --chr $chr \
           --recode \
           --out PAWS_snps_filtered_recoded.chr$chr ;
done


## beds again
for chr in $(seq 1 20); do
     plink --noweb --file PAWS_snps_filtered_recoded.chr$chr --make-bed --out PAWS_snps_filtered_recoded.chr$chr ;
done


### Data check against 1000 genomes (with duplicates excluded)
#Now we find any alignment issues and pipe them to two different files, the first will be flipped and the second will be excluded:
                    shapeit -check \
                            -B PAWS_snps_filtered_recoded.chr10\
                            -M ALL.integrated_phase1_SHAPEIT_16-06-14.nomono/genetic_map_chr10_combined_b37.txt \
                            --input-ref \
                            ALL.integrated_phase1_SHAPEIT_16-06-14.nomono/ALL.chr10.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.haplotypes.gz\
                            ALL.integrated_phase1_SHAPEIT_16-06-14.nomono/ALL.chr10.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.legend.gz\
                            ALL.integrated_phase1_SHAPEIT_16-06-14.nomono/ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.sample\
                            --output-log PAWS_snps_filtered_recoded.chr10.alignments

for chr in $(seq 1 20); do
shapeit -check \
        -B PAWS_snps_filtered_recoded.chr${chr}\
        -M ALL.integrated_phase1_SHAPEIT_16-06-14.nomono/genetic_map_chr${chr}_combined_b37.txt \
        --input-ref \
        ALL.integrated_phase1_SHAPEIT_16-06-14.nomono/ALL.chr${chr}.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.haplotypes.gz\
        ALL.integrated_phase1_SHAPEIT_16-06-14.nomono/ALL.chr${chr}.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.legend.gz\
        ALL.integrated_phase1_SHAPEIT_16-06-14.nomono/ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.sample\
        --output-log PAWS_snps_filtered_recoded.chr${chr}.alignments
done




#To obtain the sites which require flipping, we isolate the strand flip issues and pipe #to a new file.

for CHR in $(seq 1 20)
do

cat PAWS_snps_filtered_recoded.chr${CHR}.alignments.snp.strand | grep "Strand" | awk '{ print $4 }' > flip.chr${CHR}.txt

cat PAWS_snps_filtered_recoded.chr${CHR}.alignments.snp.strand | grep "Missing" | awk '{ print $4 }' > exclude.chr${CHR}.txt

done



### Now flip strand (see git hub to expand to all chr)

                        plink --noweb --bfile PAWS_snps_filtered_recoded.chr10 --flip flip.chr10.txt --make-bed --out PAWS_snps_filtered_recoded.chr10.flipped

for chr in $(seq 1 20); do
     plink --noweb --bfile PAWS_snps_filtered_recoded.chr$chr --flip flip.chr$chr.txt --make-bed --out PAWS_snps_filtered_recoded.chr$chr.flipped ;
done


### Data check against 1000 genomes (with flips fixed)
### anything that errors here is a true exclude (it is possible some of the flips can not be fixed)
for chr in $(seq 1 20); do
shapeit -check \
        -B PAWS_snps_filtered_recoded.chr${chr}.flipped\
        -M ALL.integrated_phase1_SHAPEIT_16-06-14.nomono/genetic_map_chr${chr}_combined_b37.txt \
        --input-ref \
        ALL.integrated_phase1_SHAPEIT_16-06-14.nomono/ALL.chr${chr}.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.haplotypes.gz\
        ALL.integrated_phase1_SHAPEIT_16-06-14.nomono/ALL.chr${chr}.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.legend.gz\
        ALL.integrated_phase1_SHAPEIT_16-06-14.nomono/ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.sample\
        --output-log PAWS_snps_filtered_recoded.chr${chr}.flipped.alignments
done

## exclude the true unfixable snps
for chr in $(seq 1 20); do
shapeit -check \
        -B PAWS_snps_filtered_recoded.chr${chr}.flipped\
        -M ALL.integrated_phase1_SHAPEIT_16-06-14.nomono/genetic_map_chr${chr}_combined_b37.txt \
        --input-ref \
        ALL.integrated_phase1_SHAPEIT_16-06-14.nomono/ALL.chr${chr}.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.haplotypes.gz\
        ALL.integrated_phase1_SHAPEIT_16-06-14.nomono/ALL.chr${chr}.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.legend.gz\
        ALL.integrated_phase1_SHAPEIT_16-06-14.nomono/ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.sample\
        --exclude-snp PAWS_snps_filtered_recoded.chr${chr}.flipped.alignments.snp.strand.exclude\
        --output-log PAWS_snps_filtered_recoded.chr${chr}.flipped.excluded.alignments 
done




## FINALLY the actual pre-phasing (same code as before but without the check, output line changed and using multiple threads)
for chr in $(seq 1 20); do
shapeit -B PAWS_snps_filtered_recoded.chr${chr}.flipped\
        -M ALL.integrated_phase1_SHAPEIT_16-06-14.nomono/genetic_map_chr${chr}_combined_b37.txt \
        --input-ref \
        ALL.integrated_phase1_SHAPEIT_16-06-14.nomono/ALL.chr${chr}.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.haplotypes.gz\
        ALL.integrated_phase1_SHAPEIT_16-06-14.nomono/ALL.chr${chr}.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.legend.gz\
        ALL.integrated_phase1_SHAPEIT_16-06-14.nomono/ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.sample\
        --thread 4\
        --exclude-snp PAWS_snps_filtered_recoded.chr${chr}.flipped.alignments.snp.strand.exclude\
        -O PAWS_snps_filtered_recoded.chr${chr}.phased.with.ref
done


## imputation is done in impute.R impute.sh 

Rscript impute.R



### roganize seperate files and cat to one file
ls | grep 'info\|summary\|warnings\|diplotype' | xargs -d"\n" rm 

mkdir ../final
for CHR in $(seq 1 20)
do
    ls | grep "_chr${CHR}.flipped" | xargs -d"\n" cat > ../final/PAWS_psychchip_chr${CHR}.gen
done

#move directory
cd ../final


## quality control (info score of at least 0.8, HWE p-value not lower than e-05)
for chr in $(seq 1 20)
do
qctool -g PAWS_psychchip_chr${chr}.gen -hwe 5 -info 0.8 1 -omit-chromosome -og PAWS_psychchip_chr${chr}.imputed.QC.gen
done


## Convert to ped and map (using a probabiliy threshold of 90%)
for chr in $(seq 1 20)
do
gtool -G --g PAWS_psychchip_chr${chr}.imputed.QC.gen --chr ${chr} --s ~/bootcamp/PsychChip/impute/PAWS_snps_filtered_recoded.chr${chr}.phased.with.ref.sample --ped PAWS_snps_imputed_QC.chr${chr}.ped --map PAWS_snps_imputed_QC.chr${chr}.map --threshold 0.9
done



# convert to bed and then lgen so can bring bim into R

for chr in $(seq 1 22); do
    plink --noweb --file PAWS_snps_imputed_QC.chr$chr --out PAWS_snps_imputed_QC.chr$chr --make-bed ;
done

for chr in $(seq 1 22); do
    plink --noweb --bfile PAWS_snps_imputed_QC.chr$chr --recode-lgen --out PAWS_snps_imputed_QC_recode.chr$chr ;
done

mkdir mQTL_results

## perform mQTL analysis in R
Rscript ~/bootcamp/mQTL/PAWS_mQTL.R
