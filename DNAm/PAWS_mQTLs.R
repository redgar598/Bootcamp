#Correlation of methylation and SNPs

### correlation with methylation


setwd("~/bootcamp/mQTL")

library(reshape)
library(ggplot2)
library(RColorBrewer)
library(lme4)
library(gridExtra)



## load meth data beacuse doesnt need to be repeated in each function run
load("PAWS_betas_metas_normalized_combatted_SES_genoPCs_cell_comp.RData")
load("~/Price_annotation.RData")
#filter to vCpGS

# this data may change to the VMR CpGs only
load("PAWS_meth_varprobes_2016Oct13.RData")
PAWS_Beta<-PAWS_Beta[which(rownames(PAWS_Beta)%in%rownames(PAWS_meth_varprobes)),]
PAWS_Beta<-as.data.frame(PAWS_Beta)#20494   177

#M value transformation
Mval<-function(beta) log2(beta/(1-beta))
PAWS_Mval = apply(PAWS_Beta, 1, Mval) # need mvalues for combat
PAWS_Mval = as.data.frame(PAWS_Mval)
PAWS_Mval = t(PAWS_Mval)

# Betas  20494    177





## SNPs

lapply(1:20, function(chr){

    imputed_snps<-read.table(paste0("~/bootcamp/PsychChip/PAWS_snps_imputed_QC.chr",chr,"_recode.lgen"))
    colnames(imputed_snps)<-c("Sample_ID","family_ID","SNP_Name","Allele1","Allele2")
    imputed_snps_map<-read.table(paste0("~/bootcamp/PsychChip/PAWS_snps_imputed_QC.chr",chr,"_recode.map"))
    colnames(imputed_snps_map)<-c("Chr","SNP_Name","cM","Position")
    
          
    ## recode to 0,1,2 (11=0, 12=1, 22=2) nrow(imputed_snps)
    Genotype<-sapply(1:nrow(imputed_snps), function(x) {
      if(imputed_snps$Allele1[x]==1){
        if(imputed_snps$Allele2[x]==1){0}else{
          if(imputed_snps$Allele2[x]==2){1}else{NA}}}else{
          
          if(imputed_snps$Allele1[x]==2){
            if(imputed_snps$Allele2[x]==2){2}else{
              if(imputed_snps$Allele2[x]==1){1}else{NA}}}else{NA}}      
        })
    
    
    imputed_snps$Genotype<-Genotype
    imputed_snps_mini<-imputed_snps[, c(1,3,6)]
    #many snp names where . and need to be removed
    imputed_snps_mini_filtered<-imputed_snps_mini[which(imputed_snps_mini$SNP_Name!="."),]
    #dim(x<-duplicated(imputed_snps_mini_filtered[,1:2]))
    
    
    imputed_snps_cast<-cast(imputed_snps_mini_filtered, SNP_Name~Sample_ID) #33039   193
    imputed_snps_cast<-as.data.frame(imputed_snps_cast)
    
    save(imputed_snps_cast, file=paste0("PAWS_SNPs_Genotypes",chr,".RData"))
    
    
    
    
    #Heterozygosity SNPs
    
    Heterozygosity_levels<-lapply(1:nrow(imputed_snps_cast), function(x) {
      table(unlist(imputed_snps_cast[x,2:ncol(imputed_snps_cast)]))
      })
    hets<-sapply(1:length(Heterozygosity_levels), function(x) length(Heterozygosity_levels[[x]]))#12884 no varibility in PAWS 10924 Hets and Homos 9231 hets and 2 Homos
    
    table(hets)
    
    imputed_snps_variable<-imputed_snps_cast[which(hets>1),]
    
    ## rename mixed sample
    colnames(imputed_snps_variable)[which(colnames(imputed_snps_variable)=="10_9630002082_R10C02")]<-"186"
    colnames(imputed_snps_variable)[which(colnames(imputed_snps_variable)=="10_9630002082_R08C02")]<-"10"
    
    save(imputed_snps_variable, file=paste0("imputed_snps_variable_",chr,".RData"))
    
    
    
    #add annoation information
    
    imputed_snps_map<-imputed_snps_map[which(imputed_snps_map$SNP_Name%in%imputed_snps_variable$SNP_Name),]
    
    print(nrow(imputed_snps_map)) #16377 original SNPs, and 120436 total (104059 imputed)
    
    
    
    ## which CpGs are within 5kb of a SNP    what distance to use????????????????
    
    #filter to vCpGs on the chr of interest
    annotation$CpG<-rownames(annotation)
    annotation_simple<-annotation[,c("CpG","MAPINFO","CHR")]
    annotation_simple_chr<-annotation_simple[which(annotation_simple$CHR==chr),]
    
    PAWS_vCpGs<-as.data.frame(annotation_simple_chr[which(annotation_simple_chr$CpG%in%rownames(PAWS_Beta)),])
    
    imputed_snps_map$SNP_Name<-as.character(imputed_snps_map$SNP_Name)
    PAWS_vCpGs$CHR<-as.character(PAWS_vCpGs$CHR)
    
    
    
    ## Matching table
    CpG_SNP_realtions<-lapply(1:nrow(imputed_snps_map), function(x) {
        coor<-imputed_snps_map[x,"Position"]
        CpGs<-PAWS_vCpGs[which(PAWS_vCpGs$MAPINFO>(coor-5000) & PAWS_vCpGs$MAPINFO<(coor+5000)),]
        if(nrow(CpGs)>0){
          data.frame(SNP=imputed_snps_map[x,"SNP_Name"], CpGs=rownames(CpGs), distance=(coor-CpGs$MAPINFO))}
      })
    
    CpG_SNP_realtions<-do.call(rbind, CpG_SNP_realtions)
    save(CpG_SNP_realtions, file=paste0("CpG_SNP_realtions_chr",chr,".RData"))
    
    print(paste0(length(unique(CpG_SNP_realtions$SNP))," SNPs and ", length(unique(CpG_SNP_realtions$CpGs)), " CpGs to associated on chr ", chr))
    
    
    
  
    ## mQTLs
    
    ## match SNPS, CpGs, and meta data
     rownames(imputed_snps_variable)<-imputed_snps_variable$SNP_Name
    imputed_snps_variable$SNP_Name<-NULL
    
    imputed_snps_variable<-imputed_snps_variable[,which(colnames(imputed_snps_variable)%in%PAWS_meta$PAWSG_ID)]
    imputed_snps_variable<-imputed_snps_variable[,match(PAWS_meta$PAWSG_ID, colnames(imputed_snps_variable))]
    
    PAWS_Mval<-PAWS_Mval[,which(colnames(PAWS_Mval)%in%PAWS_meta$SAMPLE_ID)]
    PAWS_Mval<-PAWS_Mval[,match(PAWS_meta$SAMPLE_ID, colnames(PAWS_Mval))]
    
    
    
    mQTL_snp_cpg<-lapply(1:length(unique(CpG_SNP_realtions$SNP)), function(snp){ #
      #print(snp)
      snp=unique(CpG_SNP_realtions$SNP)[snp]
      snp_genotypes<-imputed_snps_variable[which(rownames(imputed_snps_variable)==snp),]
      CpGs<-CpG_SNP_realtions[which(CpG_SNP_realtions$SNP==snp),]
    
      CpGs_methylation<-PAWS_Mval[which(rownames(PAWS_Mval)%in%CpGs$CpGs),]
      
      cpg_num<-if(class(CpGs_methylation)=="numeric"){1}else{nrow(CpGs_methylation)}
      
      CpG_mods<-lapply(1:cpg_num, function(cpg){
          #add geno and DNAm to meta for lm
          PAWS_meta_ex<-PAWS_meta
          if(class(CpGs_methylation)=="numeric"){PAWS_meta_ex$DNAm<-CpGs_methylation}else{PAWS_meta_ex$DNAm<-CpGs_methylation[cpg,]}
          PAWS_meta_ex$SNP<-unlist(snp_genotypes)
        
        mod<-lm(DNAm~AGE+as.factor(ReportedSex)+geno_PC1+geno_PC2+SNP, data=PAWS_meta_ex)
        pval<-if("SNP"%in%rownames(coef(summary(mod)))){coef(summary(mod))["SNP",4]}else{NA}
        aic<-AIC(mod)
        coef<-coef(mod)["SNP"]
        
        data.frame(AIC=aic, pval=pval, coef=coef)
        })
    
    CpG_mods<-do.call(rbind, CpG_mods)
    if(class(CpGs_methylation)=="numeric"){CpG_mods$CpG<-CpGs$CpGs}else{CpG_mods$CpG<-rownames(CpGs_methylation)}
    CpGs_mods<-merge(CpGs, CpG_mods, by.x="CpGs", by.y="CpG")
      })
    
    
    mQTL_models<-do.call(rbind, mQTL_snp_cpg)
    
    save(mQTL_models, file=paste0("mQTL_models_chr",chr,".RData")) 
    
    })






