setwd("/big_data/bootcamp/PAWS/mQTL")

library(reshape)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(reshape2)
library(dplyr)  


## load meth data beacuse doesnt need to be repeated in each function run
load("PAWS_betas_metas_normalized_combatted_SES_genoPCs_cell_comp.RData")
load("Price_annotation.RData")
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


# 
# ## Combine CpG_relations and SNP data
# 
# load("mQTL_models_5kb_chr1.RData")
# 
# mQTL_models_all_chr<-mQTL_models
# 
# lapply(2:20, function(chr){
#   load(paste0("mQTL_models_5kb_chr",chr,".RData"))
#   mQTL_models_all_chr<<-rbind(mQTL_models_all_chr, mQTL_models)
# })
# save(mQTL_models_all_chr, file="mQTL_models_5kb_all_chr.RData" )
# 
# setwd("/big_data/bootcamp/PAWS/mQTL")
# load("imputed_snps_variable_1.RData")
# 
# imputed_snps_variable_all_chr<-imputed_snps_variable
# 
# invisible(lapply(2:20, function(chr){
#   print(chr)
#   load(paste0("imputed_snps_variable_",chr,".RData"))
#   imputed_snps_variable_all_chr<<-rbind(imputed_snps_variable_all_chr, imputed_snps_variable)
#   rm(imputed_snps_variable)
# }))
# 
# save(imputed_snps_variable_all_chr, file="imputed_snps_variable_all_chr.RData" )
# 
# 


### BEST AIC
load("imputed_snps_variable_all_chr.RData")
load("mQTL_models_5kb_all_chr.RData")

best_model<-mQTL_models_all_chr %>% 
  group_by(CpGs) %>% 
  slice(which.min(AIC)) %>%
  as.data.frame()

colnames(best_model)<-c("CpG", "SNP","distance","AIC_mQTL","pval_mQTL","coef_mQTL")



#### interaction model
rownames(imputed_snps_variable_all_chr)<-imputed_snps_variable_all_chr$SNP_Name
imputed_snps_variable_all_chr$SNP_Name<-NULL

imputed_snps_variable_all_chr<-imputed_snps_variable_all_chr[,which(colnames(imputed_snps_variable_all_chr)%in%PAWS_meta$PAWSG_ID)]
imputed_snps_variable_all_chr<-imputed_snps_variable_all_chr[,match(PAWS_meta$PAWSG_ID, colnames(imputed_snps_variable_all_chr))]

PAWS_Mval<-PAWS_Mval[,which(colnames(PAWS_Mval)%in%PAWS_meta$SAMPLE_ID)]
PAWS_Mval<-PAWS_Mval[,match(PAWS_meta$SAMPLE_ID, colnames(PAWS_Mval))]

print(paste0(nrow(best_model), " models to run"))

GxE_snp_cpg<-lapply(1:nrow(best_model), function(x){ #
  #print(snp)
  snp=best_model$SNP[x]
  snp_genotypes<-imputed_snps_variable_all_chr[which(rownames(imputed_snps_variable_all_chr)==snp),]
  CpG<-best_model$CpG[x]
  
  CpG_methylation<-PAWS_Mval[which(rownames(PAWS_Mval)==CpG),]
  
    #add geno and DNAm to meta for lm
    PAWS_meta_ex<-PAWS_meta
    PAWS_meta_ex$SNP<-unlist(snp_genotypes)
    PAWS_meta_ex$DNAm<-unlist(CpG_methylation)
    
    mod<-lm(DNAm~PredictedCD34+PredictedBuccal+AGE+as.factor(ReportedSex)+geno_PC1+geno_PC2+Sarah_comp_SES+SNP+SNP*Sarah_comp_SES, data=PAWS_meta_ex)
    pval_GxE<-if("Sarah_comp_SES:SNP"%in%rownames(coef(summary(mod)))){coef(summary(mod))["Sarah_comp_SES:SNP",4]}else{NA}
    aic_GxE<-AIC(mod)
    coef_GxE<-coef(mod)["Sarah_comp_SES:SNP"]
    
    data.frame(CpG=CpG,AIC_GxE=aic_GxE, pval_GxE=pval_GxE, coef_GxE=coef_GxE)
  })
  
GxE_snp_cpg<-do.call(rbind, GxE_snp_cpg)

best_model<-merge(best_model, GxE_snp_cpg, by="CpG")

## Add EWAS
load("/big_data/bootcamp/PAWS/PAWS_EWAS_results_2016Oct13.RData")
PAWS_EWAS_results<-as.data.frame(t(as.data.frame(PAWS_EWAS_results)))
colnames(PAWS_EWAS_results)<-c("coef_EWAS","Std.Error","t value","pval_EWAS","AIC_EWAS")
PAWS_EWAS_results$CpG<-unlist(rownames(PAWS_EWAS_results))

best_model<-merge(best_model, PAWS_EWAS_results[,c(6,5,4,1)], by="CpG")
best_model<-best_model[,c(1,2,3,10:12,4:9)]

save(best_model, file="PAWS_best_model_EWAS_mQTL_MGE.RData")
