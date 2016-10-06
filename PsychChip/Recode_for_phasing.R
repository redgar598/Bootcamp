######################
## Recode for Hapmap inclusion
######################
setwd("~/bootcamp/PsychChip/")


PAWS_snps<-read.table("PAWS_Forward_FinalReport.txt", skip=10, sep="\t") #588454 SNPs per sample, need to parse
colnames(PAWS_snps)<-c("SNP_Name","Sample_ID","Allele1_TOP","Allele2_TOP","GC_Score","Allele1_fwd","Allele2_fwd","Position","Chr","B_Allele_Freq","Allele1_AB","Allele2_AB")

            #annotation<-read.csv("Psych_array_annotation.csv", sep="\t", header=T)
            
            #PsychChip_annotation<-read.table("PsychArray_A_annotated.txt",sep="\t", header=T) #588454 SNPs per sample, need to parse

# read in the TOP?BOT strand info
Psychchip_manifest<-read.csv("PsychArray-B.csv", skip=7,  header=T) #558741/571078 have a TOP BOT id :)
Psychchip_manifest_simple<-Psychchip_manifest[,c(2,3,4,5,7,9,10,11,16)]

            # annotation[grep("1KG_1_109440678",annotation$Name),]
            # Psychchip_manifest_simple[grep("1KG_1_109440678",Psychchip_manifest_simple$Name),]


### Explore one test SNP
PsychChip_annotation[which(PsychChip_annotation$Name=="rs1000016"),]
pawsexample<-PAWS_snps[which(PAWS_snps$SNP_Name=="rs1000016"),]
Psychchip_manifest_simple[which(Psychchip_manifest_simple$Name=="rs1000016"),]

## only want SNPs with a top bot (558741/571078 have a TOP BOT id :))
Psychchip_manifest_simple$SNP<-as.character(Psychchip_manifest_simple$SNP)
Psychchip_manifest_simple<-Psychchip_manifest_simple[grep("/", Psychchip_manifest_simple$SNP),]
Psychchip_manifest_simple<-Psychchip_manifest_simple[grep("A|T|G|C", Psychchip_manifest_simple$SNP),] # 558741 SNPs

## filter PAWS data to just those with top bottom
PAWS_snps_filtered<-PAWS_snps[which(PAWS_snps$SNP_Name%in%Psychchip_manifest_simple$Name),]


### recode method
Psychchip_manifest_simple_TOP<-Psychchip_manifest_simple[which(Psychchip_manifest_simple$IlmnStrand=="TOP"),] # 282263 SNPs
PAWS_snps_filtered_TOP<-PAWS_snps_filtered[which(PAWS_snps_filtered$SNP_Name%in%Psychchip_manifest_simple_TOP$Name),]
SNPs_inPAWS_TOP<-as.character(unique(PAWS_snps_filtered_TOP$SNP_Name))# 281375 SNPs

PAWS_snps_filtered_TOP_merge<-merge(PAWS_snps_filtered_TOP, Psychchip_manifest_simple_TOP, by.x="SNP_Name", by.y="Name")


#pull the first and second letter in [/] dependent on A or B call
PAWS_snps_filtered_TOP_merge$Allele1<-sapply(1:nrow(PAWS_snps_filtered_TOP_merge), function(x){
  if(PAWS_snps_filtered_TOP_merge$Allele1_AB[x]=="A"){
    strsplit(as.character(PAWS_snps_filtered_TOP_merge$SNP[x]), "")[[1]][2]}else{
      if(PAWS_snps_filtered_TOP_merge$Allele1_AB[x]=="B"){
        strsplit(as.character(PAWS_snps_filtered_TOP_merge$SNP[x]), "")[[1]][4]}else{
          "NoTOPBOT"}}})

PAWS_snps_filtered_TOP_merge$Allele2<-sapply(1:nrow(PAWS_snps_filtered_TOP_merge), function(x){
  if(PAWS_snps_filtered_TOP_merge$Allele2_AB[x]=="A"){
    strsplit(as.character(PAWS_snps_filtered_TOP_merge$SNP[x]), "")[[1]][2]}else{
      if(PAWS_snps_filtered_TOP_merge$Allele2_AB[x]=="B"){
        strsplit(as.character(PAWS_snps_filtered_TOP_merge$SNP[x]), "")[[1]][4]}else{
          "NoTOPBOT"}}})

save(PAWS_snps_filtered_TOP_merge, file="PAWS_recoded_TOP_only.RData")




### recode method
Psychchip_manifest_simple_BOT<-Psychchip_manifest_simple[which(Psychchip_manifest_simple$IlmnStrand=="BOT"),] # 276478 SNPs
PAWS_snps_filtered_BOT<-PAWS_snps_filtered[which(PAWS_snps_filtered$SNP_Name%in%Psychchip_manifest_simple_BOT$Name),]
SNPs_inPAWS_BOT<-as.character(unique(PAWS_snps_filtered_BOT$SNP_Name))# 275569 SNPs

PAWS_snps_filtered_BOT_merge<-merge(PAWS_snps_filtered_BOT, Psychchip_manifest_simple_BOT, by.x="SNP_Name", by.y="Name")


#pull the first and second letter in [/]
PAWS_snps_filtered_BOT_merge$Allele1<-sapply(1:nrow(PAWS_snps_filtered_BOT_merge), function(x){
  if(PAWS_snps_filtered_BOT_merge$Allele1_AB[x]=="A"){
    strsplit(as.character(PAWS_snps_filtered_BOT_merge$SNP[x]), "")[[1]][4]}else{
      if(PAWS_snps_filtered_BOT_merge$Allele1_AB[x]=="B"){
        strsplit(as.character(PAWS_snps_filtered_BOT_merge$SNP[x]), "")[[1]][2]}else{
          "NoTOPBOT"}}})

PAWS_snps_filtered_BOT_merge$Allele2<-sapply(1:nrow(PAWS_snps_filtered_BOT_merge), function(x){
  if(PAWS_snps_filtered_BOT_merge$Allele2_AB[x]=="A"){
    strsplit(as.character(PAWS_snps_filtered_BOT_merge$SNP[x]), "")[[1]][4]}else{
      if(PAWS_snps_filtered_BOT_merge$Allele2_AB[x]=="B"){
        strsplit(as.character(PAWS_snps_filtered_BOT_merge$SNP[x]), "")[[1]][2]}else{
          "NoTOPBOT"}}})

save(PAWS_snps_filtered_BOT_merge, file="PAWS_recoded_BOT_only.RData")




####### Combine
load("PAWS_recoded_BOT_only.RData")
load("PAWS_recoded_TOP_only.RData")

PAWS_snps_filtered_recoded<-rbind(PAWS_snps_filtered_TOP_merge, PAWS_snps_filtered_BOT_merge)
PAWS_snps_filtered_recoded<-PAWS_snps_filtered_recoded[,c(1,2,5,8:14,21,22)]
## order by sample ID
PAWS_snps_filtered_recoded<-PAWS_snps_filtered_recoded[order(PAWS_snps_filtered_recoded$Sample_ID),] #556944 SNPs

save(PAWS_snps_filtered_recoded, file="PAWS_snps_filtered_recoded.RData")


#########################################################################################################
# Quality Control
## want to filter probes based on Low GC score as QC step

library(reshape)
load("PAWS_snps_filtered_recoded.RData")


## this was all already run in PAWS_SNP_final_report_Loading
# PAWS_snps_GC_Score<-PAWS_snps_filtered_recoded[,c(1,2,3)]
# 
# PAWS_snps_GC_Score_casted<-cast(PAWS_snps_GC_Score, SNP_Name~Sample_ID)
# save(PAWS_snps_GC_Score_casted, file="PAWS_snps_GCScore_recoded.RData")
# 
# load("PAWS_snps_GCScore_recoded.RData")
# 
# # Bad probes (samples have already been Qc'ed above)
# Avg_GC_Score<-rowMeans(PAWS_snps_GC_Score_casted[,2:ncol(PAWS_snps_GC_Score_casted)], na.rm=T)
# save(Avg_GC_Score, file="Avg_GC_Score_recoded.RData")
# 
# 
# #10% GC
# #p10_SNP - 10th percentile GenCall score over all samples for this SNP (like GS did for sample)
# quantile(as.numeric(PAWS_snps_GC_Score_casted[2:ncol(PAWS_snps_GC_Score_casted),1]), c(0.1)) #10% of row
# p10_SNP<-as.vector(sapply(1:nrow(PAWS_snps_GC_Score_casted), function(SNP) quantile(as.numeric(PAWS_snps_GC_Score_casted[SNP,2:ncol(PAWS_snps_GC_Score_casted)]), c(0.1), na.rm=T)))
# save(p10_SNP, file="p10_SNP_recoded.RData")
# 
# Filter_p10<-as.character(PAWS_snps_GC_Score_casted$SNP_Name[which(p10_SNP>0.4 | p10_SNP<0.355)]) # rm 193 Probes w low 10%GC
# Filter_avg<-as.character(PAWS_snps_GC_Score_casted$SNP_Name[which(Avg_GC_Score<0.1)]) # rm 128 Probes w low  average GC 
# 
# #Filter SNPS on NA chr
# PAWS_snps_chr<-unique(PAWS_snps_filtered_recoded[,c(1,5)])
# Filter_NA<-as.character(PAWS_snps_chr$SNP_Name[which(is.na(PAWS_snps_chr$Chr.x))]) #32350
# 
# 
# Filter<-unique(c(Filter_avg, Filter_p10,Filter_NA))# total 32458


load("Quality_control_removed_SNPs.RData")

#########################################################################################################
#PLINK ped amd map

# filter QC'ed
PAWS_snps_recoded_filtered<-PAWS_snps_filtered_recoded[which(!(PAWS_snps_filtered_recoded$SNP_Name%in%Filter)),]# 521067 SNPs
save(PAWS_snps_recoded_filtered, file="PAWS_snps_recoded_filtered.RData")

## Make map and ped from this paws file
load("PAWS_snps_recoded_filtered.RData")

### change "NoTOPBOT" to 0
PAWS_snps_recoded_filtered$Allele1<-sapply(1:nrow(PAWS_snps_recoded_filtered), function(x) if(PAWS_snps_recoded_filtered$Allele1[x]=="NoTOPBOT"){0}else{PAWS_snps_recoded_filtered$Allele1[x]})
PAWS_snps_recoded_filtered$Allele2<-sapply(1:nrow(PAWS_snps_recoded_filtered), function(x) if(PAWS_snps_recoded_filtered$Allele2[x]=="NoTOPBOT"){0}else{PAWS_snps_recoded_filtered$Allele2[x]})

save(PAWS_snps_recoded_filtered, file="PAWS_snps_recoded_filtered.RData")


# Parse MAP for PLINK
#only going to do a 3 column map file so need to specify plink --file mydata --map3
PAWS_snps_filtered_map<-PAWS_snps_recoded_filtered[,c(1,4,5)]
PAWS_snps_filtered_map_unique<-unique(PAWS_snps_filtered_map) # 524486 SNPs
PAWS_snps_filtered_map_unique$cM<-0
PAWS_snps_filtered_map_unique<-PAWS_snps_filtered_map_unique[,c(3,1,4,2)]
PAWS_snps_filtered_map_unique<-PAWS_snps_filtered_map_unique[with(PAWS_snps_filtered_map_unique, order(Position)),]

### Fix the chromosome names in the MAP file because Genome studio annotated the chr wrong (See Rachel, Mina or Sarah G. for an explination)
### the annotation file was obtain from the Illumina website for the Psych array
annotation<-read.csv("Psych_array_annotation.csv", sep="\t", header=T)
length(which(PAWS_snps_filtered_map_unique$SNP_Name%in%annotation$Name))

## note discrepency
annotation[annotation$Name=="exm1272993",]
PAWS_snps_filtered_map_unique[PAWS_snps_filtered_map_unique$SNP_Name=="exm1272993",]

PAWS_snps_filtered_map_unique2<-merge(PAWS_snps_filtered_map_unique, annotation[,c("Chr","Name","MapInfo")], by.x="SNP_Name", by.y="Name")
PAWS_snps_filtered_map_unique2<-PAWS_snps_filtered_map_unique2[match(PAWS_snps_filtered_map_unique$SNP_Name, PAWS_snps_filtered_map_unique2$SNP_Name),]
PAWS_snps_filtered_map_unique2<-PAWS_snps_filtered_map_unique2[,c("Chr","SNP_Name","cM","MapInfo")]

write.table(as.matrix(PAWS_snps_filtered_map_unique2), "PAWS_snps_filtered_recoded_map.txt", sep="\t", row.names = F,col.names = F,quote=FALSE)





# Parse PED for PLINK
library(reshape2)
PAWS_snps_recoded_filtered<-PAWS_snps_recoded_filtered[with(PAWS_snps_recoded_filtered, order(Sample_ID, Position)),]
PAWS_snps_recoded_filtered_ped<-PAWS_snps_recoded_filtered[,c(1,2,11,12)]





#shape them datas to PED
PAWS_snps_recoded_filtered_ped_melt <- melt(PAWS_snps_recoded_filtered_ped, id.vars = c("SNP_Name","Sample_ID"))

#maintain order
PAWS_snps_recoded_filtered_ped_melt_dcast<-dcast(PAWS_snps_recoded_filtered_ped_melt,  
                                                 Sample_ID~factor(SNP_Name,levels=PAWS_snps_recoded_filtered_ped_melt$SNP_Name[1:nrow(PAWS_snps_recoded_filtered_ped_melt)]) + variable)



#Sex meta
meta<-read.csv("PAWS-GdatasetforKoborlabJuly2014-ID_ordered.csv")

# add replicate to meta data
meta_10rep<-subset(meta, PAWSG_ID=="10")
meta_10rep$PAWSG_ID<-"10_9630002082_R08C02"
meta$PAWSG_ID[which(meta$PAWSG_ID=="10")]<-"10_9630002082_R10C02"
meta<-rbind(meta, meta_10rep)
# 193 is not in meta data, can identify from genotyping later
meta_193<-subset(meta, PAWSG_ID=="192")
meta_193[1,]<-NA
meta_193$PAWSG_ID<-"193"
meta_193$Gender_final<-193

meta<-rbind(meta, meta_193)

meta_ordered<-meta[which(meta$PAWSG_ID%in%PAWS_snps_recoded_filtered_ped_melt_dcast$Sample_ID),]
meta_ordered<-meta_ordered[match(PAWS_snps_recoded_filtered_ped_melt_dcast$Sample_ID, meta_ordered$PAWSG_ID),]



# other ped columns
PAWS_snps_recoded_filtered_ped_melt_dcast$Fam_ID<-"0"
PAWS_snps_recoded_filtered_ped_melt_dcast$Pat_ID<-"0"
PAWS_snps_recoded_filtered_ped_melt_dcast$Mat_ID<-"0"
Gender<-sapply(1:nrow(meta_ordered), function(x) if(meta_ordered$Gender_final[x]==0){1}else{if(meta_ordered$Gender_final[x]==1){2}else{0}})

PAWS_snps_recoded_filtered_ped_melt_dcast$Sex<-as.character(Gender)
PAWS_snps_recoded_filtered_ped_melt_dcast$Affection<-"0"
PAWS_snps_recoded_filtered_ped_melt_dcast<-PAWS_snps_recoded_filtered_ped_melt_dcast[,c("Fam_ID","Sample_ID","Pat_ID","Mat_ID","Sex","Affection", 
                                                                                        colnames(PAWS_snps_recoded_filtered_ped_melt_dcast)[2:(ncol(PAWS_snps_recoded_filtered_ped_melt_dcast)-5)])]
save(PAWS_snps_recoded_filtered_ped_melt_dcast, file="PAWS_snps_filtered_ped_melt_dcast_recoded.RData")
write.table(as.matrix(PAWS_snps_recoded_filtered_ped_melt_dcast), "PAWS_snps_filtered_recoded_ped.txt", sep="\t", row.names = F,col.names = F,quote=FALSE)


# 
# 
# 
# ############################################################# recode the hap map map file to makes the physical cooridnates of the psychip
# hapmap_map<-read.table("hapmap3_r2_b36_fwd.consensus.qc.poly.map", header=F)
# psychchip_map<-read.table("PAWS_snps_filtered_MAP_recoded.txt", header=F)
# 
# hapmap_map[which(hapmap_map$V2%in%psychchip_map$V2[1:10]),]
# psychchip_map[1:10,]
# 
# # share an rs ID
# length(which(psychchip_map$V2%in%hapmap_map$V2))#182280/524486 (30%)
# length(which(hapmap_map$V2%in%psychchip_map$V2))#182280/1440616
# shared<-psychchip_map$V2[which(psychchip_map$V2%in%hapmap_map$V2)]
# 
# hapmap_map_onpsych<-hapmap_map[which(hapmap_map$V2%in%psychchip_map$V2),]
# hapmap_map_not_onpsych<-hapmap_map[which(!(hapmap_map$V2%in%psychchip_map$V2)),]
# 
# 
# hapmap_map_onpsych_merge<-merge(hapmap_map_onpsych, psychchip_map, by="V2")
# 
# chr_diff<-sapply(1:nrow(hapmap_map_onpsych_merge), function(x)if(hapmap_map_onpsych_merge$V3.x[x]==hapmap_map_onpsych_merge$V3.y[x]){x})
# 
# ### chromosome and physical position take from psychchip
# hapmap_map_onpsych_psych_change<-hapmap_map_onpsych_merge[,c(5,1,3,7)]
# colnames(hapmap_map_onpsych_psych_change)<-c("V1","V2","V3","V4")
# 
# hapmap_map_psypos_change<-rbind(hapmap_map_onpsych_psych_change, hapmap_map_not_onpsych)
# 
# hapmap_map_psypos_change<-hapmap_map_psypos_change[match(hapmap_map$V2, hapmap_map_psypos_change$V2),]
# 
# write.table(as.matrix(hapmap_map_psypos_change), "hapmap3_r2_b36_fwd.consensus.qc.poly.map.recoded", sep="\t", row.names = F,col.names = F,quote=FALSE) # names changed to hapmap3_r2_b36_fwd.consensus.qc.poly.map for PLINK
# 
# write.table(shared, file="list.snps", sep="\t", col.names=F, row.names=F, quote=F )
