#!/usr/bin/Rscript
#read in vars

setwd("~/bootcamp/PsychChip/impute")
DATA <- "PAWS_snps_filtered_recoded"

bound <- read.table("~/bootcamp/PsychChip/impute/chromosome_boundaries.txt", h = F)

for(i in 1:20){
  
  chr <- i
  
  n_int <- round(bound[i,3]/5000000)
  
  for(j in 1:(n_int-1)){
    
    l_range <- j*5000000
    u_range <- (j+1)*5000000
    
    option <- ""
    
    if(j == 1){
      
      l_range <- 1
      u_range <- (j+1)*5000000
      
      option <- "-allow_large_regions"
      
    }
    
    if(j == (n_int-1)){
      l_range <- j*5000000
      u_range <- bound[i,3]
      
      option <- "-allow_large_regions"
    }
    #command <- paste0("qsub -a 1326 -N ", DATA, "_", chr, "_impute_", l_range, "_", u_range, " app/impute.sh ", chr, " ", l_range, " ", u_range, " ", option)
    command <- paste0("./impute.sh ", chr, " ", l_range, " ", u_range, " ", option)
    print(command)
    system(command)
    
  }
  
}
