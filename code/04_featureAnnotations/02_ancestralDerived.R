# assign at each variant whether RM or BY is derived

library(tidyverse)
library(parallel)


load("PeterGenotypesRMAlt_200505.RData") # entire file

# let's use strain AMH ("EN14S01") in Peter 2018 as the Taiwanese outgroup

# returns:
# TRUE if the Taiwanese strain has the reference allele
# FALSE if Taiwan has the RM allele
# NA if Taiwan has neither the ref nor the RM allele
cl <- makeCluster(24, type="FORK")
clusterSetRNGStream(cl)
RMDerived <- parApply(cl, PeterGenotypesRMAlt, 1, function(x){
#RMDerived <- apply(PeterGenotypesRMAlt[1:100,], 1, function(x){
    ret <- NA
    # the as.numeric is essential here because the pos column in variantsDF becomes a CHARACTER??? within the loop, as well as in PeterGenotypes
    thisGT <- data.frame(t(x))
    #print(thisGT[1:3])
    if(nrow(thisGT) == 1){
    thisINFO <- thisGT$INFO
    # which allele does RM have?
    RMAlleles <- str_split(str_split(thisGT$AAA, ":")[[1]][1], "/")[[1]]
    RMAltAllele <- as.numeric(RMAlleles[which(!RMAlleles %in% c("0"))][1])
    print(RMAlleles)
    #print(RMAltAllele)
    
    # which allele does the Taiwanese strain have?
    TaiwanAlleles <- str_split(str_split(thisGT$AMH, ":")[[1]][1], "/")[[1]]
    print(TaiwanAlleles)
    # if Taiwan or RM is heterozygous, don't make a call
    if(length(unique(TaiwanAlleles)) > 1){return(NA)}
    if(length(unique(RMAlleles)) > 1){return(NA)}
    
    # if Taiwan is homozygous REF, RM is derived
    if(unique(TaiwanAlleles) == "0"){return(TRUE)}
    
    # if Taiwan has the RM allele, BY is derived (given we ruled out a 0/1 genotype above)
    if(RMAltAllele %in% as.numeric(TaiwanAlleles)){return(FALSE)}
    else{return(NA)}
    }
})
#RMDerived
stopCluster(cl)

summary(RMDerived)
#   Mode   FALSE    TRUE    NA's
#logical   21195   22827    4307
# this looks pretty good?

# attach variant info:
RMDerived_withInfo <- data.frame(PeterGenotypesRMAlt[,1:5], RMDerived)

#save(RMDerived_withInfo, file="R_RMDerived_withInfo_200718.RData")
