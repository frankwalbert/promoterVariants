# we had earlier made improved allele frequencies and ancestral derived assignments
# these were not yet correctly integrated with the variant effect table
# do so here


library(tidyverse)

load("R_resultsAndAnnotations_200329.RData")

# first, remove my TFBSs – Kaushik redid all of them
resultsAndAnnotations <- resultsAndAnnotations[,1:44]

# the correct (!) population frequency
load("R_PeterAlleleFreqs_200501.RData")

resultsAndAnnotations$populationFreq <- PeterAlleleFreqs[resultsAndAnnotations$variantID]

# add MAF
MAF <- resultsAndAnnotations$populationFreq
MAF[MAF > 0.5 & !is.na(MAF)] <- 1- MAF[MAF > 0.5 & !is.na(MAF)]
resultsAndAnnotations$MAF <- MAF

# DAF
load("R_RMDerived_withInfo_200718.RData")
DAFTable <- merge(resultsAndAnnotations, RMDerived_withInfo, by.x=c("chr", "pos", "ref"), by.y=c("CHROM", "POS", "REF"), all.x=TRUE)
DAF <- DAFTable$populationFreq
DAF[(!DAFTable$RMDerived) & (!is.na(DAFTable$RMDerived))] <- 1 - DAF[(!DAFTable$RMDerived) & (!is.na(DAFTable$RMDerived))]
DAF[is.na(DAFTable$RMDerived)] <- NA
DAFTable$DAF <- DAF
# bring back into same order of columns (rows are in different order due to the merge operation)
DAFTable <- DAFTable[,c(colnames(resultsAndAnnotations), "DAF", "RMDerived")]

resultsAndAnnotations <- DAFTable


#save(resultsAndAnnotations, file="R_resultsAndAnnotations_200719.RData")
# this is the final file that then got used by Kaushik for his downstream analyses
