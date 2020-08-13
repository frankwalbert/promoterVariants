# there had been an issue with the original extraction of population frequencies from Peter et al, which resulted in too much missing data
# this gets fixed here
# see 05_populationFrequencyUpdate_200719.R for how these got entered into the data

library(tidyverse)
library(parallel)

designedOligos <- read.table("jointDesignFinal27k_150612_ReverseComplementWhenALargerT.txt", head=TRUE, sep="\t", stringsAsFactors=FALSE)
variants <- unique(unlist(apply(designedOligos, 1, function(x){paste(strsplit(x["variantIDs"], ";")[[1]], x["strand"], sep="_")})))
variants <- variants[!str_detect(variants, "hunter")]
variants <- variants[!str_detect(variants, "sharon")]
# now at 10,821 variants

variantsDF <- data.frame(t(sapply(variants, function(x){strsplit(x, "_")[[1]]})), stringsAsFactors=FALSE)
colnames(variantsDF) <- c("chr", "pos", "ref", "alt", "strand")
variantsDF$pos <- as.integer(variantsDF$pos)
# in cases with commas in the alt allele, we only made the first
variantsDF$alt[str_detect(variantsDF$alt, ",")] <- sapply(variantsDF$alt[str_detect(variantsDF$alt, ",")], function(x){str_split(x, ",")[[1]][1]})

# having this column becomes important when merging below:
variantsDF$variantID <- rownames(variantsDF)


load("R_PeterGenotypes_noGTColumns_RMcolumn_181211.RData") # stripped of individual genotypes
PeterGenotypes <- PeterGenotypesRM

PeterGenotypes[PeterGenotypes$CHROM == "chrI" & PeterGenotypes$POS == 108747,]

# DO NOT PARALLELIZE if have the big vcf loaded!
cl <- makeCluster(4, type="FORK")
clusterSetRNGStream(cl)
PeterAlleleFreqs <- parApply(cl, variantsDF, 1, function(x){
#PeterAlleleFreqs <- apply(variantsDF[variantsDF$chr == "chrI" & variantsDF$pos == 108747,], 1, function(x){
#PeterAlleleFreqs <- apply(variantsDF[variantsDF$chr == "chrII" & variantsDF$pos == 718018,], 1, function(x){
#PeterAlleleFreqs <- apply(variantsDF, 1, function(x){
    ret <- NA
    print(x)
    # the as.numeric is essential here because the pos column in variantsDF becomes a CHARACTER??? within the loop, as well as in PeterGenotypes
    thisGT <- PeterGenotypes[PeterGenotypes$CHROM == x["chr"] & as.numeric(PeterGenotypes$POS) == as.numeric(x["pos"]), ]
    print(thisGT)
    if(nrow(thisGT) == 1){
    thisINFO <- thisGT[,"INFO"]
    # which allele does RM have?
    RMAlleles <- str_split(str_split(thisGT$AAA, ":")[[1]][1], "/")[[1]]
    RMAltAllele <- as.numeric(RMAlleles[which(!RMAlleles %in% c("0"))][1])
    print(RMAlleles)
    print(RMAltAllele)
    print(thisINFO)
        thisAFField <- strsplit(as.character(thisINFO), ";")[[1]][2]
        if(str_detect(thisAFField, "AF=")){
            ret <- as.numeric(strsplit(strsplit(thisAFField, "=")[[1]][2], ",")[[1]][RMAltAllele])
        }
    }
    ret
})
stopCluster(cl)
# 2740 out of 10821 of the original variant frequencies were NA
# after this reimplementation, 202 NAs (mostly from variants not in Peter I guess)
# note that there is no check that what they call as the RM allele is the same as what we put on the array.
# maybe check if the shared values are the same as before (when we did require the alleles to be the same)

#save(PeterAlleleFreqs, file="R_PeterAlleleFreqs_200501.RData")
