# compute changes to TATA boxes and start codons

library(GenomicRanges)
library(tidyverse)
library(Biostrings)
library(magrittr)

# the reference genome
genome <- readDNAStringSet("sacCer3.fa")

# TATA box from Basehoar 2004
# consensus is TATA(A/T)A(A/T)(A/G)
# https://www.sciencedirect.com/science/article/pii/S0092867404002053?
# https://ccg.vital-it.ch/pwmtools/pwmscan.php

TATAs <- c(
    "TATAAAAA",
    "TATATAAA",
    "TATAAATA",
    "TATATATA",
    "TATAAAAG",
    "TATATAAG")
# not "TATATATG" and "TATAAATG" since those have a start codon

hasAConsensusSequence <- function(consensus, sequence){
    if(nchar(sequence) < nchar(consensus[1])){return(NA)}
    sum(str_detect(sequence, consensus)) > 0
}

countConsensiInSequence <- function(consensus, sequence){
    if(nchar(sequence) < nchar(consensus[1])){return(NA)}
    str_count(sequence, paste0("(?=", consensus, ")"))
}


# all designed variants are in the design:
# we can annotate all of them; no need to restrict to TSS/UpStream here
designedOligos <- read.table("jointDesignFinal27k_150612_ReverseComplementWhenALargerT.txt", head=TRUE, sep="\t", stringsAsFactors=FALSE)
variants <- unique(unlist(apply(designedOligos, 1, function(x){paste(strsplit(x["variantIDs"], ";")[[1]], x["strand"], sep="_")})))
variants <- variants[!str_detect(variants, "hunter")]
variants <- variants[!str_detect(variants, "sharon")]
# now at 10,821 variants

variantsDF <- data.frame(t(sapply(variants, function(x){strsplit(x, "_")[[1]]})), stringsAsFactors=FALSE)
colnames(variantsDF) <- c("chr", "pos", "ref", "alt", "strand")
variantsDF$SNP <- nchar(variantsDF[,3]) == 1 & nchar(variantsDF[,4]) == 1
# in cases with commas in the alt allele, we only made the first
variantsDF$alt[str_detect(variantsDF$alt, ",")] <- sapply(variantsDF$alt[str_detect(variantsDF$alt, ",")], function(x){str_split(x, ",")[[1]][1]})
variantsGRanges <- GRanges(seqnames = variantsDF$chr, ranges = IRanges(start=as.numeric(variantsDF$pos), end=as.numeric(variantsDF$pos) + nchar(variantsDF$ref) - 1), strand = variantsDF$strand, ref=variantsDF$ref, alt=variantsDF$alt)




###############################################
# pull ref & alt sequences that are 2x consensus + 1, centered on variant
# ask if either contains a consensus
# report if there is a difference in consensus

getConsensusDifferences <- function(variants, consensus){
    cl <- makeCluster(24, type="FORK")
    clusterSetRNGStream(cl)
    retList <- parLapply(cl, (1:length(variants)), function(thisVariant){
    #retList <- lapply((1:length(variants)), function(thisVariant){
        #t(sapply((1:length(variants)), function(thisVariant){
        #print(thisVariant)
        varSplit <- strsplit(variants[thisVariant], "_")[[1]]
        allele1 <- varSplit[3]
        allele2 <- varSplit[4]
        # some alleles have commas...
        if(str_detect(allele2, ",")){allele2 <- strsplit(allele2, ",")[[1]][1]}
        # extract the sequence left and right of, but excluding, the reference allele for the length of the consensus
        extractedRefSeqLeft <- as.character(subseq(genome[seqnames(variantsGRanges)[thisVariant]], start=start(variantsGRanges)[thisVariant] - nchar(consensus[1]) + 1, end=start(variantsGRanges)[thisVariant] - 1))
        extractedRefSeqRight <- as.character(subseq(genome[seqnames(variantsGRanges)[thisVariant]], start=start(variantsGRanges)[thisVariant] + nchar(allele1), end=start(variantsGRanges)[thisVariant] + nchar(allele1) + nchar(consensus[1]) - 2))
        # stick the ref and alt alleles in between
        refSeq <- paste(extractedRefSeqLeft, allele1, extractedRefSeqRight, sep="")
        altSeq <- paste(extractedRefSeqLeft, allele2, extractedRefSeqRight, sep="")
        #print(nchar(extractedRefSeq))
        #print(nchar(altSeq))
        ## note that my variant annotation is always + strand, so we can revComp the plus
        if(as.character(strand(variantsGRanges[thisVariant])) == "-"){
            refSeq <- as.character(reverseComplement(DNAString(refSeq)))
            altSeq <- as.character(reverseComplement(DNAString(altSeq)))
        }
        #print(extractedRefSeq)
        #print(altSeq)
        # report the difference in the number of consensus sequences (e.g., TATATATA can be mutated quite a bit and still be a TATA box)
        # must sum because some consensi have multiple entries (e.g. six TATA boxes)
        refConsensi <- sum(countConsensiInSequence(consensus, refSeq))
        altConsensi <- sum(countConsensiInSequence(consensus, altSeq))
        #print(c(refHasConsensus, refHasConsensus))
        ret <- c(refConsensi, altConsensi, altConsensi - refConsensi)
        names(ret) <- c("refConsensusCount", "altConsensusCount", "consenusCountDifference")
        return(ret)
    })
    stopCluster(cl)
    ret <- do.call(rbind, retList)
    rownames(ret) <- variants
    ret
}
TATA_differences <- getConsensusDifferences(variants, TATAs)
ATG_differences <- getConsensusDifferences(variants, "ATG")

length(which(TATA_differences[,3] != 0))
# 201
head(TATA_differences[which(TATA_differences[,3] != 0),])
#TATA_differences["chrII_119244_CA_C_+",]
# symmetry?
length(which(TATA_differences[,3] > 0))
# 94
length(which(TATA_differences[,3] < 0))
# 107
# fairly similar; no rampant annotation bias


length(which(ATG_differences[,3] != 0))
# 769
head(ATG_differences[which(ATG_differences[,3] != 0),])
# symmetry?
length(which(ATG_differences[,3] > 0))
# 382
length(which(ATG_differences[,3] < 0))
# 387
# looks good; does not look like there's a big reference bias

#save(TATA_differences, file="R_TATA_differences_181211.RData")
#save(ATG_differences, file="R_ATG_differences_181211.RData")
