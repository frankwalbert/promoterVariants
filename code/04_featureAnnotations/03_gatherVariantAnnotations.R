# assemble various annotations for each variant into one table

library(stringr)
library(data.table) # for reading the gzipped vcf and deadling with the processed & saved R object
library(R.utils) # for reading the gzipped vcf
library(parallel)
library(GenomicRanges)
library(tidyverse)
library(magrittr)



# load the design:
designedOligos <- read.table("jointDesignFinal27k_150612_ReverseComplementWhenALargerT.txt", head=TRUE, sep="\t", stringsAsFactors=FALSE)
variants <- unique(unlist(apply(designedOligos, 1, function(x){paste(strsplit(x["variantIDs"], ";")[[1]], x["strand"], sep="_")})))
variants <- variants[!str_detect(variants, "hunter")]
variants <- variants[!str_detect(variants, "sharon")]
# now at 10,821 variants

# we will build up this data frame by attaching one annotation after the other:
variantsDF <- data.frame(t(sapply(variants, function(x){strsplit(x, "_")[[1]]})), stringsAsFactors=FALSE)
colnames(variantsDF) <- c("chr", "pos", "ref", "alt", "strand")
variantsDF$pos <- as.integer(variantsDF$pos)
# in cases with commas in the alt allele, we only made the first
variantsDF$alt[str_detect(variantsDF$alt, ",")] <- sapply(variantsDF$alt[str_detect(variantsDF$alt, ",")], function(x){str_split(x, ",")[[1]][1]})

# having this column becomes important when merging below:
variantsDF$variantID <- rownames(variantsDF)
# make sure to never use the rownames below, as they will disapear or be permuted by the merges!

# SNP or not
variantsDF$SNP <- nchar(variantsDF$ref) == 1 & nchar(variantsDF$alt) == 1

# how many bp difference
variantsDF$delta_BP <- nchar(variantsDF$alt) - nchar(variantsDF$ref)

# other simple features
# allele characteristics
variantsDF$delta_A <- str_count(variantsDF$alt, "A") - str_count(variantsDF$ref, "A")
variantsDF$delta_C <- str_count(variantsDF$alt, "C") - str_count(variantsDF$ref, "C")
variantsDF$delta_G <- str_count(variantsDF$alt, "G") - str_count(variantsDF$ref, "G")
variantsDF$delta_T <- str_count(variantsDF$alt, "T") - str_count(variantsDF$ref, "T")


# TATAs and ATGs
load("R_TATA_differences_181211.RData")
variantsDF$delta_TATA <- TATA_differences[variantsDF$variantID, "consenusCountDifference"]

load("R_ATG_differences_181211.RData")
variantsDF$delta_ATG <- ATG_differences[variantsDF$variantID, "consenusCountDifference"]

# conservation of the site (frequency of allele in population from 1k genomes)
# loading this takes several minutes
# get this file from Peter et al (it is ~6GB)
# if cannot for some reason, we provide various processed versions in the repository
PeterGenotypes <- fread("1011Matrix.gvcf.gz", sep="\t", header=FALSE, stringsAsFactors=FALSE)
# a 1754867 x 1020 matrix
# fix chromsome names
PeterGenotypes[,1][PeterGenotypes[,1] == "chromosome1"] <- "chrI"
PeterGenotypes[,1][PeterGenotypes[,1] == "chromosome2"] <- "chrII"
PeterGenotypes[,1][PeterGenotypes[,1] == "chromosome3"] <- "chrIII"
PeterGenotypes[,1][PeterGenotypes[,1] == "chromosome4"] <- "chrIV"
PeterGenotypes[,1][PeterGenotypes[,1] == "chromosome5"] <- "chrV"
PeterGenotypes[,1][PeterGenotypes[,1] == "chromosome6"] <- "chrVI"
PeterGenotypes[,1][PeterGenotypes[,1] == "chromosome7"] <- "chrVII"
PeterGenotypes[,1][PeterGenotypes[,1] == "chromosome8"] <- "chrVIII"
PeterGenotypes[,1][PeterGenotypes[,1] == "chromosome9"] <- "chrIX"
PeterGenotypes[,1][PeterGenotypes[,1] == "chromosome10"] <- "chrX"
PeterGenotypes[,1][PeterGenotypes[,1] == "chromosome11"] <- "chrXI"
PeterGenotypes[,1][PeterGenotypes[,1] == "chromosome12"] <- "chrXII"
PeterGenotypes[,1][PeterGenotypes[,1] == "chromosome13"] <- "chrXIII"
PeterGenotypes[,1][PeterGenotypes[,1] == "chromosome14"] <- "chrXIV"
PeterGenotypes[,1][PeterGenotypes[,1] == "chromosome15"] <- "chrXV"
PeterGenotypes[,1][PeterGenotypes[,1] == "chromosome16"] <- "chrXVI"

# chop off the individual genotypes
PeterGenotypes <- PeterGenotypes[,1:9]
# fix column names, which we could have loaded given it DID apparently load that #CHROM entry rather than comment it out, but whatever
colnames(PeterGenotypes) <- as.character(PeterGenotypes[1,])
PeterGenotypes <- PeterGenotypes[2:nrow(PeterGenotypes),]
colnames(PeterGenotypes)[1] <- "CHROM"

# save this for later, faster loading
#save(PeterGenotypes, file="R_PeterGenotypes_181211.RData")
load("../1011_yeast_variants/R_PeterGenotypes_181211.RData")
# this is 99.9Mb compared to the 5.8Gb vcf.gz

# alternatively, a version that retains the RM genotypes so we can look up what they are
PeterGenotypesRM <- PeterGenotypes[,1:10]
save(PeterGenotypesRM, file="PeterGenotypesRM_200501.RData")

# for every variant in the table, find its frequency
# this requires a single match and only one alternative allele in the Peter data
# one could figure out how to deal with variants with multiple alleles (i.e., what is their frequency?), but skip for now
# DO NOT PARALLELIZE if have the big vcf loaded
cl <- makeCluster(24, type="FORK")
clusterSetRNGStream(cl)
PeterAlleleFreqs <- parApply(cl, variantsDF, 1, function(x){
#PeterAlleleFreqs <- apply(variantsDF, 1, function(x){
    ret <- NA
    #print(x)
    # the as.numeric is essential here because the pos column in variantsDF becomes a CHARACTER??? within the loop, as well as in PeterGenotypes
    thisINFO <- PeterGenotypes[PeterGenotypes$CHROM == x["chr"] & as.numeric(PeterGenotypes$POS) == as.numeric(x["pos"]) & PeterGenotypes$REF == x["ref"] & PeterGenotypes$ALT == x["alt"], "INFO"]
    #print(thisINFO)
    if(nrow(thisINFO) == 1){
        thisAFField <- strsplit(as.character(thisINFO), ";")[[1]][2]
        if(str_detect(thisAFField, "AF=")){
            ret <- as.numeric(strsplit(thisAFField, "=")[[1]][2])
        }
    }
    ret
})
stopCluster(cl)
# 2740 out of 10821 are NA
#save(PeterAlleleFreqs, file="../1011_yeast_variants/R_PeterAlleleFreqs_181211.RData")

# IMPORTANT: we will later replace this with a better frequency estimate
# in file 04_populationFrequencyUpdate_200501.R
# and add minor (and derived) allele frequency
# in file 05_populationFrequencyUpdate_200719.R
variantsDF$populationFreq <- PeterAlleleFreqs[variantsDF$variantID]


# essential gene etc; from our Albert/Bloom 2018 eLife paper
# FIRST ASSIGN EACH VARIANT TO ITS GENE
# use designedOligos for this
variantsToGenes <- do.call(rbind, lapply((1:nrow(designedOligos)), function(x){
    theseVars <- paste(strsplit(designedOligos[x, "variantIDs"], ";")[[1]], designedOligos[x, "strand"], sep="_")
    thisGene <- strsplit(designedOligos[x, "gene"], "_")[[1]][1]
    cbind(theseVars, thisGene)
}))
variantsToUniqueGenes <- sapply(variantsDF$variantID, function(x){unique(variantsToGenes[variantsToGenes[,1] == x,2])[1]})
# when run without that [1] in the line just above, this becomes a list
#length(which(sapply(variantsToUniqueGenes, length)>1))
# 42 variants were designed for more than one gene
# for these, pick the "first"; see code above

# and record how many genes a variant tags
numberOfGenesPerVariant <- sapply(variantsDF$variantID, function(x){length(unique(variantsToGenes[variantsToGenes[,1] == x,2]))})

variantsDF$gene <- variantsToUniqueGenes[variantsDF$variantID]
variantsDF$numberOfGenes <- numberOfGenesPerVariant[variantsDF$variantID]


# how many genes per variant?


# this is the same table used in Albert/Bloom 2018
load("R_externalGeneData_170501.RData")
variantsDF <- merge(variantsDF, externalGeneData, by.x="gene", by.y=0, all.x=TRUE, sort=FALSE)
# THIS MERGE CHANGES THE ORDER OF THE ROWS AND REMOVED THE ROWNAMES!


# distance to the associated gene
geneAnnotation = read.table("ensemblGenes_ensembl83_160307_MOD.txt", sep="\t", stringsAsFactors=FALSE, head=TRUE)
rownames(geneAnnotation) = geneAnnotation[,1]

distFromGene <- sapply(variantsDF$variantID, function(y){
    ret <- NA
    x <- variantsDF[variantsDF$variantID == y,]
    if(x$strand == "+"){
        ret <- as.numeric(geneAnnotation[x$gene, "start"]) - as.numeric(x$pos)
    }
    if(x$strand == "-"){
        ret <- as.numeric(x$pos) - as.numeric(geneAnnotation[x$gene, "end"])
    }
    ret
})
variantsDF$distFromGene <- distFromGene[variantsDF$variantID]


# in a nucleosome?
nucleosomes <- read.table("nucleosome_data_SacCer3.bed", stringsAsFactors=FALSE, sep="\t", header=TRUE)
# the nucleosomes may overlap! find the best score per variant
nucleosomeGRanges <- GRanges(seqnames = nucleosomes$chromosome, ranges = IRanges(start=nucleosomes$start, end=nucleosomes$end), score=nucleosomes$score)

variantsGRanges <- GRanges(seqnames = variantsDF$chr, ranges = IRanges(start=as.numeric(variantsDF$pos), end=as.numeric(variantsDF$pos) + nchar(variantsDF$ref) - 1), strand = variantsDF$strand, ref=variantsDF$ref, alt=variantsDF$alt, id=variantsDF$variantID)
variantNucOverlaps <- findOverlaps(variantsGRanges, nucleosomeGRanges)

# this has redundant lookups of variantID to make sure everything lines up
nucleosomeScorePerVariant <- sapply(variantsDF$variantID, function(x){
    positionInDF <- which(variantsDF$variantID == x)
    max(nucleosomeGRanges[subjectHits(variantNucOverlaps[queryHits(variantNucOverlaps) == positionInDF,])]$score)[1]
})
#nucleosomeScorePerVariant["chrI_106011_G_GT_-"]
# what should we do with the variants that are not covered by a nucleosome?
summary(nucleosomeScorePerVariant[!is.infinite(nucleosomeScorePerVariant)])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.220   1.690   2.590   3.161   4.010  29.160
# could set to zero; does not seem like it does it justice?
# set to NA, then have separate column that just says "is this in a nucleosome at all?"

nucleosomeScorePerVariant[is.infinite(nucleosomeScorePerVariant)] <- NA
#save(nucleosomeScorePerVariant, file="../nucleosomes/R_nucleosomeScorePerVariant_181212.RData")

variantsDF$nucleosomeScore <- nucleosomeScorePerVariant[variantsDF$variantID]
variantsDF$inANucleosome <- !is.na(variantsDF$nucleosomeScore)

# SAVE THE FINAL PRODUCT
#save(variantsDF, file="R_variantsDF_181212.RData")
