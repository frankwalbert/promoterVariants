# this script merges the oligo based results from our experiments with the variant annotations
# note that the produt of this file (R_resultsAndAnnotations_200329.RData) will later be augmented with updated allele frequencies
# see the feature annotation folder for this
library(tidyverse)

# LOAD THE COMPLETED VARIANT ANNOTATION
load("R_variantsDF_181212.RData")
# you'll find that this has TFBS information which is not explaiend by the code in this repo
# those were earlier attempts that were then replaced by Kaushik's analyses

# oligo-based effects
load("R_designedVariantsList_200329_noControls.RData")

#for each variant, carve out the most significant oligo effect *per strand*
bestOligoEffects <- do.call(rbind.data.frame, sapply(variantsDF$variantID, function(thisVarStrand){
    thisVarSplit <- strsplit(thisVarStrand, "_")[[1]]
    thisStrand <- thisVarSplit[length(thisVarSplit)]
    thisVar <- paste(thisVarSplit[1:(length(thisVarSplit)-1)], collapse="_")
    theseOligoEffects <- designedVariantsList[[thisVar]]
    theseOligoEffects <- theseOligoEffects[theseOligoEffects[,"strand"] == thisStrand,]
    # only do if we *have* results for this variant
    # catch same variant multiple times on same strand here, with the "[1]", which should pick just one of the cases if they have the same p-value (i.e., if the two lines correspond to the same test)
    # note that this issue means that the same set of upstream oligos got desiged & printed twice. In the annotation, they cannot be disambiguated and therefore become the same result?
    if(!is.null(theseOligoEffects)){
        if(nrow(theseOligoEffects) > 0){
            return(theseOligoEffects[which(theseOligoEffects$P.Value == min(theseOligoEffects$P.Value))[1],])
        }
    }
}))
# this has 6362 rows


# fuse with annotations
resultsAndAnnotations <- merge(bestOligoEffects, variantsDF, by.x=0, by.y="variantID")
# tests & reduce redundancy
summary(resultsAndAnnotations$gene.x == resultsAndAnnotations$gene.y)
#   Mode   FALSE    TRUE
#logical      21    6341

# what are those FALSEs?
resultsAndAnnotations[resultsAndAnnotations$gene.x != resultsAndAnnotations$gene.y, 1:20]
# upstream (mostly?) variants of multiple genes where in the annotation we picked another than the one that came back most signficant in the oligo results
# we could remove gene info from those so that we're at least not confused
# better still, add a column to indicate collisions & keep them in the data
# in fact, even cleaner to go back into the annotation and add a multi-gene column there
# then, here, add a warning of gene inconsistency
# because there might be multiple genes tagged for a variant but the best oligo result happens to be the same as the gene we picked

resultsAndAnnotations$geneConflict <- resultsAndAnnotations$gene.x == resultsAndAnnotations$gene.y
# adjust names for the two columns
colnames(resultsAndAnnotations)[which(colnames(resultsAndAnnotations) == "gene.x")] <- "gene"
colnames(resultsAndAnnotations)[which(colnames(resultsAndAnnotations) == "gene.y")] <- "otherGene"

summary(resultsAndAnnotations$strand.x == resultsAndAnnotations$strand.y)
#    Mode    TRUE
# logical    6362

# all good; remove one & adjust other
resultsAndAnnotations <- resultsAndAnnotations[,-which(colnames(resultsAndAnnotations) == "strand.y")]
colnames(resultsAndAnnotations)[which(colnames(resultsAndAnnotations) == "strand.x")] <- "strand"

# any other collisions?
colnames(resultsAndAnnotations)[which(str_detect(colnames(resultsAndAnnotations), "\\.x"))]
# none for .x, and none for .y

colnames(resultsAndAnnotations)[which(colnames(resultsAndAnnotations) == "Row.names")] <- "variantID"
# make a prettier order of columns?
resultsAndAnnotations <- resultsAndAnnotations[,c(1:10, 15, 29, 371, 11:14, 16:28, 30:370)]

#save(resultsAndAnnotations, file="R_resultsAndAnnotations_200329.RData")

# confirming how many unique variants this corresponds to:
notYetUniqueVars <- sapply(resultsAndAnnotations$variantID, function(x){paste0(strsplit(x, "_")[[1]][1:4], collapse="_")})
length(unique(notYetUniqueVars))
# 5832
# this is the expected number of tested variants

# and how many significant variants?
length(which(resultsAndAnnotations$adj.P.Val < 0.05))
# 459 (more then 451 due to strandedness)
