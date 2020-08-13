# process all replicates, from reads to barcode counts
# this process can take a lot of memory and time. Don't try on your laptop
# we ran this on a machine with 256Gb RAM
# this pulls in the three files created from each annotation run
# the excel file that gets pulled in early has their paths – adjust accordingly

library(readxl)
library(ShortRead)
library(Biostrings)
library(reshape2)
library(plyr)
library(dplyr)

# read the sampleOverview and annotation files
# these have paths to the various raw etc files
samples <- read_excel("sampleOverview.xlsx")
annotationRuns <- read_excel("annotationRunOverview.xlsx")

# which samples to run?
# exclude samples in which the RNA failed
samples <- samples[samples$RNA_OK,]

# make a list of the annotation files referred to by the samples
annotationFilesList <- lapply(unique(paste(samples$library, ":", samples$annotationVersion, sep="")), function(i){list(library=NA, version=NA, bcCounts=NA, bcCountsAssigned=NA, readOligos=NA)})
names(annotationFilesList) <- unique(paste(samples$library, ":", samples$annotationVersion, sep=""))

# and load the actual annotation file into that list
# this will be a huge object!
for(i in 1:length(annotationFilesList)){
    print(i)
    #print(file.exists(annotationRuns[annotationRuns$annotation_ID == names(annotationFilesList)[i],]$path_bcCounts))
    #print(file.exists(annotationRuns[annotationRuns$annotation_ID == names(annotationFilesList)[i],]$path_bcCountsAssigned))
    #print(file.exists(annotationRuns[annotationRuns$annotation_ID == names(annotationFilesList)[i],]$path_readOligos))
    load(annotationRuns[annotationRuns$annotation_ID == names(annotationFilesList)[i],]$path_bcCounts)
    load(annotationRuns[annotationRuns$annotation_ID == names(annotationFilesList)[i],]$path_bcCountsAssigned)
    load(annotationRuns[annotationRuns$annotation_ID == names(annotationFilesList)[i],]$path_readOligos)
    annotationFilesList[[i]][["bcCounts"]] <- bcCountsTopOligoOnly
    annotationFilesList[[i]][["bcCountsAssigned"]] <- bcCountsAssigned
    annotationFilesList[[i]][["readOligos"]] <- readOligos
    annotationFilesList[[i]][["version"]] <- annotationRuns[annotationRuns$annotation_ID == names(annotationFilesList)[i],]$annotationVersion
    annotationFilesList[[i]][["library"]] <- annotationRuns[annotationRuns$annotation_ID == names(annotationFilesList)[i],]$library
}


# generate tables for the barcodes from DNA and RNA
# CAREFUL: this will get huge
# after reading the reads and the annotations above, the R session is at ≤50Gb
barcodes = list()
for (i in 1:nrow(samples)){
    print(i)
    #print(paste(samples[i, "readPath"], "/", "RNA_", samples[i, "sampleNumberWithinBatch"], ".fq.gz", sep=""))
    #print(file.exists(paste(samples[i, "readPath"], "/", "RNA_", samples[i, "sampleNumberWithinBatch"], ".fq.gz", sep="")))
    DNA = as.character(subseq(sread(readFastq(paste(samples[i, "readPath"], sep=""), pattern=paste("DNA_", samples[i, "sampleNumberWithinBatch"], ".fq.gz", sep=""))), start=1, width=20))
    RNA = as.character(subseq(sread(readFastq(paste(samples[i, "readPath"], sep=""), pattern=paste("RNA_", samples[i, "sampleNumberWithinBatch"], ".fq.gz", sep=""))), start=1, width=20))
    barcodes[[i]] = list(t(t(DNA)))
    barcodes[[i]][[2]] = t(t(RNA))
}
names(barcodes) <- samples$experiment_ID

# how many reads?

t(round(sapply(barcodes, function(x){sapply(x, nrow)})/1e6, 2))
#UpStream_November2016_upstream_2018_1            20.70  18.10
#UpStream_November2016_upstream_2018_2             0.00  21.16
#UpStream_November2016_upstream_2018_3            15.90  18.80
#UpStream_November2016_upstream_2018_4             0.07  19.17 <- be sure to null this guy's DNA, there are some scattered matches
#UpStream_November2016_upstream_2018_A            15.00  21.24
#UpStream_November2016_upstream_2018_B            13.56  22.56
#TSS_scaleUpJune2016_replicates_2018_A             0.00  67.99
#TSS_scaleUpJune2016_replicates_2018_B             0.00  58.03
#TSS_scaleUpJune2016_replicates_2017_A            27.06  35.43
#TSS_scaleUpJune2016_replicates_2017_B            30.45  21.60
#TSS_scaleUpJune2016_replicates_2016_1_all        49.79  55.21
#TSS_scaleUpJune2016_replicates_2016_2_all        84.49  51.33
#TSS_scaleUpJune2016_replicates_2016_3_all        77.67  48.58
#TSS_scaleUpJune2016_replicates_2016_4_all        60.62  43.76
#TSS_scaleUpJune2016_replicates_2016_5_all        65.65  62.29
#TSS_scaleUpJune2016_replicates_2016_6_all        58.43  18.94
#TSS_scaleUpJune2016_replicates_2016_A1           30.13  62.71
#TSS_scaleUpJune2016_replicates_2016_A2           33.07  64.53
#TSS_scaleUpJune2016_replicates_2016_B1           28.29  56.71
#TSS_scaleUpJune2016_replicates_2016_B2           30.98  56.26

# count barcodes in each
# (this should really be mclapply, but we'll run it only once. hopefully.)
barcodeCounts = lapply(barcodes, function(x){
    ret = lapply(x, table)
    names(ret) = c("DNA", "RNA")
    ret
})

# split the counts into the different annotation groups
# and process each annotation group separately

# put them into one table
bcCountTables <- lapply(unique(paste(samples$library, ":", samples$annotationVersion, sep="")), function(x){
    print(x)
    thisAnno <- strsplit(x, ":")[[1]][2]
    thisLibrary <- strsplit(x, ":")[[1]][1]
    thisWhich <- which(samples$library == thisLibrary & samples$annotationVersion == thisAnno)
    molten = melt(barcodeCounts[thisWhich])
    bcCountTable = dcast(molten, Var1 ~ ..., fill=0)
    names(bcCountTable)[1] <- "barcodeSeq"
    rownames(bcCountTable) = bcCountTable[,1]
    bcCountTable = bcCountTable[,2:ncol(bcCountTable)]
    bcCountTable
})
names(bcCountTables) <- unique(paste(samples$library, ":", samples$annotationVersion, sep=""))

#save(bcCountTables, file="R_bcCountTables_180616.RData")

# map to annotated barcodes (all, not yet restricted to designed oligos)
countsMappedToAnnotatedBarcodes <- lapply(names(bcCountTables), function(x){
    hiseqCountsAssignedToBarcodes <- merge(annotationFilesList[[x]][["bcCounts"]], cbind(rownames(bcCountTables[[x]]), bcCountTables[[x]]), by.x = "bcSeq", by.y = 1)
    rownames(hiseqCountsAssignedToBarcodes) <- hiseqCountsAssignedToBarcodes[,"bcSeq"]
    hiseqCountsAssignedToBarcodes
})
names(countsMappedToAnnotatedBarcodes) <- names(bcCountTables)

#save(countsMappedToAnnotatedBarcodes, file="R_countsMappedToAnnotatedBarcodes_180616.RData")


# map to designed oligos
# this is a new merge, it is not building on top of countsMappedToAnnotatedBarcodes
countsMappedToDesignedOligos <- lapply(names(bcCountTables), function(x){
    hiseqCountsAssignedToDesignedOligos <- merge(annotationFilesList[[x]][["bcCountsAssigned"]], cbind(rownames(bcCountTables[[x]]), bcCountTables[[x]]), by.x = "bcSeq", by.y = 1)
    #rownames(hiseqCountsAssignedToDesignedOligos) <- hiseqCountsAssignedToDesignedOligos[,"bcSeq"]
    hiseqCountsAssignedToDesignedOligos
})
names(countsMappedToDesignedOligos) <- names(bcCountTables)

#save(countsMappedToDesignedOligos, file="R_countsMappedToDesignedOligos_180616.RData")

# write for GEO
# all barcodes
countsUpstreamPrint <- countsMappedToAnnotatedBarcodes[["UpStream:November2016"]][,c(1, 7, 9, 11:18)]
colnames(countsUpstreamPrint) <- c("barcodeSequence", "Upstream_1_DNA", "Upstream_3_DNA", "Upstream_A_DNA", "Upstream_B_DNA", "Upstream_1_RNA", "Upstream_2_RNA", "Upstream_3_RNA", "Upstream_4_RNA", "Upstream_A_RNA", "Upstream_B_RNA")
#write.table(countsUpstreamPrint, file="Upstream_countedMappedBarcodes_200806.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

countsTSSPrint <- countsMappedToAnnotatedBarcodes[["TSS:scaleUpJune2016"]][,c(1, 7:18, 21:34)]
colnames(countsTSSPrint) <- c("barcodeSequence", "TSS_2016_1_DNA", "TSS_2016_2_DNA", "TSS_2016_3_DNA", "TSS_2016_4_DNA", "TSS_2016_5_DNA", "TSS_2016_6_DNA", "TSS_2016_A1_DNA", "TSS_2016_A2_DNA", "TSS_2016_B1_DNA", "TSS_2016_B2_DNA", "TSS_2017_A_DNA", "TSS_2017_B_DNA", "TSS_2016_1_RNA", "TSS_2016_2_RNA", "TSS_2016_3_RNA", "TSS_2016_4_RNA", "TSS_2016_5_RNA", "TSS_2016_6_RNA", "TSS_2016_A1_RNA", "TSS_2016_A2_RNA", "TSS_2016_B1_RNA", "TSS_2016_B2_RNA", "TSS_2017_A_RNA", "TSS_2017_B_RNA", "TSS_2018_A_RNA", "TSS_2018_B_RNA")
#write.table(countsTSSPrint, file="TSS_countedMappedBarcodes_200806.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

# barcodes for designed oligos only (smaller, more useful, less info)
countsUpstreamDesignPrint <- countsMappedToDesignedOligos[["UpStream:November2016"]][,c(1, 7, 20, 22, 24:31)]
colnames(countsUpstreamDesignPrint) <- c("barcodeSequence", "oligoID", "Upstream_1_DNA", "Upstream_3_DNA", "Upstream_A_DNA", "Upstream_B_DNA", "Upstream_1_RNA", "Upstream_2_RNA", "Upstream_3_RNA", "Upstream_4_RNA", "Upstream_A_RNA", "Upstream_B_RNA")
#write.table(countsUpstreamDesignPrint, file="Upstream_countedBarcodesMappedToDesign_200806.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

countsTSSDesignPrint <- countsMappedToDesignedOligos[["TSS:scaleUpJune2016"]][,c(1, 7, 20:31, 34:47)]
colnames(countsTSSDesignPrint) <- c("barcodeSequence", "oligoID", "TSS_2016_1_DNA", "TSS_2016_2_DNA", "TSS_2016_3_DNA", "TSS_2016_4_DNA", "TSS_2016_5_DNA", "TSS_2016_6_DNA", "TSS_2016_A1_DNA", "TSS_2016_A2_DNA", "TSS_2016_B1_DNA", "TSS_2016_B2_DNA", "TSS_2017_A_DNA", "TSS_2017_B_DNA", "TSS_2016_1_RNA", "TSS_2016_2_RNA", "TSS_2016_3_RNA", "TSS_2016_4_RNA", "TSS_2016_5_RNA", "TSS_2016_6_RNA", "TSS_2016_A1_RNA", "TSS_2016_A2_RNA", "TSS_2016_B1_RNA", "TSS_2016_B2_RNA", "TSS_2017_A_RNA", "TSS_2017_B_RNA", "TSS_2018_A_RNA", "TSS_2018_B_RNA")
#write.table(countsTSSDesignPrint, file="TSS_countedBarcodesMappedToDesign_200806.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


###################
#diagnose all samples

# readOligos to get distributions of barcodes per designed oligos

# reads
# barcodes
# reads that do not tag an annotated barcode (and are therefore useless; this could be improved by allowing in off-by-ones but is usually quite low anyway)
# barcodes that map to annotation
# barcodes that tag a designed oligo
# designed oligos
# designed oligos found in this sample

sampleStats <- sapply(samples$experiment_ID, function(x){
    # which batch is the sample in?
    print(x)
    thisBatch <- paste(filter(samples, experiment_ID==x)$library, filter(samples, experiment_ID==x)$annotationVersion, sep=":")
    print(thisBatch)
    ret <- sapply(c("DNA", "RNA"), function(i){
        thisColumn <- paste(i, x, sep="_")
        print(thisColumn)
        reads <- sum(bcCountTables[[thisBatch]][,thisColumn])
        barcodes <- length(which(bcCountTables[[thisBatch]][,thisColumn] > 0))
        readsWithoutBarcodes <- reads - sum(countsMappedToAnnotatedBarcodes[[thisBatch]][,thisColumn])
        annotatedBarcodes <- length(which(countsMappedToAnnotatedBarcodes[[thisBatch]][,thisColumn] > 0))
        barcodesWithDesignedOligos <- length(which(countsMappedToDesignedOligos[[thisBatch]][,thisColumn] > 0))
        designedOligosInSample <- length(unique(countsMappedToDesignedOligos[[thisBatch]][which(countsMappedToDesignedOligos[[thisBatch]][,thisColumn] > 0), "oligoID"]))
        designedOligosTotal <- nrow(annotationFilesList[[thisBatch]][["readOligos"]])
        c(reads, barcodes, readsWithoutBarcodes, annotatedBarcodes, barcodesWithDesignedOligos, designedOligosInSample, designedOligosTotal)
    })
    names(ret) <- paste(c(rep("DNA", 7), rep("RNA", 7)), c("reads", "barcodes", "readsWithoutBarcodes", "annotatedBarcodes", "barcodesWithDesignedOligos", "designedOligosInSample", "designedOligosTotal"), sep="_")
    ret
})
#write.table(sampleStats, file="sampleStats_180620.txt", sep="\t", quote=FALSE)


#################
# diagnose each batch

# oligos designed
# oligos represented (at all, by > 5, >10, >50, >100, >500, >1000 barcodes)
# oligos represented (at all, by > 5, >10, >50, >100, >500, >1000 counts after summing barcodes)
# barcodes annotated to designed oligos
# barcodes recovered in experiments
# median (etc) counts (after summing barcodes) per oligo in DNA and RNA
# barcodes recovered that tag an error oligo (how much potential is there for an analysis of these)

batchStats <- sapply(names(annotationFilesList), function(x){
    designedOligos <- filter(annotationRuns, annotation_ID==x)$numberOligosDesigned
    annotatedOligos <- filter(annotationRuns, annotation_ID==x)$numberOligosCovered
    recoveredOligos <- length(unique(countsMappedToDesignedOligos[[x]][, "oligoID"]))
    annotatedDesignBarcodes <- filter(annotationRuns, annotation_ID==x)$numberBarcodesAssignedToDesign
    recoveredDesignBarcodes <- length(unique(countsMappedToDesignedOligos[[x]][, "bcSeq"]))
    recoveredErrorBarcodes <- length(unique(countsMappedToAnnotatedBarcodes[[x]][, "bcSeq"])) - recoveredDesignBarcodes
    ret <- c(designedOligos, annotatedOligos, recoveredOligos, annotatedDesignBarcodes, recoveredDesignBarcodes, recoveredErrorBarcodes)
    names(ret) <- c("designedOligos", "annotatedOligos", "recoveredOligos", "annotatedDesignBarcodes", "recoveredDesignBarcodes", "recoveredErrorBarcodes")
    ret
})


