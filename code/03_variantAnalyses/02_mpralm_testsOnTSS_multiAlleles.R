# statistical tests on the TSS library
# creates a "toptab" object that will be used downstream

# currently need to run on R 3.5.0 on MSI
library(mpra)
library(tidyverse)
library(qvalue)

# this holds all the counts from both libraries
load("../R_countsMappedToDesignedOligos_180616.RData")

# pull out the correct one
dat <- countsMappedToDesignedOligos[["TSS:scaleUpJune2016"]]
# for printing supplement:
datSummed4Print <- dat %>% group_by(oligoID) %>% summarize_at(20:ncol(dat), sum)
# write out datSummed for supplement
# in resulting file, deleted DNA counts for 2018 A&B to indictae that these were replaced
#write.table(datSummed4Print, file="datSummed_TSS_200810.txt", sep="\t", quote=FALSE, row.names=FALSE)


length(unique(countsMappedToDesignedOligos[["TSS:scaleUpJune2016"]]$oligoBlock))
# 2022
# combine technical replicates
dat <- dat %>% mutate(
    DNA_TSS_scaleUpJune2016_replicates_2016_A = DNA_TSS_scaleUpJune2016_replicates_2016_A1 + DNA_TSS_scaleUpJune2016_replicates_2016_A2,
    DNA_TSS_scaleUpJune2016_replicates_2016_B = DNA_TSS_scaleUpJune2016_replicates_2016_B1 + DNA_TSS_scaleUpJune2016_replicates_2016_B2,
    RNA_TSS_scaleUpJune2016_replicates_2016_A = RNA_TSS_scaleUpJune2016_replicates_2016_A1 + RNA_TSS_scaleUpJune2016_replicates_2016_A2,
    RNA_TSS_scaleUpJune2016_replicates_2016_B = RNA_TSS_scaleUpJune2016_replicates_2016_B1 + RNA_TSS_scaleUpJune2016_replicates_2016_B2
) %>% select(-matches(".*2016_..$") ) # dplyr::select is necessary if have MASS loaded...
dat <- dat[,c(1, 7, 8, 16, 20:ncol(dat))] # 3 would be annotation count; we don't need them here


# for final TSS replicates, need to replace the failed DNA:
dat[,"DNA_TSS_scaleUpJune2016_replicates_2018_A"] <- dat[,"DNA_TSS_scaleUpJune2016_replicates_2017_A"]
dat[,"DNA_TSS_scaleUpJune2016_replicates_2018_B"] <- dat[,"DNA_TSS_scaleUpJune2016_replicates_2017_B"]

# our comparison will be all-BY vs each of the others
# so let's condense all allele strings that are all 1s to just "1"
# to save some columns in the resulting matrix
dat$alleleString[!str_detect(dat$alleleString, "2")] <- "1"

# throw out oligos with zero counts
# first sum the barcodes
datSummed <- dat %>% group_by(oligoID) %>% summarize_at(5:ncol(dat), sum)
# which oligos have no zero counts ever, anywhere?
keepThese <- datSummed$oligoID[which(!apply(datSummed[,2:ncol(datSummed)], 1, function(x){0 %in% x}))]
length(keepThese)
# keep 4900 out of 6463. oof. lose 1/4. could revisit
dat <- dat[dat$oligoID %in% keepThese,]
# down to 5,316,181 from 5,345,528

# because of the mpralm fluke in which it apparently works with oligos that HAVE NO DATA, we need to keep only data in which all to-be-compared oligos are actually present!
presentOligos <- dat %>% group_by(oligoBlock) %>% summarize(nOligo = length(unique(oligoID)))
# here we cannot filter for "complete" blocks, nor do we want to
# because the TSS blocks have different numbers of oligos, and we can do pairwise tests for which both oligos are there even if another is missing
# but we can require the "1" oligo to be there
firstOligos <- dat %>% group_by(oligoBlock) %>% summarize(firstOligo = sort(unique(alleleString))[1])
length(which(firstOligos[,2] == "1"))
# 1495 out of 1623
# we started with 2,022 oligoBlocks in the initial counts. ~1/4 loss again
hasBYOligo <- firstOligos %>% filter(firstOligo == "1")

dat <- dat[dat$oligoBlock %in% hasBYOligo$oligoBlock,]
# down to 5,286,534

# now we melt:
# note that we want to keep the Fraser oligos so that the results can be compared across libraries
meltDat <- gather(dat, key="sample", value="count", 5:ncol(dat)) %>% mutate(molecule = str_sub(sample, 1, 3)) %>% mutate(sample = str_sub(sample, start=5)) %>% mutate(sample = paste(sample, alleleString, sep = ";allele_")) %>% filter(oligoBlock != "EranControl")

# we want a structure where every block gets split into several two way comparisons
# each comparing all_BY to one other allele
# we can call the oligoBlocks "block_otherAllele" to differentiate them, and drag along block ID as something to group on
# loop across all the other allele strings in a given oligo block
otherAlleles <- unique(meltDat$alleleString)
otherAlleles <- otherAlleles[-which(otherAlleles=="1")]

allBYCounts <- lapply(c("DNA", "RNA"), function(x){
        ret <- meltDat %>% filter(molecule==x) %>% filter(alleleString == "1")
        ret
})
names(allBYCounts) <- c("DNA", "RNA")

pairwiseCounts <- lapply(c("DNA", "RNA"), function(x){
    ret1 <- lapply(otherAlleles, function(thisAllele){
        #ret1 <- lapply("111111112", function(thisAllele){
        print(c(x, thisAllele))
        theseAlleleCounts <- meltDat %>% filter(molecule==x) %>% filter(alleleString == thisAllele)
        theseallBYCounts <- allBYCounts[[x]] %>% filter(oligoBlock %in% theseAlleleCounts$oligoBlock) # this pulls allBY alleles only from those oligo blocks that have the given alt allele
        ret <- rbind(theseAlleleCounts, theseallBYCounts) %>% spread(key=sample, value=count)
        # replace the alt allele in the column names with a common word so they can be glued together later
        colnames(ret)[str_detect(colnames(ret), paste(";allele_", thisAllele, sep=""))] <- str_replace(colnames(ret)[str_detect(colnames(ret), paste(";allele_", thisAllele, sep=""))], paste(";allele_", thisAllele, sep=""), ";otherAllele")
        # need to expand the block ID with the allele being compared to allBY here:
        ret$oligoBlock <- paste(ret$oligoBlock, thisAllele, sep="_")
        ret
    })
    names(ret1) <- otherAlleles
    ret1
})
names(pairwiseCounts) <- c("DNA", "RNA")
save(pairwiseCounts, file="R_pairwiseCounts_200326.RData")

# glue the many matrixes together into the 'dna' and 'rna' objects
DNACounts <- bind_rows(pairwiseCounts[["DNA"]])
RNACounts <- bind_rows(pairwiseCounts[["RNA"]])
dim(DNACounts)
# 6,560,693        29
# after filter: 6,471,168

dna <- DNACounts[,c(6:ncol(DNACounts))]
rna <- RNACounts[,c(6:ncol(RNACounts))]
eid <- DNACounts$oligoBlock

# note that design has to be a numeric matrix. Not sure I like that...
designFactors <- data.frame(intercept=1, allele=str_sub(str_extract(colnames(dna), ";(.)*"), 2))
#designFactors[,"year"] <- sapply(colnames(dna), function(x){strsplit(x, "_")[[1]][4]}) # year
# decided to leave out year as a covariate. results are extremely highly correlated with and without year
#cor.test(-log10(tt$P.Value), -log10(toptab[rownames(tt),]$P.Value))
# cor=0.9946429
# and for logFC: cor=0.9994852
# also now consistent with Upstream

design <- apply(designFactors, 2, function(x){as.numeric(as.factor(x))})

sampleBlocks <- str_sub(str_extract(colnames(dna), "(.)*;"), end = -2)

mpraset <- MPRASet(DNA = dna, RNA = rna, eid = eid, eseq = NULL, barcode = NULL)

pdf("mpralm_200326.pdf")
mpralm_fit <- mpralm(object = mpraset, design = design, aggregate = "sum", normalize = TRUE, model_type = "corr_groups", block=sampleBlocks, plot = TRUE)
dev.off()

toptab <- topTable(mpralm_fit, coef = 2, number = Inf)
toptab <- toptab[order(toptab$P.Value),]
#save(toptab, file="R_toptab_200326.RData")
# this will be used downstream

pdf("pVals.pdf")
hist(toptab$P.Value, breaks=40)
dev.off()

# how many variant effects in this?
# toptab contains all-allele swaps – need to remove those
# how many significant without the controls:
length(which((!str_detect(rownames(toptab), "unter")) & toptab$"adj.P.Val" < 0.1))
# 351 at FDR10, 243 at FDR5
# DO NOT USE IN PAPER, CONTAINS the all-RM SWAPS (i.e., not single variants necessarily)

# qvalue
1- qvalue(toptab$"P.Value"[!str_detect(rownames(toptab), "unter")])$pi0
# 0.3202846

# again, without the all-allele swaps
toptab1 <- toptab[!str_detect(rownames(toptab), "unter"),]
altAlleles <- str_split_fixed(rownames(toptab1), "_", n=4)[,4]
toptab1 <- toptab1[!(lapply(str_split(altAlleles, ""), function(x){length(unique(x))}) == 1 & nchar(altAlleles) > 1),]
dim(toptab1)
# 2427 <= this is the number of pairwise variant tests that were actually performed

# in how many genes?
length(unique(sapply(rownames(toptab1), function(x){str_split(x, "_")[[1]][1]})))
# 1429

# how many unique variants is this?
length(rownames(toptab1)[which(toptab1$"adj.P.Val" < 0.05)])
# 166 at FDR5, 248 at FDR10
# USE THIS IN PAPER, NOT THE NUMBERS ABOVE THAT INCLUDE THE all-RM ALLELES

# pi1?
1-qvalue(toptab1$P.Value)$pi0
# 0.2644029


###############################
# plot to make sure everything checks out now, after the filters:

#TSS
load("R_readOligos.RData")
oligoDesign <- readOligos
# Upstream
load("R_readOligos_UpStream.RData")
oligoDesign <- rbind(oligoDesign, readOligos)

# first sum all the counts per oligo
# above, we made meltDat:
meltDatForPlots <- gather(dat, key="sample", value="count", 5:ncol(dat)) %>% mutate(molecule = str_sub(sample, 1, 3)) %>% mutate(sample = str_sub(sample, start=5)) %>% filter(oligoBlock != "EranControl")

meltDatSummed <- meltDatForPlots %>% group_by(oligoID, sample, molecule) %>% summarize(count=sum(count)) %>% spread(key=molecule, value=count)
# from https://stackoverflow.com/questions/49404461/dividing-columns-by-particular-values-using-dplyr
meltDatSummedNorm <- meltDatSummed %>% group_by(sample) %>% mutate_at(vars(DNA, RNA), funs(./ sum(.))) %>% mutate(Expression = log2(RNA/DNA))


plotSingleVariantTSS <- function(thisTest){
    # need to split out the oligo block and the allele string
    thisTestSplit <- str_split(thisTest, "_")[[1]]
    thisOligoBlock <- paste0(thisTestSplit[1:(length(thisTestSplit)-1)], collapse="_")
    thisAltString <- thisTestSplit[length(thisTestSplit)]
    thisRefString <- paste0(rep("1", nchar(thisAltString)), collapse="")
    theseAlleleStrings <- c(thisRefString, thisAltString)
    theseOligoIDs <- oligoDesign[oligoDesign$oligoBlock == thisOligoBlock & oligoDesign$alleleString %in% theseAlleleStrings,]$oligoID
    theseVariants <- oligoDesign[oligoDesign$oligoBlock == thisOligoBlock, ]$variantIDs[1]
    thisVariant <- str_split(theseVariants, ";")[[1]][str_locate(thisAltString, "2")[1]]
    
    # now pull the data from the tables above
    datForPlot <- meltDatSummedNorm %>% filter(oligoID %in% theseOligoIDs)
    datForPlot$plotOligoName <- paste(datForPlot$oligoID, oligoDesign[datForPlot$oligoID,]$alleleString, sep="\nAllele: ")

    pdf(paste0(thisTest, ".pdf"), width=7, height=5, useDingbats=FALSE)
    theme_set(theme_bw())
    p <- ggplot(datForPlot, aes(x = plotOligoName, y=Expression)) +
    geom_boxplot(outlier.shape = NA) +
    #geom_jitter(position=position_jitter(width=0.2), aes(color=sample)) +
    geom_point(aes(color=sample)) +
    geom_line(aes(group=sample, colour=sample),alpha=0.3) +
    labs(title="", x = paste(thisTest, thisVariant, sep="\nVariant: "), y = "log2(RNA / DNA)") +
    scale_colour_viridis_d()
    print(p)
    dev.off()
}

for(i in 1:10){
    plotSingleVariantTSS(rownames(toptab1)[i])
}
# none of the top ten look like artifacts, missing data, etc

# SFA1:
plotSingleVariantTSS("YDL168W_5UTR_0_21")
