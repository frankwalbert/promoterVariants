# currently need to run on R 3.5.0 on MSI
library(mpra)
library(tidyverse)
library(qvalue)
library(data.table) # for reading the big gvcf <= may need to update this or can get errors during reading of file


load("R_countsMappedToDesignedOligos_180616.RData")

# combine technical replicates
dat <- countsMappedToDesignedOligos[["TSS:scaleUpJune2016"]]
dat <- dat %>% mutate(
    DNA_TSS_scaleUpJune2016_replicates_2016_A = DNA_TSS_scaleUpJune2016_replicates_2016_A1 + DNA_TSS_scaleUpJune2016_replicates_2016_A2,
    DNA_TSS_scaleUpJune2016_replicates_2016_B = DNA_TSS_scaleUpJune2016_replicates_2016_B1 + DNA_TSS_scaleUpJune2016_replicates_2016_B2,
    RNA_TSS_scaleUpJune2016_replicates_2016_A = RNA_TSS_scaleUpJune2016_replicates_2016_A1 + RNA_TSS_scaleUpJune2016_replicates_2016_A2,
    RNA_TSS_scaleUpJune2016_replicates_2016_B = RNA_TSS_scaleUpJune2016_replicates_2016_B1 + RNA_TSS_scaleUpJune2016_replicates_2016_B2
) %>% select(-matches(".*2016_..$") )
dat <- dat[,c(1, 7, 8, 16, 20:ncol(dat))] # 3 would be annotation count; we don't need them here

# for final TSS replicates, need to replace the failed DNA:
dat[,"DNA_TSS_scaleUpJune2016_replicates_2018_A"] <- dat[,"DNA_TSS_scaleUpJune2016_replicates_2017_A"]
dat[,"DNA_TSS_scaleUpJune2016_replicates_2018_B"] <- dat[,"DNA_TSS_scaleUpJune2016_replicates_2017_B"]

# our comparison will only be among 11, 12, 21, 22
dat <- dat[dat$alleleString %in% c("11", "12", "21", "22"), ]
# down to 1,624,588 from 5,345,528

# throw out oligos with zero counts
# first sum the barcodes
datSummed <- dat %>% group_by(oligoID) %>% summarize_at(5:ncol(dat), sum)
# which oligos have no zero counts ever, anywhere?
keepThese <- datSummed$oligoID[which(!apply(datSummed[,2:ncol(datSummed)], 1, function(x){0 %in% x}))]
length(keepThese)
# keep 1538 out of 2047
dat <- dat[dat$oligoID %in% keepThese,]
# down to 1,604,700

# need to keep only data in which all four oligos are actually present!
presentOligos <- dat %>% group_by(oligoBlock) %>% summarize(nOligo = length(unique(oligoID)))
completeOligoBlocks <- presentOligos %>% filter(nOligo == 4)
# 342 such blocks (was 483 when not throwing out zero counts)
dat <- dat[dat$oligoBlock %in% completeOligoBlocks$oligoBlock,]
# 1,574,894 rows


# now we melt:
meltDat <- gather(dat, key="sample", value="count", 5:ncol(dat)) %>% mutate(molecule = str_sub(sample, 1, 3)) %>% mutate(sample = str_sub(sample, start=5)) %>% mutate(sample = paste(sample, alleleString, sep = ";allele_")) %>% filter(oligoBlock != "EranControl")


# we want a structure where every block gets split into one model with four oligos
# pull out the data into DNA and RNA matrices, one for each of the four two-variant alleles
epistasisCounts <- lapply(c("DNA", "RNA"), function(x){
    theseAlleleCounts <- meltDat %>% filter(molecule==x) # %>% filter(alleleString == thisAllele)
        ret <- theseAlleleCounts %>% spread(key=sample, value=count)
        # replace the alt allele in the column names with a common word so they can be glued together later
        #colnames(ret)[str_detect(colnames(ret), paste(";allele_", thisAllele, sep=""))] <- str_replace(colnames(ret)[str_detect(colnames(ret), paste(";allele_", thisAllele, sep=""))], paste(";allele_", thisAllele, sep=""), ";otherAllele")
        # need to expand the block ID with the allele being compared to allBY here:
        #ret$oligoBlock <- paste(ret$oligoBlock, thisAllele, sep="_")
        ret
})
names(epistasisCounts) <- c("DNA", "RNA")


dna <- epistasisCounts[["DNA"]][,c(6:ncol(epistasisCounts[["DNA"]]))]
rna <- epistasisCounts[["RNA"]][,c(6:ncol(epistasisCounts[["RNA"]]))]
eid <- epistasisCounts[["DNA"]]$oligoBlock
length(unique(eid))
# 537 tests; NO! now 483 when only considering complete blocks
# and 342 when throwing out any oligos with one zero count anywhere

# note that design has to be a numeric matrix
designFactors <- data.frame(intercept=1, allele1=substr(str_sub(str_extract(colnames(dna), ";(.)*"), 2), 8, 8), allele2=substr(str_sub(str_extract(colnames(dna), ";(.)*"), 2), 9, 9))

design <- apply(designFactors, 2, function(x){as.numeric(as.factor(x))})
design <- cbind(design, model.matrix(~designFactors[,2] * designFactors[,3])[,4])
colnames(design)[ncol(design)] <- "interaction"
# model.matrix(~designFactors[,2] * designFactors[,3] + designFactors[,4])

sampleBlocks <- str_sub(str_extract(colnames(dna), "(.)*;"), end = -2)

mpraset <- MPRASet(DNA = dna, RNA = rna, eid = eid, eseq = NULL, barcode = NULL)

pdf("mpralm_200329_withHunter.pdf")
mpralm_fit <- mpralm(object = mpraset, design = design, aggregate = "sum", normalize = TRUE, model_type = "corr_groups", block=sampleBlocks, plot = TRUE)
dev.off()

topEpistasis <- topTable(mpralm_fit, coef = 4, number = Inf)
# sort by p-value
topEpistasis <- topEpistasis[order(topEpistasis$P.Value),]
#save(topEpistasis, file="R_topEpistasis_200326.RData")

pdf("pVals.pdf")
hist(topEpistasis$P.Value, breaks=40)
dev.off()

head(topEpistasis)
#                    logFC     AveExpr         t      P.Value  adj.P.Val
#YKR083C_5UTR_0  0.4952915  0.07625026  3.789837 0.0003836318 0.1312021
#YHR146W_5UTR_0  0.4274082  0.01341127  3.302595 0.0017095287 0.1730139
#YJL047C_5UTR_0  0.7418944 -1.02192933  3.249347 0.0019996991 0.1730139
#YGR286C_5UTR_0  0.7498502 -0.76486057  3.243675 0.0020332106 0.1730139
#YDL111C_5UTR_0 -0.5575386 -0.47685121 -3.168636 0.0025294425 0.1730139
#YLL035W_5UTR_0  0.6845004 -0.27946706  2.912272 0.0052195538 0.2577198
# not much to write home about


# how many significant without the controls:
length(which((!str_detect(rownames(topEpistasis), "unter")) & topEpistasis$"adj.P.Val" < 0.1))
# 0 at FRD10, out of 342; 5 at 20% FDR, 12 at 30%FDR

# qvalue
1- qvalue(topEpistasis$"P.Value"[!str_detect(rownames(topEpistasis), "unter")])$pi0
# pi1 = 0.3378885


#########################
# plots of individual examples
# four oligos: for each, show the replicates as points in a boxplot/line plot
# use symbols to indicate the replicate (and connect by lines)
# plot log2(RNA/DNA ratio) per replicate
# print the oligo ID below each boxplot

# for a given test, what are the oligos?
# grab them from the design
# use the recoded allele strings:
load("R_readOligos.RData")
oligoDesign <- readOligos
load("R_readOligos_UpStream.RData")
oligoDesign <- rbind(oligoDesign, readOligos)

# first sum all the counts per oligo
# above, we made meltDat:
meltDatForPlots <- gather(dat, key="sample", value="count", 5:ncol(dat)) %>% mutate(molecule = str_sub(sample, 1, 3)) %>% mutate(sample = str_sub(sample, start=5)) %>% filter(oligoBlock != "EranControl")

meltDatSummed <- meltDatForPlots %>% group_by(oligoID, sample, molecule) %>% summarize(count=sum(count)) %>% spread(key=molecule, value=count)
# from https://stackoverflow.com/questions/49404461/dividing-columns-by-particular-values-using-dplyr
meltDatSummedNorm <- meltDatSummed %>% group_by(sample) %>% mutate_at(vars(DNA, RNA), funs(./ sum(.))) %>% mutate(Expression = log2(RNA/DNA))

plotTSSEpistasis <- function(thisTest){
    #print(oligoDesign[oligoDesign$oligoBlock == thisTest,])
    theseAlleleStrings <- oligoDesign[oligoDesign$oligoBlock == thisTest,]$alleleString
    theseOligoIDs <- oligoDesign[oligoDesign$oligoBlock == thisTest,]$oligoID[order(theseAlleleStrings)] # the order turns them into 11, 12, 21, 22
    #print(theseOligoIDs)
    theseVariants <- oligoDesign[oligoDesign$oligoBlock == thisTest, ]$variantIDs[1]
    
    # now pull the data from the tables above
    datForPlot <- meltDatSummedNorm %>% filter(oligoID %in% theseOligoIDs)
    datForPlot$plotOligoName <- paste(datForPlot$oligoID, oligoDesign[datForPlot$oligoID,]$alleleString, sep="\nAllele: ")
    
    #print(data.frame(datForPlot))
    pdf(paste0(thisTest, ".pdf"), width=10, height=5, useDingbats=FALSE)
    theme_set(theme_bw())
    p <- ggplot(datForPlot, aes(x = plotOligoName, y=Expression)) +
    geom_boxplot(outlier.shape = NA) +
    #geom_jitter(position=position_jitter(width=0.2))
    geom_point(aes(color=sample)) +
    geom_line(aes(group=sample, colour=sample),alpha=0.3) +
    labs(title="", x = paste(thisTest, theseVariants, sep="\nVariants: "), y = "log2(RNA / DNA)") +
    scale_colour_viridis_d()+
    annotate("text",  x=1, y = max(datForPlot$Expression),
        label = paste("Adjusted interaction p-value: ", round(topEpistasis[thisTest,]$adj.P.Val, 2),
        "\n",
        "Interaction p-value: ", round(topEpistasis[thisTest,]$P.Value, 4),
        sep=""),
        vjust="inward", hjust="inward")

    print(p)
    dev.off()
}

plotTSSEpistasis("YKR083C_5UTR_0")

for(i in 1:5){
    plotTSSEpistasis(rownames(topEpistasis)[i])
}


# load the TSS results for comparison
load("R_toptab_200326_TSS.RData")
# now look up the results from the epistasis genes:

# DAD2
toptab[sapply(rownames(toptab), function(x){strsplit(x, "_")[[1]][1]}) == "YKR083C",]
#                        logFC     AveExpr          t      P.Value   adj.P.Val      B
#YKR083C_5UTR_0_22  0.48089883  0.194026245  5.8872357 1.817933e-06 0.0001149464  4.961418
#YKR083C_5UTR_0_12 -0.05928132 -0.075326024 -0.7185521 4.779026e-01 0.7450050717 -6.489776
#YKR083C_5UTR_0_21  0.06614310 -0.009953151  0.7127010 4.814646e-01 0.7452403290 -6.519399
# joint effect STRONGER than the addition

#RTT101
toptab[sapply(rownames(toptab), function(x){strsplit(x, "_")[[1]][1]}) == "YJL047C",]
#                       logFC    AveExpr         t      P.Value    adj.P.Val   B
#YJL047C_5UTR_0_21 -0.7129598 -0.7984431 -5.969949 1.439711e-06 9.703358e-05  5.2467238
#YJL047C_5UTR_0_12 -0.7466149 -0.8111372 -5.069409 1.849656e-05 8.019580e-04  2.8774533
#YJL047C_5UTR_0_22 -0.7164694 -0.8295717 -3.585360 1.160780e-03 2.001685e-02 -0.8927859
# either of the two variants is enough for the full effect <- nice
