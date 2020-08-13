# statistical tests on the Upstream library
# AND correlations among samples, with Sharon2012, and the genome

library(mpra)
library(tidyverse)
library(qvalue)
library(gplots) # heatmap
library(GGally) # ggpairs
library(RColorBrewer)
library(parallel)
library(corrplot)
#library(MASS)
library(viridis)

load("R_countsMappedToDesignedOligos_180616.RData")

############################
# correlations among samples

datForSampleChecks <- countsMappedToDesignedOligos[["UpStream:November2016"]]
# group by barcodes; who cares about their individual counts
datForSampleChecks <- datForSampleChecks[,c("oligoID", "count", colnames(datForSampleChecks)[20:31])]
colnames(datForSampleChecks)[2] <- "annotationCount"
# replace the two failed DNAs
datForSampleChecks[,"DNA_UpStream_November2016_upstream_2018_2"] <- datForSampleChecks[,"DNA_UpStream_November2016_upstream_2018_1"]
datForSampleChecks[,"DNA_UpStream_November2016_upstream_2018_4"] <- datForSampleChecks[,"DNA_UpStream_November2016_upstream_2018_3"]

datByOligo <- datForSampleChecks %>% group_by(oligoID) %>% summarize_all(sum)

# what about the RNA/DNA ratios?
# let's see if we can do this using tidy:
datByOligoGathered <- datByOligo %>% dplyr::select(-c("annotationCount"))
datByOligoGathered <- datByOligoGathered %>% gather(key="sample", value="readCount", 2:ncol(datByOligoGathered))
datByOligoGathered <- mutate(datByOligoGathered, molecule = str_sub(sample, 1, 3), sample=str_sub(sample, 5))

datByOligoSpread <- datByOligoGathered %>% spread(key=molecule, value=readCount) %>% mutate(logRatio = log2(RNA / DNA)) %>% dplyr::select(-one_of(c("DNA", "RNA"))) %>% spread(key=sample, value=logRatio)
colnames(datByOligoSpread) <- gsub("_all", "", colnames(datByOligoSpread))
colnames(datByOligoSpread) <- sapply(colnames(datByOligoSpread), function(x){paste0(strsplit(x, "_")[[1]][(length(strsplit(x, "_")[[1]]) - 1):length(strsplit(x, "_")[[1]])], collapse="_")})
datByOligoSpread <- datByOligoSpread %>% mutate_all(funs(replace(., is.infinite(.), NA)))
datByOligoSpread <- datByOligoSpread %>% mutate_all(funs(replace(., is.nan(.), NA)))
# save for across library comparisons
#save(datByOligoSpread, file="R_datByOligoSpread_Upstream_200718.RData")


# prettier with corrplot:
pdf("ratioPairwise_corrplot.pdf", width=10, height=10)
corrplot.mixed(cor(datByOligoSpread[,2:ncol(datByOligoSpread)], use="complete", method="s"),
    upper="number",
    lower="ellipse"
)
dev.off()

# get the pairwise correlations for analyses:
pairwiseCors <- cor(datByOligoSpread[,2:ncol(datByOligoSpread)], method="s", use="pairwise.complete.obs")
summary(pairwiseCors[lower.tri(pairwiseCors)])
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.7927  0.8166  0.8292  0.8462  0.8808  0.9032

# compare them to Eran 2012 as well:
# need the design:
oligoDesign <- read.table("jointDesignFinal27k_150612_ReverseComplementWhenALargerT.txt", sep="\t", head=TRUE, stringsAsFactors=FALSE)
rownames(oligoDesign) <- oligoDesign$oligoID
# for Upstream, the Eran controls are:
EranControls <- oligoDesign[str_detect(oligoDesign$oligoID, "singleVariantTiles") & oligoDesign$oligoBlock == "EranControl",]

# pull the Sharon 2012 data
eranData = read.table("Sharon_2012_nbt_S2_mod.txt", sep="\t", stringsAsFactors=FALSE, head=TRUE)
eranData = eranData[complete.cases(eranData),]
eranData[,1] = paste("ID_sharon2012", eranData[,1], sep="_")
# now the IDs in this table should match the entries in the "variantIDs" column of our design

# attach whether the oligo contains a GAL4 binding site
hasGal4 = str_detect(eranData$Description, "Name=GAL4")
hasStrongGal4 = str_detect(eranData$Description, "Name=GAL4_S_p")
eranData <- cbind(eranData, hasGal4, hasStrongGal4)

# join Eran and our data
eranToOurDesign <- merge(eranData, EranControls, by.x="LibraryID", by.y="variantIDs")
eranAndOurs <- merge(eranToOurDesign, datByOligoSpread, by.x="oligoID", by.y="oligoID")

# now step through each sample and plot:
pdf("eranComparison.pdf")
for(i in 23:ncol(eranAndOurs)){
    plot(
        (log2(eranAndOurs[,"ExpressionReplicate1"]) + log2(eranAndOurs[,"ExpressionReplicate2"]))/2,
        eranAndOurs[,i],# xlim=c(-2,2), ylim=c(-2,2),
        xlab="Sharon 2012 expression", ylab="This experiment", main=colnames(eranAndOurs)[i], col="lightgrey"
    )
    points(
        ((log2(eranAndOurs[,"ExpressionReplicate1"]) + log2(eranAndOurs[,"ExpressionReplicate2"]))/2)[!eranAndOurs$hasStrongGal4 & !eranAndOurs$hasGal4],
        eranAndOurs[,i][!eranAndOurs$hasStrongGal4 & !eranAndOurs$hasGal4],
        pch=19, col="blue"
    )
    thisCor <- cor.test(
        ((log2(eranAndOurs[,"ExpressionReplicate1"]) + log2(eranAndOurs[,"ExpressionReplicate2"]))/2)[!eranAndOurs$hasStrongGal4 & !eranAndOurs$hasGal4],
        eranAndOurs[,i][!eranAndOurs$hasStrongGal4 & !eranAndOurs$hasGal4], method="s")
        legend("topleft", legend=c(paste("rho = ", round(thisCor$est, 2), sep=""), paste("p = ", round(thisCor$p.value, 6), sep="")))
        print(thisCor)
    abline(lm(eranAndOurs[,i][!eranAndOurs$hasStrongGal4 & !eranAndOurs$hasGal4] ~ ((log2(eranAndOurs[,"ExpressionReplicate1"]) + log2(eranAndOurs[,"ExpressionReplicate2"]))/2)[!eranAndOurs$hasStrongGal4 & !eranAndOurs$hasGal4]), col="grey", lty=2, lwd=2)
}
dev.off()
# they are all significantly and positively correlated!
# all p < 2.2e-16
# this library does have the HIS3 promoter, and it replicates better than the TSS one
# very nice

# average of all samples against Eran:
pdf("eranComparison_average.pdf")
    plot(
        (log2(eranAndOurs[,"ExpressionReplicate1"]) + log2(eranAndOurs[,"ExpressionReplicate2"]))/2,
    rowMeans(eranAndOurs[,23:ncol(eranAndOurs)]),
        xlab="Sharon 2012 expression", ylab="This experiment", main="average of our samples", col="lightgrey"
    )
    points(
        ((log2(eranAndOurs[,"ExpressionReplicate1"]) + log2(eranAndOurs[,"ExpressionReplicate2"]))/2)[!eranAndOurs$hasStrongGal4 & !eranAndOurs$hasGal4],
        rowMeans(eranAndOurs[,23:ncol(eranAndOurs)])[!eranAndOurs$hasStrongGal4 & !eranAndOurs$hasGal4],
        pch=19, col="blue"
    )
    thisCor <- cor.test(
        ((log2(eranAndOurs[,"ExpressionReplicate1"]) + log2(eranAndOurs[,"ExpressionReplicate2"]))/2)[!eranAndOurs$hasStrongGal4 & !eranAndOurs$hasGal4],
        rowMeans(eranAndOurs[,23:ncol(eranAndOurs)])[!eranAndOurs$hasStrongGal4 & !eranAndOurs$hasGal4], method="s")
        legend("topleft", legend=c(paste("rho = ", round(thisCor$est, 2), sep=""), paste("p = ", round(thisCor$p.value, 6), sep="")))
#        print(thisCor)
    abline(lm(rowMeans(eranAndOurs[,23:ncol(eranAndOurs)])[!eranAndOurs$hasStrongGal4 & !eranAndOurs$hasGal4] ~ ((log2(eranAndOurs[,"ExpressionReplicate1"]) + log2(eranAndOurs[,"ExpressionReplicate2"]))/2)[!eranAndOurs$hasStrongGal4 & !eranAndOurs$hasGal4]), col="grey", lty=2, lwd=2)
dev.off()
# rho=0.54, p<2.2e-16

# compare against gene expression of the native gene
# make average of all oligos for the given gene (assuming variant effects are minor compared to between-gene variation)
# attach the gene to each oligo, take the mean
oligoExpressionByGene <- datByOligoSpread %>%
    add_column(gene= str_split_fixed(oligoDesign[datByOligoSpread$oligoID,]$gene, "_", n=2)[,1]) %>%
    select(-oligoID) %>%
    filter(gene != "control") %>%
    group_by(gene) %>%
    summarise_all(funs(mean(., na.rm = TRUE)))

# load the eLife2018 data
eLifeExpression <- read.table("../../annotations/AlbertBloom_2018/SI_Data_01_expressionValues.txt", head=TRUE, sep="\t")
geneExpression <- colMeans(eLifeExpression)

expressionComparison <- add_column(oligoExpressionByGene, oligoExpression=rowMeans(oligoExpressionByGene[,2:ncol(oligoExpressionByGene)]) , geneExpression = geneExpression[oligoExpressionByGene$gene])

cor.test(expressionComparison$geneExpression, expressionComparison$oligoExpression, method="s")
#rho=0.234, p<2.2e-16
# huh, even UpStream has signal about relative gene expression. very nice.

pdf("oligosVsGeneExpression.pdf", width=5, height=5)
plot(expressionComparison$geneExpression, expressionComparison$oligoExpression, col="#00000044", xlab="mRNA expression (log2(TPM))", ylab="oligo expression (log2(RNA/DNA))")
#abline(lm(expressionComparison$oligoExpression ~ expressionComparison$geneExpression), col="grey", lwd=2, lty=2)
dev.off()

# more pretty
# from https://slowkow.com/notes/ggplot2-color-by-density/
# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
expressionComparisonComplete <- expressionComparison[complete.cases(expressionComparison),]
expressionComparisonComplete$plottingDensity <- get_density(expressionComparisonComplete$geneExpression, expressionComparisonComplete$oligoExpression, n = 100)

pdf("oligosVsGeneExpression_densities.pdf", width=7, height=5, useDingbats=FALSE)
theme_set(theme_bw(base_size = 16))
ggplot(expressionComparisonComplete, aes(x=geneExpression, y=oligoExpression, color = plottingDensity)) + geom_point() + scale_color_viridis() +
    geom_smooth(method=lm) +
    xlab("mRNA expression (log2(TPM))") +
    ylab("Oligo expression (log2(RNA/DNA))")
dev.off()




###########################
# statistical tests
load("R_countsMappedToDesignedOligos_180616.RData")

dat <- countsMappedToDesignedOligos[["UpStream:November2016"]]
dat <- dat[,c(1, 7, 8, 16, 20:31)] # '3' is count; might need to remove; is remove start with '5' below
dat[,"DNA_UpStream_November2016_upstream_2018_2"] <- dat[,"DNA_UpStream_November2016_upstream_2018_1"]
dat[,"DNA_UpStream_November2016_upstream_2018_4"] <- dat[,"DNA_UpStream_November2016_upstream_2018_3"]

dat$alleleString[!str_detect(dat$alleleString, "2")] <- "1"

# throw out oligos with zero counts
# first sum the barcodes
datSummed <- dat %>% group_by(oligoID) %>% summarize_at(5:ncol(dat), sum)
# which oligos have no zero counts ever, anywhere?
keepThese <- datSummed$oligoID[which(!apply(datSummed[,2:ncol(datSummed)], 1, function(x){0 %in% x}))]
length(keepThese)
# keep 9253 out of 9592; very little loss so far; this cost 1/4 of data for epistasis (which has 12 samples, not 6)
dat <- dat[dat$oligoID %in% keepThese,]
# down to 5,642,447 from 5,644,368

# write out datSummed for supplement
# in resulting file, delete DNA counts for 2 & 4 to indictae that these were replaced
write.table(datSummed, file="datSummed_Upstream_200810.txt", sep="\t", quote=FALSE, row.names=FALSE)


# because of the mpralm fluke in which it apparently works with oligos that HAVE NO DATA, we need to keep only data in which all four oligos are actually present!
presentOligos <- dat %>% group_by(oligoBlock) %>% summarize(nOligo = length(unique(oligoID)))
completeOligoBlocks <- presentOligos %>% filter(nOligo == 2)
# 4469 such blocks, out of 4582
dat <- dat[dat$oligoBlock %in% completeOligoBlocks$oligoBlock,]
# 4,208,841 rows


meltDat <- gather(dat, key="sample", value="count", 5:ncol(dat)) %>% mutate(molecule = str_sub(sample, 1, 3)) %>% mutate(sample = str_sub(sample, start=5)) %>% mutate(sample = paste(sample, alleleString, sep = ";allele_")) %>% filter(oligoBlock != "EranControl")

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
        theseallBYCounts <- allBYCounts[[x]] %>% filter(oligoBlock %in% theseAlleleCounts$oligoBlock)
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
#save(pairwiseCounts, file="R_pairwiseCounts_200326.RData")


# glue the many matrixes together into the 'dna' and 'rna' objects
DNACounts <- bind_rows(pairwiseCounts[["DNA"]])
RNACounts <- bind_rows(pairwiseCounts[["RNA"]])
dim(DNACounts)
# 4,208,841, used to be: 4,214,870

dna <- DNACounts[,c(6:ncol(DNACounts))]
rna <- RNACounts[,c(6:ncol(RNACounts))]
eid <- DNACounts$oligoBlock
design <- data.frame(intercept=1, allele=str_sub(str_extract(colnames(dna), ";(.)*"), 2))
design[,2] <- as.numeric(design[,2])
sampleBlocks <- str_sub(str_extract(colnames(dna), "(.)*;"), end = -2)

mpraset <- MPRASet(DNA = dna, RNA = rna, eid = eid, eseq = NULL, barcode = NULL)

pdf("mpralm_200326.pdf")
mpralm_fit <- mpralm(object = mpraset, design = design, aggregate = "sum", normalize = TRUE, model_type = "corr_groups", block=sampleBlocks, plot = TRUE)
dev.off()

toptab <- topTable(mpralm_fit, coef = 2, number = Inf)
toptab <- toptab[order(toptab$P.Value),]
#save(toptab, file="R_toptab_200326.RData")



# how many significant without the controls:
length(which((!str_detect(rownames(toptab), "unter")) & toptab$"adj.P.Val" < 0.1))
# 453 FDR10, 300 FDR5

# qvalue
1- qvalue(toptab$"P.Value"[!str_detect(rownames(toptab), "unter")])$pi0
# 0.3102795 same!

# how many unique variants is this?

length(unique(apply(str_split_fixed(rownames(toptab)[which((!str_detect(rownames(toptab), "unter")) & toptab$"adj.P.Val" < 0.05)], "_", n=6)[,2:5], 1, function(x){str_c(x, collapse="_")})))
# 420 unique at FDR10, 300 at FDR5
# after filter:
# 442 at FDR10, 293 at FDR5

# how many tests were actually performed (after filter)?
nrow(toptab[!str_detect(rownames(toptab), "unter"),])
#4467
# for how many genes?
length(unique(sapply(rownames(toptab[!str_detect(rownames(toptab), "unter"),]), function(x){str_split(x, "_")[[1]][1]})))
# 1824


# pi1?
1-qvalue(toptab$P.Value)$pi0
#0.3103615


###############################
# plots for single variants:

# TSS
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


# OLE1:
plotSingleVariantTSS("YGL055W_chrVII_398081_A_G_72_2")
