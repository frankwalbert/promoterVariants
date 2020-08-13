# correlations among TSS samples, as well as to Sharon 2012 and to genomic gene expression

library(mpra)
library(tidyverse)
library(qvalue)
library(gplots) # heatmap
library(GGally) # ggpairs
library(RColorBrewer)
library(corrplot)
#library(MASS) will cause "select" to collide
library(viridis)

load("R_countsMappedToDesignedOligos_180616.RData")


# for each sample:
# correlation plot DNA with annotation
# barcode and oligo level
# correlations DNA/RNA, barcode and annotation

datForSampleChecks <- countsMappedToDesignedOligos[["TSS:scaleUpJune2016"]]
# group by barcodes; who cares about their individual counts
datForSampleChecks <- datForSampleChecks[,c("oligoID", "count", colnames(datForSampleChecks)[20:47])]
colnames(datForSampleChecks)[2] <- "annotationCount"
datForSampleChecks[,"DNA_TSS_scaleUpJune2016_replicates_2018_A"] <- datForSampleChecks[,"DNA_TSS_scaleUpJune2016_replicates_2017_A"]
datForSampleChecks[,"DNA_TSS_scaleUpJune2016_replicates_2018_B"] <- datForSampleChecks[,"DNA_TSS_scaleUpJune2016_replicates_2017_B"]

datByOligo <- datForSampleChecks %>% group_by(oligoID) %>% summarize_all(sum)
datByOligoMerged <- datByOligo # historic reasons for this
datByOligo <- datByOligoMerged %>% select(-matches("RNA"))


# but what about the RNA/DNA ratios?
# are they really as poorly correlated as indicated by the fold changes (which are for alleles, but let's start with oligos for now)?
# let's see if we can do this using tidy:
datByOligoGathered <- datByOligoMerged %>% select(-one_of(c("annotationCount")))
datByOligoGathered <- datByOligoGathered %>% gather(key="sample", value="readCount", 2:ncol(datByOligoGathered))
datByOligoGathered <- mutate(datByOligoGathered, molecule = str_sub(sample, 1, 3), sample=str_sub(sample, 5))

datByOligoSpread <- datByOligoGathered %>% spread(key=molecule, value=readCount) %>% mutate(logRatio = log2(RNA / DNA)) %>% select(-one_of(c("DNA", "RNA"))) %>% spread(key=sample, value=logRatio)
colnames(datByOligoSpread) <- gsub("_all", "", colnames(datByOligoSpread))
colnames(datByOligoSpread) <- sapply(colnames(datByOligoSpread), function(x){paste0(strsplit(x, "_")[[1]][(length(strsplit(x, "_")[[1]]) - 1):length(strsplit(x, "_")[[1]])], collapse="_")})
datByOligoSpread <- datByOligoSpread %>% mutate_all(funs(replace(., is.infinite(.), NA)))
datByOligoSpread <- datByOligoSpread %>% mutate_all(funs(replace(., is.nan(.), NA)))
# save for across library comparisons
#save(datByOligoSpread, file="R_datByOligoSpread_TSS_200718.RData")


# NOTE: to get the cor to plot, there cannot be Inf in the table
pdf("ratioPairwise_corrplot.pdf", width=10, height=10)
corrplot.mixed(cor(datByOligoSpread[,10:ncol(datByOligoSpread)], use="complete", method="s"),
    upper="number",
    lower="ellipse"
)
dev.off()


# get the pairwise correlations for analyses:
pairwiseCors <- cor(datByOligoSpread[,2:ncol(datByOligoSpread)], method="s", use="pairwise.complete.obs")
summary(pairwiseCors[lower.tri(pairwiseCors)])
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.5508  0.6721  0.7357  0.7284  0.7766  0.9273

# compare them to Eran 2012 as well:
# need the design:
oligoDesign <- read.table("jointDesignFinal27k_150612_ReverseComplementWhenALargerT.txt", sep="\t", head=TRUE, stringsAsFactors=FALSE)
rownames(oligoDesign) <- oligoDesign$oligoID
# for TSS, the Eran controls are:
EranControls <- oligoDesign[str_detect(oligoDesign$oligoID, "firstTileUTR") & oligoDesign$oligoBlock == "EranControl",]

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
    abline(lm(eranAndOurs[,i][!eranAndOurs$hasStrongGal4 & !eranAndOurs$hasGal4] ~ ((log2(eranAndOurs[,"ExpressionReplicate1"]) + log2(eranAndOurs[,"ExpressionReplicate2"]))/2)[!eranAndOurs$hasStrongGal4 & !eranAndOurs$hasGal4]), col="grey", lty=2, lwd=2)
}
dev.off()
# they are all significantly and positively correlated!
# and this is although TSS did not have the HIS3 promoter

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
# rho=0.403, p=8e-7

# compare against gene expression of the native gene
# make average of all oligos for the given gene (assuming variant effects are minor compared to between-gene variation)
# attach the gene to each oligo, take the mean
oligoExpressionByGene <- datByOligoSpread %>%
    #select(-c(colnames(datByOligoSpread)[2:9])) %>%
    add_column(gene= str_split_fixed(oligoDesign[datByOligoSpread$oligoID,]$gene, "_", n=2)[,1]) %>%
    select(-oligoID) %>%
    filter(gene != "control") %>%
    group_by(gene) %>%
    summarise_all(funs(mean(., na.rm = TRUE)))

# load the eLife2018 data
eLifeExpression <- read.table("AB18_SI_Data_01_expressionValues.txt", head=TRUE, sep="\t")
geneExpression <- colMeans(eLifeExpression)

expressionComparison <- add_column(oligoExpressionByGene, oligoExpression=rowMeans(oligoExpressionByGene[,2:ncol(oligoExpressionByGene)]) , geneExpression = geneExpression[oligoExpressionByGene$gene])

cor.test(expressionComparison$geneExpression, expressionComparison$oligoExpression, method="s")
#rho=0.409, p<2.2e-16


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

