# combine the TSS and Upstream results
# do reproducibility across them
# make the volcano plot

library(tidyverse)
library(ggExtra)
library(ggrepel)

#################
# make a variant-centric list
# for each variant, catalogues results from wherever it was assayed

# load the oligo design
# TSS
load("R_readOligos.RData")
oligoDesign <- readOligos
# this for UpStream:
load("R_readOligos_UpStream.RData")
oligoDesign <- rbind(oligoDesign, readOligos)

# write out
#save(oligoDesign, file="R_oligoDesign_TSS_UpStream_190924.RData")
#write.table(oligoDesign, file="oligoDesign_TSS_UpStream_190924.txt", quote=FALSE, sep="\t", row.names=FALSE)

# first, make a list of all variants in the design

designedVariants <- unique(unlist(lapply(oligoDesign$variantIDs, function(x){strsplit(x, ";")[[1]]})))
length(designedVariants)
# 7207 for TSS & Upstream

length(which(str_detect(designedVariants, "sharon")))
# 200

# we can remove the Sharon oligos
designedVariants <- designedVariants[which(!str_detect(designedVariants, "sharon"))]
# now 9607

# for every variant, make a list results:
# whole line from toptab, as well as:
# comparison (i.e., library and oligo pair) these results come from
# strand for this comparison (can be different for bidirectional promoters)
# gene we assigned this test to

# it may be easiest to step through toptab and add results to the variant table
# load the TSS toptabs:
topTabs <- vector("list", 2)
names(topTabs) <- c("TSS", "UpStream")

load("../mpralm_TSS_multiAllele_Tests/R_toptab_200326.RData")
topTabs[["TSS"]] <- toptab
rm(toptab)
load("../mpralm_UpStream_Tests/R_toptab_200326.RData")
topTabs[["UpStream"]] <- toptab
rm(toptab)

designedVariantsList <- vector("list", length(designedVariants))
names(designedVariantsList) <- designedVariants

for (j in names(topTabs)){
    for(i in (1:nrow(topTabs[[j]]))){
        splitRow <- strsplit(rownames(topTabs[[j]])[i], "_")[[1]]
        print(splitRow)
        thisGene <- splitRow[1]
        thisAllele <- splitRow[length(splitRow)]
        #print(thisAllele)
        # if this is a block comparison, do nothing; we'll treat those separately; this is just for effects of single variants
        if((!str_detect(thisAllele, "1")) & nchar(thisAllele) > 1){next()}
        # oligoBlock is the rowname less the allele tag
        thisOligoBlock <- paste(splitRow[1:(length(splitRow)-1)], collapse="_")
        # variants in this block (split them into a vector in case there are multiple)
        theseVars <- strsplit(oligoDesign[oligoDesign$oligoBlock==thisOligoBlock, "variantIDs"][1], ";")[[1]]
        theseVarPositionsInOligo <- strsplit(oligoDesign[oligoDesign$oligoBlock==thisOligoBlock, "variantPositions"][1], ";")[[1]]
        thisNumberVariantsInOligo <- oligoDesign[oligoDesign$oligoBlock==thisOligoBlock, "numberVariants"][1]
        thisOligoStrand <- oligoDesign[oligoDesign$oligoBlock==thisOligoBlock, "strand"][1]
        #print(theseVars)
        # which variant in the block is flipped?
        #print(thisAllele)
        #print(str_locate(thisAllele, "2")[,"start"])
        whereIsTwo <- str_locate(thisAllele, "2")[,"start"]
        thisVariantInBlock <- theseVars[whereIsTwo]
        #print(thisVariantInBlock)
        thisVarPos <- theseVarPositionsInOligo[whereIsTwo]
        #print(thisVarPos)
        # now put all this into the result list
        # rownames business is to avoid weirdness during rbind (if same test re-occurs, as it can for the controls)
        resultDF <- data.frame(rownames(topTabs[[j]])[i], topTabs[[j]][i,], j, thisGene, thisOligoStrand, thisNumberVariantsInOligo, thisVarPos, whereIsTwo, stringsAsFactors=FALSE)
        names(resultDF) <- c("testName", names(resultDF[2:7]), "library", "gene", "strand", "numberVarsInOligo", "varPosInOligo", "varNumberInOligo")
        rownames(resultDF) <- c()
        #print(resultDF)
        # make sure to ADD not overwrite
        #length(designedVariantsList[[thisVariantInBlock]]) + 1
        designedVariantsList[[thisVariantInBlock]] <- rbind(designedVariantsList[[thisVariantInBlock]], resultDF)
    }
}

# how many variants were measured multiple times?
#length(designedVariantsList[names(which(unlist(sapply(designedVariantsList, nrow))>1))])
length(which(unlist(sapply(designedVariantsList, nrow))>1))
# 981
# good, that's a lot of comparisons

# how many variants were actually tested (after oligo dropout)
# 7007 entries (incl two controls)
length(which(sapply(designedVariantsList, length)>1))
# 5834
# DO NOT USE THIS NUMBER; it is inflated by the five variants tested twice on same strand same library

#save(designedVariantsList, file="R_designedVariantsList_200329.RData")
# can throw out the controls:
#designedVariantsList <- designedVariantsList[!str_detect(names(designedVariantsList), "hunter")]
#save(designedVariantsList, file="R_designedVariantsList_200329_noControls.RData")


###############################
# correlation of fold changes across replicated variants
datForEffectRepro <- c()
multiTestVariants <- designedVariantsList[names(which(unlist(sapply(designedVariantsList, nrow))>1))]

# 4/5/2020:
# update to always put + strand or TSS first, in the cases where those differ but the other not
for(i in (1:length(multiTestVariants))){
    temp <- multiTestVariants[[i]]
    apply(combn(nrow(temp), 2), 2, function(comb){
        sameStrand <- as.numeric(temp[comb[1], "strand"] == temp[comb[2], "strand"])
        sameLibrary <- as.numeric(temp[comb[1], "library"] == temp[comb[2], "library"])
        if(sameLibrary==1 & sameStrand==0){comb <- comb[c(which(temp[comb, "strand"] == "+"), which(temp[comb, "strand"] == "-"))]} # added 4/5/20
        if(sameLibrary==0 & sameStrand==1){comb <- comb[c(which(temp[comb, "library"] == "TSS"), which(temp[comb, "library"] == "UpStream"))]} # added 4/5/20
        datForEffectRepro <<- rbind(datForEffectRepro, c(temp[comb, "logFC"], temp[comb, "adj.P.Val"], sameStrand, sameLibrary, i))
    })
}
colnames(datForEffectRepro) <- c("logFC1", "logFC2", "adjP1", "adjP2", "sameStrand", "sameLibrary", "mTVIndex")
datForEffectRepro <- data.frame(datForEffectRepro)

cor.test(datForEffectRepro[,1], datForEffectRepro[,2])
#r=0.2, p=8e-12
#rho=0.11, p=0.0001

pdf("effectCorrelations_all.pdf")
        thisOne <- datForEffectRepro[,1]
        thisTwo <- datForEffectRepro[,2]
        plot(thisOne, thisTwo, main="single variant effects", xlab="logFC test 1", ylab="logFC test 2", col="#00000022")
        points(thisOne[datForEffectRepro[,3] <= 0.1], thisTwo[datForEffectRepro[,3] <= 0.1], col="#FF000044", pch=19)
        points(thisOne[datForEffectRepro[,4] <= 0.1], thisTwo[datForEffectRepro[,4] <= 0.1], col="#0000FF44", pch=19)
        abline(0, 1, col="grey", lty=2, lwd=2)
        abline(h=0, col="grey", lty=2, lwd=2)
        abline(v=0, col="grey", lty=2, lwd=2)
        legend("topleft", box.lty=0, legend=
        c(
        paste("cor = ", round(cor(thisOne, thisTwo, method="s"), 2), sep=" "),
        paste("cor significant = ", round(cor(thisOne[datForEffectRepro[,3] <= 0.1 | datForEffectRepro[,4] <= 0.1], thisTwo[datForEffectRepro[,3] <= 0.1 | datForEffectRepro[,4] <= 0.1], method="s"), 2), sep=" ")
        )
        )
dev.off()

sigLevel <- 0.05
pdf("effectCorrelations_byLibraryAndStrand.pdf")
for (i in 0:1){
    for(j in 0:1){
        subDat <- datForEffectRepro[datForEffectRepro[,"sameStrand"] == i & datForEffectRepro[,"sameLibrary"] == j,]
    thisOne <- subDat[,1]
    thisTwo <- subDat[,2]
    plot(thisOne, thisTwo, main="single variant effects", xlab="logFC test 1", ylab="logFC test 2", col="#00000022", sub=paste("same strand: ", i, "; same library: ", j, sep=""))
    if(i == 1 & j == 1){print(subDat)}
    points(thisOne[subDat[,3] <= sigLevel], thisTwo[subDat[,3] <= sigLevel], col="#FF000044", pch=19)
    points(thisOne[subDat[,4] <= sigLevel], thisTwo[subDat[,4] <= sigLevel], col="#0000FF44", pch=19)
    abline(0, 1, col="grey", lty=2, lwd=2)
    abline(h=0, col="grey", lty=2, lwd=2)
    abline(v=0, col="grey", lty=2, lwd=2)
    
    pAll <- cor.test(thisOne, thisTwo, method="s")$p.value
    if(pAll < 0.01){pValFormattedAll <- formatC(pAll, format = "e", digits = 1)}
    else{pValFormattedAll <- round(pAll, 2)}

    pSig <- cor.test(thisOne[subDat[,3] <= sigLevel | subDat[,4] <= sigLevel], thisTwo[subDat[,3] <= sigLevel | subDat[,4] <= sigLevel], method="s")$p.value
    if(pSig < 0.01){pValFormattedSig <- formatC(pSig, format = "e", digits = 1)}
    else{pValFormattedSig <- round(pSig, 2)}

    legend("topleft", box.lty=0, legend=
        c(
            paste("cor all = ", round(cor(thisOne, thisTwo, method="s"), 2), sep=" "),
            paste("p all = ", pValFormattedAll, sep=" "),
            paste("cor significant = ", round(cor(thisOne[subDat[,3] <= sigLevel | subDat[,4] <= sigLevel], thisTwo[subDat[,3] <= sigLevel | subDat[,4] <= sigLevel], method="s"), 2), sep=" "),
            paste("p significant = ", pValFormattedSig, sep=" "),        paste("n all = ", nrow(subDat)),
            paste("n sig = ", length(which(subDat[,3] <= sigLevel | subDat[,4] <= sigLevel)))
        )
    )
    }
}
dev.off()
# the five variants on the same strand same library are:
# 587: YJL149W_chrX_136809_GT_G & YJL150W_chrX_136809_GT_G => overlapping same strand ORFs, YJL150W seems not expressed
# 1284, 1323, 1335, 1337 all five are in the Upstream library for overlapping ORFs on the same strand


pdf("effectCorrelations_byLibraryOrStrand.pdf")
for (i in c("sameStrand", "sameLibrary")){
    for(j in 0:1){
        subDat <- datForEffectRepro[datForEffectRepro[,i] == j,]
        thisOne <- subDat[,1]
        thisTwo <- subDat[,2]
        plot(thisOne, thisTwo, main="single variant effects", xlab="logFC test 1", ylab="logFC test 2", col="#00000022", sub=paste(i, ": ", j, sep=""))
        points(thisOne[subDat[,3] <= 0.1], thisTwo[subDat[,3] <= 0.1], col="#FF000044", pch=19)
        points(thisOne[subDat[,4] <= 0.1], thisTwo[subDat[,4] <= 0.1], col="#0000FF44", pch=19)
        abline(0, 1, col="grey", lty=2, lwd=2)
        abline(h=0, col="grey", lty=2, lwd=2)
        abline(v=0, col="grey", lty=2, lwd=2)
        
        pAll <- cor.test(thisOne, thisTwo, method="s")$p.value
        if(pAll < 0.01){pValFormattedAll <- formatC(pAll, format = "e", digits = 1)}
        else{pValFormattedAll <- round(pAll, 2)}
        
        pSig <- cor.test(thisOne[subDat[,3] <= 0.1 | subDat[,4] <= 0.1], thisTwo[subDat[,3] <= 0.1 | subDat[,4] <= 0.1], method="s")$p.value
        print(pSig)
        if(pSig < 0.01){pValFormattedSig <- formatC(pSig, format = "e", digits = 1)}
        else{pValFormattedSig <- round(pSig, 2)}
        
        legend("topleft", box.lty=0, legend=
        c(
        paste("cor all = ", round(cor(thisOne, thisTwo, method="s"), 2), sep=" "),
        paste("p all = ", pValFormattedAll, sep=" "),
        paste("cor significant = ", round(cor(thisOne[subDat[,3] <= 0.1 | subDat[,4] <= 0.1], thisTwo[subDat[,3] <= 0.1 | subDat[,4] <= 0.1], method="s"), 2), sep=" "),
        paste("p significant = ", pValFormattedSig, sep=" "),
        paste("n all = ", nrow(subDat)),
        paste("n sig = ", length(which(subDat[,3] <= 0.1 | subDat[,4] <= 0.1)))
        )
        )
    }
}
dev.off()
# correlation for sig is rho=0.33, p=5e-7

# do FETs on same direction
for (i in 0:1){
    for(j in 0:1){
        print(paste0("samestrand: ", i, "sameLibrary: ", j))
        subDat <- datForEffectRepro[datForEffectRepro[,"sameStrand"] == i & datForEffectRepro[,"sameLibrary"] == j,]
        #subDat <- datForEffectRepro[datForEffectRepro[,"sameLibrary"] == j,]
        FETDirection <- fisher.test(
            cbind(
                c(length(which(subDat$logFC1 > 0 & subDat$logFC2 > 0)), length(which(subDat$logFC1 > 0 & subDat$logFC2 < 0))),
                c(length(which(subDat$logFC1 < 0 & subDat$logFC2 > 0)), length(which(subDat$logFC1 < 0 & subDat$logFC2 < 0)))
            )
        )
        print(FETDirection)

        FETMatForSig <- cbind(
            c(length(which(subDat$adjP1 > 0.05 & subDat$adjP2 > 0.05)), length(which(subDat$adjP1 > 0.05 & subDat$adjP2 < 0.05))),
            c(length(which(subDat$adjP1 < 0.05 & subDat$adjP2 > 0.05)), length(which(subDat$adjP1 < 0.05 & subDat$adjP2 < 0.05)))
        )

        FETSignificance <- fisher.test(FETMatForSig)
        print(FETMatForSig)
        print(FETSignificance)
        }
}


# do FETs on same direction AMONG SIGNIFICANT VARS
for (i in 0:1){
    for(j in 0:1){
        print(paste0("samestrand: ", i, "sameLibrary: ", j))
        subDat <- datForEffectRepro[datForEffectRepro[,"sameStrand"] == i & datForEffectRepro[,"sameLibrary"] == j,]
        subDat <- subDat[subDat$adjP1 < 0.05 | subDat$adjP2 < 0.05,]
        #subDat <- datForEffectRepro[datForEffectRepro[,"sameLibrary"] == j,]
        FETDirection <- fisher.test(
            cbind(
                c(length(which(subDat$logFC1 > 0 & subDat$logFC2 > 0)), length(which(subDat$logFC1 > 0 & subDat$logFC2 < 0))),
                c(length(which(subDat$logFC1 < 0 & subDat$logFC2 > 0)), length(which(subDat$logFC1 < 0 & subDat$logFC2 < 0)))
            )
        )
        print(FETDirection)
        }
}



#####################
# how many sig variants, using the bestP if run several times?

bestPs <- t(sapply(designedVariantsList, function(x){
    if(is.null(x)){return(c(NA, NA, NA, NA, NA, NA))}
    minPValRow <- x[which(x$P.Value == min(x$P.Value)),][1,]
    minAdjPRow <- x[which(x$adj.P.Val == min(x$adj.P.Val)),][1,]
    TSSOrUSminP <- 0
    TSSOrUSminAdjP <- 0
    if(minPValRow$library == "UpStream"){TSSOrUSminP <- 1}
    if(minAdjPRow$library == "UpStream"){TSSOrUSminAdjP <- 1}
    c(minPValRow$P.Value, minAdjPRow$adj.P.Val, minPValRow$logFC, minAdjPRow$logFC, TSSOrUSminP, TSSOrUSminAdjP)
}))
colnames(bestPs) <- c("minPVal", "minAdjPVal", "FCMinPVal", "FCMinAdjPVal", "TSSOrUSminP", "TSSOrUSminAdjP")
bestPs <- data.frame(bestPs)
bestPs <- bestPs[complete.cases(bestPs),]
bestPs$variant <- rownames(bestPs)

# 5832 were run
length(which(!is.na(bestPs[,"minPVal"])))
# 5832 with data
length(which(bestPs[,"minPVal"] < 0.05))
# 1289
length(which(bestPs[,"minAdjPVal"] < 0.05))
# 451 at 5%FDR


# their magnitude of effects:
# I'm 90% sure the "log" in MPRA package is log2
summary(2^(abs(bestPs[,"FCMinAdjPVal"][bestPs[,"minAdjPVal"] < 0.05])))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#1.090   1.241   1.347   1.464   1.516   8.862
length(which(abs(bestPs[,"FCMinAdjPVal"][bestPs[,"minAdjPVal"] < 0.05]) >= 1))
# 29, out of:
length(abs(bestPs[,"FCMinAdjPVal"][bestPs[,"minAdjPVal"] < 0.05]))
# 451=> 6%

# how many come from the TSS vs UpStream library?
summary(bestPs[,"TSSOrUSminAdjP"] == 1)
#   Mode   FALSE    TRUE
#logical    1998    3834
# more from Upstream, but not exclusively so

# plot these fold changes
pdf("foldChanges.pdf", width=10, height=5)
par(mfrow=c(1,2))
hist(abs(bestPs[,"FCMinAdjPVal"][bestPs[,"minAdjPVal"] < 0.05 & !is.na(bestPs[,2])]), breaks=100, xlab="Absolute log2 fold change", main="Effect sizes")
abline(v = median(abs(bestPs[,"FCMinAdjPVal"][bestPs[,"minAdjPVal"] < 0.05 & !is.na(bestPs[,2])])), col="grey", lwd=2, lty=2)

hist(bestPs[,"FCMinAdjPVal"][bestPs[,"minAdjPVal"] < 0.05 & !is.na(bestPs[,2])], breaks=100, xlab="Log2 fold change", main="Effect sizes")

dev.off()


######################
# volcano plot

plotCol <- rep("#00000022", nrow(bestPs))
plotCol[bestPs[,"TSSOrUSminAdjP"] == 1] <- "#0575BC22"
plotCol[bestPs[,"TSSOrUSminAdjP"] == 0 & bestPs[,"minAdjPVal"] < 0.05] <- "black"
plotCol[bestPs[,"TSSOrUSminAdjP"] == 1 & bestPs[,"minAdjPVal"] < 0.05] <- "blue"
plotCex = rep(0.5, nrow(bestPs))
plotCex[bestPs[,"minAdjPVal"] < 0.05] <- 1
plotPch = rep(19, nrow(bestPs))
plotPch[bestPs[,"minAdjPVal"] < 0.05] <- 1

plot(bestPs[,"FCMinAdjPVal"], -log10(bestPs[,"minPVal"]), col=plotCol, pch=plotPch, cex=plotCex, xlab="log2 fold change", ylab="-log10(p-value)")

# examples:
# OLE1:
designedVariantsList[["chrVII_398081_A_G"]]
#                        testName      logFC  AveExpr         t      P.Value  adj.P.Val         B  library    gene strand numberVarsInOligo
#1 YGL055W_chrVII_398081_A_G_72_2 -0.5002607 1.790179 -4.387375 0.0004216676 0.01070758 0.1378976 UpStream YGL055W      +                 1
# expected direction; significant


# SFA1 Amanda https://pubs.acs.org/doi/full/10.1021/acssynbio.6b00264 and https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-13-479
# MSN2/4 motif with the C, not with the T; C should be higher
designedVariantsList[["chrIV_159428_T_C"]]
#                      testName     logFC   AveExpr        t      P.Value    adj.P.Val         B  library    gene strand numberVarsInOligo
#1             YDL168W_5UTR_0_21 0.6198072 0.3725458 8.585917 1.256107e-09 1.588452e-07 11.997778      TSS YDL168W      +                 2
#2 YDL168W_chrIV_159428_T_C_72_2 0.5706761 1.0974803 6.223105 1.023640e-05 8.109029e-04  3.540509 UpStream YDL168W      +                 1
# nice

exampleVars <- c("chrVII_398081_A_G", "chrIV_159428_T_C")
names(exampleVars) <- c("YGL055W", "YDL168W")

# prettier with ggplot:
# using, in part, http://www.sthda.com/english/wiki/ggplot2-texts-add-text-annotations-to-a-graph-in-r-software

geneAnnotation = read.table("ensemblGenes_ensembl83_160307_MOD.txt", sep="\t", stringsAsFactors=FALSE, head=TRUE)
rownames(geneAnnotation) = geneAnnotation[,1]
allNames <- geneAnnotation[, "geneName"]
names(allNames) <- geneAnnotation[,1]
allNames[which(allNames == "")] <- names(allNames)[which(allNames == "")]
allNamesInv <- names(allNames)
names(allNamesInv) <- allNames

bestPs <- bestPs %>% mutate(Library = ifelse(TSSOrUSminP == 1, "Upstream", "TSS"))
#bestPs$plotSize <- 0.5
#bestPs$plotSize[bestPs$minAdjPVal < 0.05] <- -log10(bestPs$minAdjPVal[bestPs$minAdjPVal < 0.05])
bestPs$MinusLog10P <- -log10(bestPs$minPVal)
bestPs$Significant <- FALSE
bestPs$Significant[bestPs$minAdjPVal < 0.05] <- TRUE
bestPs$Gene <- sapply(bestPs$variant, function(x){designedVariantsList[[x]]$gene[1]})
bestPs$Variant <- paste(allNames[bestPs$Gene], bestPs$variant, sep=":")
# smallest P that is significant adjusted P
minSigP <- bestPs$minPVal[which(bestPs$minAdjPVal == max(bestPs$minAdjPVal[bestPs$minAdjPVal < 0.05]))]

# write for suppl table
write.table(bestPs, file="bestPs_200810.txt", sep="\t", quote=FALSE, row.names=FALSE)


# how many genes?
length(unique(bestPs$Gene))
# 2503

p <- ggplot(bestPs, aes(x=FCMinAdjPVal, y=-log10(minPVal), color=Significant, size=MinusLog10P, shape=Library)) +
    geom_point() +
    scale_color_manual(values=c("#00000044", "#0575BC44")) +
    scale_shape_manual(values=c(19, 15))+
    theme_bw() +
    xlab("log2 fold change") + ylab("-log10(p value)") +
    geom_hline(yintercept=-log10(minSigP), colour="#990000", linetype="dashed") +
    geom_text_repel(data = subset(bestPs, -log10(minPVal) > 8 | abs(FCMinAdjPVal) >= 1.5), aes(label = Variant), size = 3, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"), color="black") +
    geom_text_repel(data = subset(bestPs, Gene == "YGL055W" & Significant), aes(label = Variant), size = 5, color="red", box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines")) +
    geom_text_repel(data = subset(bestPs, Gene == "YDL168W" & Significant), aes(label = Variant), size = 5, color="orange", box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines")) +
    geom_text_repel(data = subset(bestPs, Gene == "YKL096W" & Significant), aes(label = Variant), size = 5, color="#808000", box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines")) +
    geom_text_repel(data = subset(bestPs, Gene == "YML088W" & Significant), aes(label = Variant), size = 5, color="#990000", box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines")) +
    geom_point(data=filter(bestPs, Gene == "YGL055W"), color="red") +
    geom_point(data=filter(bestPs, Gene == "YDL168W"), color="orange") +
    geom_point(data=filter(bestPs, Gene == "YKL096W"), color="#808000") +
    geom_point(data=filter(bestPs, Gene == "YML088W"), color="#990000")
ggMarginal(p, margins='x', groupFill=TRUE, groupColour=TRUE)
