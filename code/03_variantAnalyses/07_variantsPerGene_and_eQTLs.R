# compare MPRA effects with eQTLs

library(tidyverse)
library(ggrepel)
library(DescTools) # spearman CIs

# mpralm results for each variant:
load("R_designedVariantsList_200329_noControls.RData")

#################
# aggregate variant results PER GENE

#TSS
load("R_readOligos.RData")
oligoDesign <- readOligos
# UpStream:
load("R_readOligos_UpStream.RData")
oligoDesign <- rbind(oligoDesign, readOligos)

# find all the genes for which variants were assayed
assayedGenes <- unique(unlist(sapply(designedVariantsList, function(x){x$gene})))
# 2632

# for every gene, pull all its variants
# loop over the variantsList
variantsPerGene <- vector("list", length=length(assayedGenes))
names(variantsPerGene) <- assayedGenes

for(variant in names(designedVariantsList)){
    if(is.null(designedVariantsList[[variant]])){next()}
    for(i in 1:nrow(designedVariantsList[[variant]])){
        variantsPerGene[[designedVariantsList[[variant]][i, "gene"]]] <- rbind(variantsPerGene[[designedVariantsList[[variant]][i, "gene"]]], data.frame(designedVariantsList[[variant]][i, ], variant, stringsAsFactors=FALSE))
    }
}

# also pull the full-swap comparisons from the TSS:
load("mpralm_TSS_multiAllele_Tests/R_toptab_200326_TSS.RData")
toptab <- toptab[which(!str_detect(rownames(toptab), "unter")),]

TSSFullSwapsPerGene <- vector("list", length=length(assayedGenes))
names(TSSFullSwapsPerGene) <- assayedGenes
for(x in rownames(toptab)){
        thisAlleleString <- strsplit(x, "_")[[1]][length(strsplit(x, "_")[[1]])]
        if(!str_detect(thisAlleleString, "1")){
            thisGene <- strsplit(x, "_")[[1]][1]
            TSSFullSwapsPerGene[[thisGene]] <- rbind(TSSFullSwapsPerGene[[thisGene]], toptab[x, ])
        }
}

# load annotation: is the variant in a nucleosome or not?
load("R_variantsDF_181212.RData")
# this has some variants twice, if that variant was run multiple times
# collapse:
nucleosomesPerVariant <- variantsDF %>%
    select(variantID, inANucleosome) %>%
    mutate(variantID_noStrand = apply(str_split_fixed(variantsDF$variantID, "_", n=5)[,1:4], 1, function(x){paste(x, collapse="_")})) %>%
    select(variantID_noStrand, inANucleosome) %>%
    group_by(variantID_noStrand) %>%
    summarize(inANucleosome=first(inANucleosome))

inANucleosomeVector <- nucleosomesPerVariant$inANucleosome
names(inANucleosomeVector) <- nucleosomesPerVariant$variantID_noStrand

# aggregate the variants per gene in several ways
waysOfSumming <- c("biggestFC", "biggestAverageFC", "bestNomSigFC", "secondNomSigFC", "thirdNomSigFC", "bestNomSigAverageFC", "secondNomSigAverageFC", "thirdNomSigAverageFC", "bestFDRFC", "secondFDRFC", "thirdFDRFC", "bestFDRAverageFC", "secondFDRAverageFC", "thirdFDRAverageFC", "sumAll", "sumNomSig", "sumFDR", "sumUS_TSSFull", "sumUS_TSSFull_nomSig", "sumUS_TSSFull_FDR", "numberSigVarsPerGene", "numberNomSigVarsPerGene")

# we want to do the below for all variants, and for those that are/aren't in a nucleosome
# because we might imagine that genomic nucleosomes would not exist on the plasmid, and that therefore "naked" variants may correspond better with eQTL effects
# coding wise, I don't want to re-write the whole thing, and I don't want to turn it into a function either
# so note that block of filter for nucleosomes, and turn on/off accordingly

effectsPerGene <- t(sapply(assayedGenes, function(x){
    #print(x)
    theseVars <- variantsPerGene[[x]]
    # important:
    # when a variant occurs multiple times, we need to pick the best p-value for it:
    # (this should now be filtered prior to getting here)
    theseVars <- data.frame(theseVars %>% group_by(variant) %>% top_n(-1, P.Value))
    
    #print(theseVars)
    # switch on/off nucleosome filter here
    #theseVars <- theseVars[!inANucleosomeVector[theseVars$variant],]
    #print(theseVars)
    ret <- rep(NA, length(waysOfSumming))
    if(nrow(theseVars) == 0){return(ret)}
    
    averageVars <- group_by(theseVars, variant) %>% summarize(meanLogFC = mean(logFC), bestP = min(P.Value), worstP = max(P.Value), bestAdjP=min(adj.P.Val), )
    
    names(ret) <- waysOfSumming
        # biggest variant effect
        ret["biggestFC"] <- theseVars$logFC[which(abs(theseVars$logFC) == max(abs(theseVars$logFC)))]
        # biggest variant effect after averaging
        ret["biggestAverageFC"] <- averageVars$meanLogFC[which(abs(averageVars$meanLogFC) == max(abs(averageVars$meanLogFC)))]
        # most significant variant, but only if at least nominally significant
        if(length(theseVars$logFC[theseVars$P.Value == min(theseVars$P.Value) & theseVars$P.Value < 0.05]) > 0){
            ret["bestNomSigFC"] <- theseVars$logFC[theseVars$P.Value == min(theseVars$P.Value) & theseVars$P.Value < 0.05]
        }
        # extract the 2nd, 3rd most significant variant, if there is one. CAREFUL! THIS CAN BE A SECOND TEST FOR THE SAME VARIANT!
        if(length(theseVars$logFC[theseVars$P.Value < 0.05]) > 1){
            #ret["secondNomSigFC"] <- theseVars$logFC[theseVars$P.Value == min(theseVars$P.Value[theseVars$P.Value > min(theseVars$P.Value)]) & theseVars$P.Value < 0.05]
            ret["secondNomSigFC"] <- theseVars$logFC[theseVars$P.Value == sort(theseVars$P.Value[theseVars$P.Value < 0.05])[2]]
        }
        if(length(theseVars$logFC[theseVars$P.Value < 0.05]) > 2){
            ret["thirdNomSigFC"] <- theseVars$logFC[theseVars$P.Value == sort(theseVars$P.Value[theseVars$P.Value < 0.05])[3]]
        }
        # as above, but for average
        if(length(averageVars$meanLogFC[averageVars$bestP == min(averageVars$bestP) & averageVars$bestP < 0.05]) > 0){
            ret["bestNomSigAverageFC"] <- averageVars$meanLogFC[averageVars$bestP == min(averageVars$bestP) & averageVars$bestP < 0.05]
        }
        if(length(averageVars$meanLogFC[averageVars$bestP < 0.05]) > 1){
            #ret["secondNomSigAverageFC"] <- averageVars$meanLogFC[averageVars$bestP == min(averageVars$bestP[averageVars$bestP > min(averageVars$bestP)]) & averageVars$bestP < 0.05]
            ret["secondNomSigAverageFC"] <- averageVars$meanLogFC[averageVars$bestP == sort(averageVars$bestP[averageVars$bestP < 0.05])[2]]
        }
        if(length(averageVars$meanLogFC[averageVars$bestP < 0.05]) > 2){
            ret["thirdNomSigAverageFC"] <- averageVars$meanLogFC[averageVars$bestP == sort(averageVars$bestP[averageVars$bestP < 0.05])[3]]
        }
        # most significant variant, but only if 5%FDR
        if(length(theseVars$logFC[theseVars$adj.P.Val == min(theseVars$adj.P.Val) & theseVars$adj.P.Val < 0.05]) > 0){
            ret["bestFDRFC"] <- theseVars$logFC[theseVars$adj.P.Val == min(theseVars$adj.P.Val) & theseVars$adj.P.Val < 0.05]
        }
        # extract the 2nd most significant variant, if there is one. CAREFUL! THIS CAN BE A SECOND TEST FOR THE SAME VARIANT!
        if(length(theseVars$logFC[theseVars$adj.P.Val < 0.05]) > 1){
            #print(theseVars)
            # interestingly, some close p-values can result in *exactly* the same adjusted p (see YLL029W). rank by p-value but gate by adjusted!
            #ret["secondFDRFC"] <- theseVars$logFC[theseVars$P.Value == min(theseVars$P.Value[theseVars$adj.P.Val > min(theseVars$adj.P.Val)]) & theseVars$adj.P.Val < 0.05]
            #            test <- theseVars$logFC[theseVars$adj.P.Val == min(theseVars$adj.P.Val[theseVars$adj.P.Val > min(theseVars$adj.P.Val)]) & theseVars$adj.P.Val < 0.05]
            #if(length(test)>1){print(theseVars); print(test)}
            ret["secondFDRFC"] <- theseVars$logFC[theseVars$P.Value == sort(theseVars$P.Value[theseVars$adj.P.Val < 0.05])[2]]
        }
        if(length(theseVars$logFC[theseVars$adj.P.Val < 0.05]) > 2){
            ret["thirdFDRFC"] <- theseVars$logFC[theseVars$P.Value == sort(theseVars$P.Value[theseVars$adj.P.Val < 0.05])[3]]
        }
        # as above, but for average
        if(length(averageVars$meanLogFC[averageVars$bestAdjP == min(averageVars$bestAdjP) & averageVars$bestAdjP < 0.05]) > 0){
            ret["bestFDRAverageFC"] <- averageVars$meanLogFC[averageVars$bestAdjP == min(averageVars$bestAdjP) & averageVars$bestAdjP < 0.05]
        }
        if(length(averageVars$meanLogFC[averageVars$bestAdjP < 0.05]) > 1){
            #ret["secondFDRAverageFC"] <- averageVars$meanLogFC[averageVars$bestP == min(averageVars$bestP[averageVars$bestAdjP > min(averageVars$bestAdjP)]) & averageVars$bestAdjP < 0.05]
            ret["secondFDRAverageFC"] <- averageVars$meanLogFC[averageVars$bestP == sort(averageVars$bestP[averageVars$bestAdjP < 0.05])[2]]
        }
        if(length(averageVars$meanLogFC[averageVars$bestAdjP < 0.05]) > 2){
            ret["thirdFDRAverageFC"] <- averageVars$meanLogFC[averageVars$bestP == sort(averageVars$bestP[averageVars$bestAdjP < 0.05])[3]]
        }
        # sum of all variants (after averaging)
        ret["sumAll"] <- sum(averageVars$meanLogFC)
        # sum of all significant variants
        ret["sumNomSig"] <- sum(averageVars$meanLogFC[averageVars$bestP < 0.05])
        # as above but for 5%FDR
        ret["sumFDR"] <- sum(averageVars$meanLogFC[averageVars$bestAdjP < 0.05])
        #if(!is.null(TSSFullSwapsPerGene[[x]])){
        # sum of all upstream variants that are not also in TSS plus full TSS swap
        # the averageing of variant effects is always across libraries (TSS and UpStream)
        # therefore, for any variant not in TSS will have "average" == the UpStream value
        # so we can just use that
        ret["sumUS_TSSFull"] <- sum(c(theseVars$logFC[theseVars$library=="UpStream" & !theseVars$variant %in% theseVars$variant[theseVars$library=="TSS"]], TSSFullSwapsPerGene[[x]]$logFC))
        # as above, but only if nominally significant
        ret["sumUS_TSSFull_nomSig"] <- sum(c(theseVars$logFC[theseVars$library=="UpStream" & (!theseVars$variant %in% theseVars$variant[theseVars$library=="TSS"]) & theseVars$P.Value < 0.05], TSSFullSwapsPerGene[[x]]$logFC[TSSFullSwapsPerGene[[x]]$P.Value < 0.05]))
        # as above, but only if FDR 5% significant
        ret["sumUS_TSSFull_FDR"] <- sum(c(theseVars$logFC[theseVars$library=="UpStream" & (!theseVars$variant %in% theseVars$variant[theseVars$library=="TSS"]) & theseVars$adj.P.Val < 0.05], TSSFullSwapsPerGene[[x]]$logFC[TSSFullSwapsPerGene[[x]]$adj.P.Val < 0.05]))
        #}
        ret["numberSigVarsPerGene"] <- length(theseVars$logFC[theseVars$adj.P.Val < 0.05])
        ret["numberNomSigVarsPerGene"] <- length(theseVars$logFC[theseVars$P.Value < 0.05])

        # if there are no sig values, I'd rather set the result to NA. "Zero" implies a measured FC of value zero, which isn't right
        ret[ret == 0] <- NA
        ret
}))
#summary(effectsPerGene)
# 2,955 genes in this; after filter: 2632
# note that this is a larger number then the number of genes with individual variants
# ... because of genes that have a full TSS swap, but for which the single variants didn't work!
#in some cases, e.g. YEL037C, there are variants in the TSS only but the full swap failed (is not in the library)
# these are cases in which the UpStream (none in this case) plus the full swap (missing) result in NA
# there are only 24 such cases, and these all have a sum from the individual variant effects
# should be OK


#save(effectsPerGene, file="R_effectsPerGene_FDR5_200405.RData")
#save(effectsPerGene, file="R_effectsPerGene_NucleosomeFiltered_FDR5_200405.RData")
#save(effectsPerGene, file="R_effectsPerGene_NotNucleosomeFiltered_FDR5_200405.RData")


####################################
# compare to eQTL data:

# this table has both effects (ASE & local eQTLs) from the eLife paper:
load("R_allFCsWithSig_170324.RData")
allFCsWithSig <- cbind(rownames(allFCsWithSig), allFCsWithSig)
colnames(allFCsWithSig)[1] <- "geneID"
# this has 3340 genes that were present in the ASE data
# it has all eQTL info for them, but lacks eQTLs for genes that were not in the ASE data


# add data for ALL eQTLs, including from genes that were not in the ASE data
load("cisModel.effects.RData") # 2*this is "eQTL_fc" in "allFCsWithSig"
load("R_localeQTL_170222.RData")

oligosAndeQTLs <- merge(cisModel.effects * 2, effectsPerGene, by.x = "row.names", by.y = "row.names")
oligosAndeQTLs <- merge(oligosAndeQTLs, localeQTL[,c("gene", "LOD")], by.x = "Row.names", by.y="gene", all.x=TRUE)
# "LOD" is from a model with detected local eQTLs; smallest is LOD=2.5
# local_eQTL_LOD comes from the allFCsWithSig table. it is whatever the local LOD is at the gene itself, irrespective of detection
# let's use these to filter below

oligosAndeQTLs <- merge(oligosAndeQTLs, allFCsWithSig, by.x = "Row.names", by.y = "geneID", all.x=TRUE)
oligosAndeQTLs <- oligosAndeQTLs %>% select(-eQTL_fc)
colnames(oligosAndeQTLs)[1] <- "geneID"
colnames(oligosAndeQTLs)[2] <- "eQTL_fc"
dim(oligosAndeQTLs)
# 2,300 after filter

#write for supplement (without nucleosome filter only)
#write.table(oligosAndeQTLs, file="oligosAndeQTLs_200810.txt", sep="\t", quote=FALSE, row.names=FALSE)

# are there ANY combinations/batches that correlate better?
# ggplot2
#pdf("ASE_correlations_all_gg.pdf", width=11, height=6, useDingBats=FALSE)
#pdf("ASE_correlations_all_NOTnucleosomeFiltered_gg.pdf", width=11, height=6, useDingBats=FALSE)
#pdf("ASE_correlations_all_INnucleosomeFiltered_gg.pdf", width=11, height=6, useDingBats=FALSE)
par(mfrow=c(1,2))
#for (k in waysOfSumming){
for (k in c("bestFDRFC", "secondFDRFC")){
    print(k)
        for (i in c("eQTL_fc")){
            for (j in c("topOligoEffect")){try({
                # compute stats to be added to the plots:
                lodFilter <- oligosAndeQTLs[,"local_eQTL_LOD"] >= 50 & !is.na(oligosAndeQTLs[,"local_eQTL_LOD"])
                corAll <- cor.test(oligosAndeQTLs[,i], oligosAndeQTLs[,k], method="s", use="complete")
                corSig <- cor.test(oligosAndeQTLs[,i][lodFilter], oligosAndeQTLs[,k][lodFilter], method="s", use="complete")
                # were significant oligos more likely to correspond to significant local eQTLs?
                
                # do all points overlap in direction?
                datForFet <- oligosAndeQTLs[,c(i,k)]
                datForFet <- datForFet[complete.cases(datForFet),]
                fetDirection <- fisher.test(cbind(
                    c(length(which(datForFet[,1] < 0 & datForFet[,2] < 0)), length(which(datForFet[,1] < 0 & datForFet[,2] > 0))),
                    c(length(which(datForFet[,1] > 0 & datForFet[,2] < 0)), length(which(datForFet[,1] > 0 & datForFet[,2] > 0)))
                ))
                datForFeteQTL <- oligosAndeQTLs[,c(i,k)][lodFilter,]
                datForFeteQTL <- datForFeteQTL[complete.cases(datForFeteQTL),]
                fetDirectioneQTL <- fisher.test(cbind(
                    c(length(which(datForFeteQTL[,1] < 0 & datForFeteQTL[,2] < 0)), length(which(datForFeteQTL[,1] < 0 & datForFeteQTL[,2] > 0))),
                    c(length(which(datForFeteQTL[,1] > 0 & datForFeteQTL[,2] < 0)), length(which(datForFeteQTL[,1] > 0 & datForFeteQTL[,2] > 0)))
                ))

                if(corAll$p.value < 0.01){pValFormattedAll <- formatC(corAll$p.value, format = "e", digits = 1)}
                else{pValFormattedAll <- round(corAll$p.value, 3)}
                if(corSig$p.value < 0.01){pValFormattedSig <- formatC(corSig$p.value, format = "e", digits = 1)}
                else{pValFormattedSig <- round(corSig$p.value, 3)}
                if(fetDirection$p.value < 0.01){pValFormattedFet <- formatC(fetDirection$p.value, format = "e", digits = 1)}
                else{pValFormattedFet <- round(fetDirection$p.value, 3)}
                if(fetDirectioneQTL$p.value < 0.01){pValFormattedFeteQTL <- formatC(fetDirectioneQTL$p.value, format = "e", digits = 1)}
                else{pValFormattedFeteQTL <- round(fetDirectioneQTL$p.value, 3)}
                
                
                # now plot
                ggDat <- oligosAndeQTLs
                ggDat <- ggDat %>% mutate(lodFilter = local_eQTL_LOD >= 50 & !is.na(local_eQTL_LOD))
                ggDat <- ggDat %>% mutate(noeQTLButVariant = local_eQTL_LOD < 5 & !is.na(local_eQTL_LOD) & abs(ggDat[,k]) >= 0.25 )
                
                p <- ggplot(ggDat, aes_string(x=i, y=k, color="lodFilter", size="local_eQTL_LOD")) +
                geom_point() +
                theme_bw() +
                xlab("eQTL log2 fold change") + ylab("Oligo log2 fold change") + ggtitle(k) +
                scale_color_manual(values=c("#00000022", "#FF000066")) +
                geom_abline(intercept=0, slope=1, color="black", linetype="dashed") +
                geom_hline(yintercept=0, color="black", linetype="dashed") +
                geom_vline(xintercept=0, color="black", linetype="dashed") +
                #geom_text_repel(data = subset(ggDat, lodFilter), aes(label = geneName), size = 2, color="black") +
                geom_text_repel(data = subset(ggDat, noeQTLButVariant), aes(label = geneName), size = 2, color="black") +

                annotate("text", -Inf, Inf, hjust = 0, vjust = 1, label =
                    paste(
                    paste("cor all = ", round(corAll$est, 2), sep=" "),
                    paste("p all = ", pValFormattedAll, sep=" "),
                    paste("cor significant local eQTL = ", round(corSig$est, 2), sep=" "),
                    paste("p significant local eQTL = ", pValFormattedSig, sep=" "),
                    paste("number of eQTLs = ", length(which(lodFilter)), sep=" "),
                    paste("OR FET = ", round(fetDirection$est, 2), sep=" "),
                    paste("p FET = ", pValFormattedFet, sep=" "),
                    paste("OR FET eQTL = ", round(fetDirectioneQTL$est, 2), sep=" "),
                    paste("p FET eQTL = ", pValFormattedFeteQTL, sep=" "),
                    paste("directions: ", length(which(datForFeteQTL[,1] < 0 & datForFeteQTL[,2] > 0)), length(which(datForFeteQTL[,1] > 0 & datForFeteQTL[,2] > 0)), length(which(datForFeteQTL[,1] > 0 & datForFeteQTL[,2] < 0)), length(which(datForFeteQTL[,1] < 0 & datForFeteQTL[,2] < 0)), sep=" ")
                    , sep="\n")
                )

                print(p)

            })}
        }
    }
dev.off()

# get CIs for the nucleosome correlations (for plotting)

for (k in c("bestFDRFC", "sumFDR")){
    for (i in c("eQTL_fc")){
        lodFilter <- oligosAndeQTLs[,"local_eQTL_LOD"] >= 50 & !is.na(oligosAndeQTLs[,"local_eQTL_LOD"])
        #corAll <- cor.test(oligosAndeQTLs[,i], oligosAndeQTLs[,k], method="s", use="complete")
        corSig <- SpearmanRho(oligosAndeQTLs[,i][lodFilter], oligosAndeQTLs[,k][lodFilter], use="complete", conf.level=0.95)
        print(corSig)
    }
}
# NO nucleosomes:
#      rho    lwr.ci    upr.ci
#0.4879505 0.2366219 0.6780745
#      rho    lwr.ci    upr.ci
#0.5376683 0.2992501 0.7128940

# IN nucleosomes:
#       rho      lwr.ci      upr.ci
#0.31000764 -0.03145568  0.58667037
#       rho      lwr.ci      upr.ci
#0.30603514 -0.03583962  0.58378463


# plot
ggDat <- data.frame(
    name = factor(c("free, sum", "bound, sum", "free, best", "bound, best"), levels=c("free, sum", "bound, sum", "free, best", "bound, best")),
    bound = factor(c("free", "bound", "free", "bound"), levels=c("free", "bound")),
    rho = c(0.538, 0.306, 0.488, 0.31),
    rhoUpper = c(0.7128940, 0.58378463, 0.6780745, 0.58667037),
    rhoLower = c(0.2992501, -0.03583962, 0.2366219, -0.03145568)
)
pdf("nulceosomeCorsBarplots.pdf")
ggplot(ggDat) + theme_bw() +
    geom_bar( aes(x=name, y=rho, fill=bound), stat="identity", alpha=0.7) +
    geom_errorbar( aes(x=name, ymin=rhoLower, ymax=rhoUpper), width=0.4, colour="black", alpha=0.9, size=1.3)
dev.off()



############
# does number of sig variants correlate with eQTL LOD?

cor.test(oligosAndeQTLs$numberNomSigVarsPerGene, oligosAndeQTLs$LOD, method="s")
# rho=0.12, p=0.005
cor.test(oligosAndeQTLs$numberNomSigVarsPerGene, oligosAndeQTLs$local_eQTL_LOD, method="s")
# rho=0.12, p=0.003
cor.test(oligosAndeQTLs$numberNomSigVarsPerGene, abs(oligosAndeQTLs$eQTL_fc), method="s")
# rho=0.14, p=3e-5


# does the relationship persist among strong eQTLs?
cor.test(oligosAndeQTLs$numberNomSigVarsPerGene[oligosAndeQTLs$local_eQTL_LOD >= 50], abs(oligosAndeQTLs$eQTL_fc)[oligosAndeQTLs$local_eQTL_LOD >= 50], method="s")
# rho=0.20, p=0.02
cor.test(oligosAndeQTLs$numberNomSigVarsPerGene[oligosAndeQTLs$local_eQTL_LOD >= 50], oligosAndeQTLs$local_eQTL_LOD[oligosAndeQTLs$local_eQTL_LOD >= 50], method="s")
# rho=0.26, p=0.002
cor.test(oligosAndeQTLs$numberNomSigVarsPerGene[oligosAndeQTLs$local_eQTL_LOD >= 50], oligosAndeQTLs$LOD[oligosAndeQTLs$local_eQTL_LOD >= 50], method="s")
# rho=0.20, p=0.02


cor.test(oligosAndeQTLs$numberSigVarsPerGene, oligosAndeQTLs$LOD, method="s")
# rho=0.02, p=0.7
cor.test(oligosAndeQTLs$numberSigVarsPerGene, oligosAndeQTLs$local_eQTL_LOD, method="s")
# rho=0.02, p=0.8
cor.test(oligosAndeQTLs$numberSigVarsPerGene, abs(oligosAndeQTLs$eQTL_fc), method="s")
# after filter: rho=0.07, p=0.2
# not enough variation in the FDR variants it seems

# OVERALL, the filters did improve agreement with local eQTLs!





###############
# for CWP1, show the summing of effects vs the eQTL:
theseVars <- variantsPerGene[["YKL096W"]] #CWP1
theseVars <- data.frame(theseVars %>% group_by(variant) %>% top_n(-1, P.Value))
theseVars[theseVars$adj.P.Val < 0.05,][,c("logFC", "adj.P.Val", "variant")]
#       logFC   adj.P.Val           variant
#2  0.2110941 0.032794246 chrXI_260957_A_AC
#11 0.4070036 0.009716799  chrXI_260406_A_G
#12 0.3830562 0.028993177 chrXI_260785_GC_G
#13 0.7083652 0.001748200  chrXI_260911_A_G

# eQTL fc:
oligosAndeQTLs[oligosAndeQTLs$geneID == "YKL096W",c("eQTL_fc", "local_eQTL_LOD")]
#      eQTL_fc local_eQTL_LOD
#1304 2.041115       206.3517


# same for OLE1:
theseVars <- variantsPerGene[["YGL055W"]] #CWP1
theseVars <- data.frame(theseVars %>% group_by(variant) %>% top_n(-1, P.Value))
theseVars[theseVars$adj.P.Val < 0.05,][,c("logFC", "adj.P.Val", "variant")]
#       logFC  adj.P.Val           variant
#1 -0.5002607 0.01070758 chrVII_398081_A_G
oligosAndeQTLs[oligosAndeQTLs$geneID == "YGL055W",c("eQTL_fc", "local_eQTL_LOD")]
#       eQTL_fc local_eQTL_LOD
#722 -0.5133217       64.92405

# UFO1:
theseVars <- variantsPerGene[["YML088W"]] #CWP1
theseVars <- data.frame(theseVars %>% group_by(variant) %>% top_n(-1, P.Value))
theseVars[theseVars$adj.P.Val < 0.05,][,c("logFC", "adj.P.Val", "variant")]
#       logFC   adj.P.Val           variant
#5 -0.4827896 0.001822111 chrXIII_91485_G_A
#7  0.4564551 0.002559011 chrXIII_91532_T_C
oligosAndeQTLs[oligosAndeQTLs$geneID == "YML088W",c("eQTL_fc", "local_eQTL_LOD")]
# -0.03014453; LOD=1.07 <= perfect

# let's put them all in the same scale
#library(RColorBrewer)
barplot(cbind(
c(-0.5133217, 0, 0, 0), c(-0.5002607, 0, 0, 0), # OLE1
c(2.041115, 0, 0, 0), c(0.7083652, 0.4070036, 0.3830562, 0.2110941), #CWP1
c(-0.03014453, 0, 0, 0), c(-0.4827896, 0.4564551, 0, 0) #UFO1
),
col=rev(colorRampPalette(brewer.pal(4, "YlOrRd"))(4)))
