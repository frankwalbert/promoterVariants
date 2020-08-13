# correlate the Sharon/Segal control oligos within AND between libraries

library(tidyverse)
library(corrplot)


# load the two datasets (need to make these from their respective scripts)
load("R_datByOligoSpread_Upstream_200718.RData")
UpstreamDat <- datByOligoSpread
rm(datByOligoSpread)

load("R_datByOligoSpread_TSS_200718.RData")
TSSDat <- datByOligoSpread
rm(datByOligoSpread)

#TSS needs some extra wrangling
TSSDat <- TSSDat %>% select(-contains("His"))
TSSDat <- TSSDat %>% mutate('2016_A' = (TSSDat$"2016_A1" + TSSDat$"2016_A2")/2)
TSSDat <- TSSDat %>% mutate('2016_B' = (TSSDat$"2016_B1" + TSSDat$"2016_B2")/2)
TSSDat <- TSSDat %>% select(-c("2016_A1", "2016_A2", "2016_B1", "2016_B2"))

# find the oligos that have the controls
oligoDesign <- read.table("jointDesignFinal27k_150612_ReverseComplementWhenALargerT.txt", sep="\t", head=TRUE, stringsAsFactors=FALSE)
rownames(oligoDesign) <- oligoDesign$oligoID
# Eran controls:
EranControlsUpstream <- oligoDesign[str_detect(oligoDesign$oligoID, "singleVariantTiles") & oligoDesign$oligoBlock == "EranControl",]
EranControlsTSS <- oligoDesign[str_detect(oligoDesign$oligoID, "firstTileUTR") & oligoDesign$oligoBlock == "EranControl",]

# reduce each dataset to the shared controls
UpstreamDat <- UpstreamDat[UpstreamDat$oligoID %in% rownames(EranControlsUpstream),]
TSSDat <- TSSDat[TSSDat$oligoID %in% rownames(EranControlsTSS),]

# rename/add control name so we can merge on them
UpstreamDat$EranID <- oligoDesign[UpstreamDat$oligoID, "variantIDs"]
TSSDat$EranID <- oligoDesign[TSSDat$oligoID, "variantIDs"]


# now we should be able to simply merge them, and only the Sharons should be there in both
combinedDat <- merge(TSSDat, UpstreamDat, by.x="EranID", by.y="EranID", suffixes=c("_TSS", "_Upstream"))

# now let's correlate these all:
plotCor <- cor(combinedDat[,c("2016_1", "2016_2", "2016_3", "2016_4", "2016_5", "2016_6", "2016_A", "2016_B", "2017_A", "2017_B", "2018_A_TSS", "2018_B_TSS", "2018_1", "2018_2", "2018_3", "2018_4", "2018_A_Upstream", "2018_B_Upstream")], use="complete", method="s")

pdf("sharonCorrelations.pdf", height=14, width=14)
corrplot.mixed(plotCor,
    upper="number",
    lower="ellipse"
)
dev.off()

# distributions
#TSS
summary(c(plotCor[c("2016_1", "2016_2", "2016_3", "2016_4", "2016_5", "2016_6", "2016_A", "2016_B", "2017_A", "2017_B", "2018_A_TSS", "2018_B_TSS"),c("2016_1", "2016_2", "2016_3", "2016_4", "2016_5", "2016_6", "2016_A", "2016_B", "2017_A", "2017_B", "2018_A_TSS", "2018_B_TSS")]))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.5734  0.7116  0.8015  0.7972  0.8770  1.0000

#Upstream
summary(c(plotCor[c("2018_1", "2018_2", "2018_3", "2018_4", "2018_A_Upstream", "2018_B_Upstream"),c("2018_1", "2018_2", "2018_3", "2018_4", "2018_A_Upstream", "2018_B_Upstream")]))
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.9558  0.9816  0.9897  0.9873  0.9970  1.0000

#between
summary(c(plotCor[c("2016_1", "2016_2", "2016_3", "2016_4", "2016_5", "2016_6", "2016_A", "2016_B", "2017_A", "2017_B", "2018_A_TSS", "2018_B_TSS"),c("2018_1", "2018_2", "2018_3", "2018_4", "2018_A_Upstream", "2018_B_Upstream")]))
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.5939  0.6756  0.7637  0.7438  0.7953  0.8538



