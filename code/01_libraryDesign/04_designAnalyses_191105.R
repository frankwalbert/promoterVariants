# extract descriptions about the design

library(tidyverse)

jointDesignFinal27k <- read.table("jointDesignFinal27k_150612_ReverseComplementWhenALargerT.txt", head=TRUE, sep="\t", stringsAsFactors=FALSE)

#allVariants <- unlist(lapply(jointDesignFinal27k$variantIDs[jointDesignFinal27k$designType == "firstTileUTR"], function(x){
#allVariants <- unlist(lapply(jointDesignFinal27k$variantIDs[jointDesignFinal27k$designType == "singleVariantTiles"], function(x){
allVariants <- unlist(lapply(jointDesignFinal27k$variantIDs[jointDesignFinal27k$designType %in% c("firstTileUTR", "singleVariantTiles")], function(x){
    theseVariants <- strsplit(x, ";")[[1]]
}))

# remove controls
allVariantsNoControl <- allVariants[!grepl("sharon", allVariants)]
allVariantsNoControl <- allVariantsNoControl[!grepl("hunter", allVariantsNoControl)]

theseVariantsSplit <- t(sapply(unique(allVariantsNoControl), function(y){strsplit(y, "_")[[1]]}))

# was for diagnoses of controls:
#theseVariantsSplit[[which(sapply(theseVariantsSplit, length) == 3)]]

dim(theseVariantsSplit)
# TSS: 3,645
# UpStream: 4547
# both: 7,005


# also count genes:
allGenes <- unique(str_split_fixed(jointDesignFinal27k$gene[jointDesignFinal27k$designType %in% c("firstTileUTR")], "_", n=2)[,1])
allGenes <- unique(str_split_fixed(jointDesignFinal27k$gene[jointDesignFinal27k$designType %in% c("singleVariantTiles")], "_", n=2)[,1])
#allGenes <- unique(str_split_fixed(jointDesignFinal27k$gene[jointDesignFinal27k$designType %in% c("firstTileUTR", "singleVariantTiles")], "_", n=2)[,1])
length(unique(allGenes[allGenes != "control"]))
# TSS: 2172
# UpStream: 1,918
# both: 3,076

# how many SNPs?
summary(nchar(theseVariantsSplit[,3]) == nchar(theseVariantsSplit[,4]))
# FALSE 1247 TRUE 5758

summary(nchar(theseVariantsSplit[,3]) == nchar(theseVariantsSplit[,4]) & nchar(theseVariantsSplit[,4]) == 1)
# only 5741 true:

theseVariantsSplit[(nchar(theseVariantsSplit[,3]) == nchar(theseVariantsSplit[,4]) & nchar(theseVariantsSplit[,4]) != 1),]
#chrI_192551_CAAA_AAAA                                        "chrI"   "192551"  "CAAA"                    "AAAA"
#chrII_84412_AT_TT                                            "chrII"  "84412"   "AT"                      "TT"
#chrIV_1517934_CG_AG                                          "chrIV"  "1517934" "CG"                      "AG"
#chrVII_459753_GT_TT                                          "chrVII" "459753"  "GT"                      "TT"
#chrVII_1051846_CTT_TTT                                       "chrVII" "1051846" "CTT"                     "TTT"
#chrX_450783_GT_TT                                            "chrX"   "450783"  "GT"                      "TT"
#chrXI_557499_GA_AA                                           "chrXI"  "557499"  "GA"                      "AA"
#chrXII_127407_AT_TT                                          "chrXII" "127407"  "AT"                      "TT"
#chrXIV_412586_AT_GT                                          "chrXIV" "412586"  "AT"                      "GT"
#chrXVI_338769_TAAG_CAAG                                      "chrXVI" "338769"  "TAAG"                    "CAAG"
#chrIV_370968_TTATATATA_ATATATATA                             "chrIV"  "370968"  "TTATATATA"               "ATATATATA"
#chrV_266246_CT_AT                                            "chrV"   "266246"  "CT"                      "AT"
#chrVII_115467_CGA_TGA                                        "chrVII" "115467"  "CGA"                     "TGA"
#chrXII_514965_CTT_TTT                                        "chrXII" "514965"  "CTT"                     "TTT"
#chrI_202713_TA_CA                                            "chrI"   "202713"  "TA"                      "CA"
#chrII_353616_ATATATATATATATATATATGTG_GTATATATATATATATATATGTG "chrII"  "353616"  "ATATATATATATATATATATGTG" "GTATATATATATATATATATGTG"
#chrXV_354680_ACT_TCT                                         "chrXV"  "354680"  "ACT"                     "TCT"
# these are SNPs?

# size spectrum of indels:
indels <- theseVariantsSplit[(nchar(theseVariantsSplit[,3]) != nchar(theseVariantsSplit[,4])),]
indelSizes <- nchar(indels[,3]) - nchar(indels[,4])
summary(abs(indelSizes))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#1.000   1.000   1.000   2.655   2.000 132.000

# most are small:
table(abs(indelSizes))
#  1   2   3   4   5   6   7   8   9  10  11  12  14  15  16  17  18  19  20  21  24  28  31  32  33  37  38  40  42  44  52 132
#770 226  62  58  30  21   6  11   5  11   1   9   5   4   3   1   3   1   2   4   1   1   1   1   3   1   1   1   1   1   1   1

pdf("indelSizes.pdf", width=12, height=6)
barplot(table(indelSizes))
dev.off()


######
# distribution of genes and their variants:
variants <- jointDesignFinal27k$variantIDs[jointDesignFinal27k$designType %in% c("firstTileUTR", "singleVariantTiles")]
genes <- str_split_fixed(jointDesignFinal27k$gene[jointDesignFinal27k$designType %in% c("firstTileUTR", "singleVariantTiles")], "_", n=2)[,1]

genesAndVariants <- data.frame(gene = genes, variants, stringsAsFactors=FALSE)

# remove controls
genesAndVariants <- genesAndVariants[!grepl("sharon", genesAndVariants$variants),]
genesAndVariants <- genesAndVariants[!grepl("hunter", genesAndVariants$variants),]
# now 16677 oligos
# before removing variants, we were at 17093 (the sum of the two numbers in the paper)
# the 3076 genes in the paper are WITHOUT the two sets of controls

genesVariantsUnfolded <- genesAndVariants %>%
    unnest(variant = str_split(variants, ";")) %>%
    select(-variants) %>%
    distinct()
    
genesVariantsFolded <- genesVariantsUnfolded %>%
    group_by(gene) %>%
    summarize(numberVariants = n())

summary(genesVariantsFolded$numberVariants)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#1.000   1.000   2.000   2.532   3.000  38.000

# how many genes had one variant?
length(which(genesVariantsFolded$numberVariants == 1))
# 1391 (out of 3076; 45%)

pdf("variantsPerGene_TSS_UpStream.pdf", width=10, height=5)
par(mfrow=c(1, 2))
hist(genesVariantsFolded$numberVariants, breaks=50, main="Entire range", xlab="Tested variants per gene")
barplot(c(as.numeric(table(genesVariantsFolded$numberVariants[genesVariantsFolded$numberVariants <= 10])), length(which(genesVariantsFolded$numberVariants > 10))), main="At most 10 variants", xlab="Tested variants per gene")
dev.off()





##################
# variants in the genome (as opposed to the design)

#allVariants = read.table("A_forFrank051815.sort.blocked.vcf", stringsAsFactors=FALSE)
# 45543 variants in here

# now putting it through ENSEMBLv98 VEP on 11/7/19
# ran this with output "one specific consequence", if select "most severe", I don't get the gene etc info back
# restricted upstream/downstream to 1kb while running VEP
allVariants = read.table("A_forFrank051815.sort.blocked_VEP_191107_mod.txt", stringsAsFactors=FALSE, sep="\t", head=TRUE)

table(allVariants$Consequence)
length(which(allVariants$Consequence == "upstream_gene_variant"))
# 14210
alleles <- str_split_fixed(str_split_fixed(allVariants$Uploaded_variation, "_", n=2)[,2], "/", n=2)

table(abs(nchar(alleles[,1]) - nchar(alleles[,2])))
summary(abs(nchar(alleles[,1]) - nchar(alleles[,2])) == 0)
#   Mode   FALSE    TRUE
#logical    3936   41607 <- this many SNPs

# restrict to upstream:
summary(abs(nchar(alleles[allVariants$Consequence == "upstream_gene_variant",1]) - nchar(alleles[allVariants$Consequence == "upstream_gene_variant",2])) == 0)
#   Mode   FALSE    TRUE
#logical    2442   11768


# distribution of THESE (i.e. all, irrespective of whether they are in the design) variants per gene promoter:
# need gene info to get a list of all genes
allGenes <- read.table("ensemblGenes_ensembl83_160307_MOD.txt", stringsAsFactors=FALSE, sep="\t", head=TRUE)
allGenes <- allGenes[,1] 
upstreamVars <- allVariants[allVariants$Consequence == "upstream_gene_variant",]
length(unique(upstreamVars$Gene))
#3176
# how many genes have no annotated upstream variants?
length(which(!allGenes %in% upstreamVars$Gene))
#3950(!)
head(allGenes[which(!allGenes %in% upstreamVars$Gene)])

# how many genes have no variants at all? (presumbly, many of these will be in the "blocked" regions)
length(which(!allGenes %in% allVariants$Gene))
# 1914

length(unique(allVariants$Gene))
#5,214
