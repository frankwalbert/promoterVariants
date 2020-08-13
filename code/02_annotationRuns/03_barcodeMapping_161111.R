# also see the corresponding TSS code, which may be slightly better commented
# PEAR the reads:
# run from within the read files:
~/src/pear/bin/pear -f read1.fastq -r read2.fastq -v 75 -m 250 -j 20 -y 2G -c 40 -b 64 -o pear_min75_max250

# we have to juggle this because the unzipped read files are so huge (up to 90Gb each!) and we blow out our space limit
# keep one lane zipped and process the other

Assembled reads ...................: 132,654,107 / 155,274,356 (85.432%)
Discarded reads ...................: 0 / 155,274,356 (0.000%)
Not assembled reads ...............: 22,620,249 / 155,274,356 (14.568%)
# i.e., quite good


####################################
# R starts here


library(ShortRead)
library(plyr)
library(dplyr)
library(doSNOW) # for parallel ddply
library(gplots)
library(RColorBrewer)
library(Biostrings)
library(ggplot2)
library(devtools) # for for with_options to reformat allele strings in the design
library(stringr) # for for with_options to reformat allele strings in the design
library(stringdist) # to compute barcode distances

# this Upstream annotation run was lane 2 – lane 1 was an "AUG" library we did not take forward
lane = 2

allMergedReads = FastqStreamer(paste0("../../reads/SxaQSEQsXap126L", lane, "/pear_min75_max250.assembled.fastq.gz", sep=""))

invariantStretches = c("CCTGCAGGGCCGTGTGAAGCTGG", "CCTGCAGGGTTCCGCAGCCACAC")
# add TSS here so that can make consistent analyses at the bottom, when reading in those files
libraryID = c("AUG", "UpStream", "TSS")

outFile = paste0("allOligosBarcodes_lane", lane, ".txt")
#write.table("", file=outFile, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
file.create(outFile)
repeat {
    thisChunk <- sread(yield(allMergedReads))
    #print(thisChunk)
    if (length(thisChunk) == 0) { break }

    invariantStretch1 = vmatchPattern(invariantStretches[lane], thisChunk, max.mismatch=1)
    invariantStretch2 = vmatchPattern("GGCCGTAATGGCC", thisChunk, max.mismatch=1)
    #print(invariantStretch1)
    #print(invariantStretch2)
    barcodes = substr(thisChunk, 1, as.numeric(as.character(startIndex(invariantStretch1))) - 1)
    oligos = substr(thisChunk, as.numeric(as.character(endIndex(invariantStretch1))) + 1, as.numeric(as.character(startIndex(invariantStretch2))) - 1)

    barcodesOligos = cbind(barcodes, oligos)
    barcodesOligos = barcodesOligos[complete.cases(barcodesOligos),]
    barcodesOligos = barcodesOligos[nchar(barcodesOligos[,1]) == 20,]

    # now, what do do with them? write them to a big file, concatenate?
    # then operate on that file
    write.table(barcodesOligos, file=outFile, append=TRUE, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
}

# the resulting file is 24Gb and has 170,317,989 (~170M) barcode-oligo combinations


barcodesOligos = read.table(outFile, sep="\t")
barcodeOligoCounts = plyr::count(data.frame(barcodesOligos))
sum(barcodeOligoCounts[,3]) == nrow(barcodesOligos)
# TRUE
# numbers of reads:
# UpStream: 124,558,155

barcodeOligoCounts = barcodeOligoCounts[order(barcodeOligoCounts[,"freq"], decreasing=TRUE),]
names(barcodeOligoCounts) = c("barcodes", "oligos", "freq")

# FIX THE  FACTORS
barcodeOligoCounts[,1] <- as.character(barcodeOligoCounts[,1])
barcodeOligoCounts[,2] <- as.character(barcodeOligoCounts[,2])

#save(barcodeOligoCounts, file=paste0("R_barcodeOligoCounts_", libraryID[lane], "_161202.RData"))
# this file has one row per barcode/oligo combo and the count of the combo
# has all detected oligos, knows nothing yet aout the design or oligos going to the same barcode

nrow(barcodeOligoCounts)
median(barcodeOligoCounts[,3])
length(which(barcodeOligoCounts[,3] > 10))
length(which(barcodeOligoCounts[,3] == 1))
length(unique(barcodeOligoCounts[,2]))
length(unique(barcodeOligoCounts[,1]))

# UpStream:
# 53,991,723 rows. good.
# most frequent barcode/oligo: 1306. good (as in, nothing took over)
# median count: 1
# 778,720 combos more than 10
# 29,542,846 seen once.
# unique oligos: 14,600,594
# unique barcodes: 46,349,032


pdf(paste0("barcodeHistogramRaw_", libraryID[lane], ".pdf", sep=""), width=9, height=9)
par(mfrow=c(2,2))
hist(barcodeOligoCounts[,3], breaks=100, xlab="count of barcode/oligo combination", main="all")
hist(barcodeOligoCounts[,3][barcodeOligoCounts[,3] < 100], breaks=100, xlab="count of barcode/oligo combination", main="< 100")
hist(barcodeOligoCounts[,3][barcodeOligoCounts[,3] < 100 & barcodeOligoCounts[,3] > 1], breaks=100, xlab="count of barcode/oligo combination", main="1 < x < 100")
hist(barcodeOligoCounts[,3][barcodeOligoCounts[,3] < 100 & barcodeOligoCounts[,3] > 10], breaks=100, xlab="count of barcode/oligo combination", main="10 < x < 100")
dev.off()

# for every barcode, count how many oligos they tag, the distribution of these oligos, and their sequences
bcCounts = ddply(barcodeOligoCounts, .(barcodeOligoCounts$barcodes), .fun = function(x) {
    c(sum(x[,3]), nrow(x), paste(x[,3], collapse=","), paste(x[,2], collapse=","))
}, .progress="text", .parallel=FALSE)

bcCounts[,1] <- as.character(bcCounts[,1])
bcCounts[,2] <- as.integer(bcCounts[,2])
bcCounts[,3] <- as.integer(bcCounts[,3])
bcCounts = bcCounts[order(bcCounts[,2], decreasing=TRUE),]
colnames(bcCounts) = c("bcSeq", "count", "oligoCount", "oligoDistribution", "oligoSeqs")

# aggregate oligos tagged by a given barcode
# this get pretty intense: ran for ~24 to ~36h and took up to 77Gb RAM
doFunc <- function(x){
    ret = data.frame(count=sum(x[,3]), oligoCount=nrow(x), oligoDistribution=paste(x[,3], collapse=","), oligoSeqs=paste(x[,2], collapse=","), fractionCountsPerTopOligo=max(x[,3])/sum(x[,3]), stringsAsFactors=FALSE)
    ret
}

bcCountsDplyr = barcodeOligoCounts %>% group_by(barcodes) %>% do(doFunc(.))
bcCountsDplyr <- bcCountsDplyr[order(bcCountsDplyr$count, decreasing=TRUE),]
#save(bcCountsDplyr, file=paste0("R_bcCountsDplyr_", libraryID[lane], "_161202.RData", sep=""))

nrow(bcCountsDplyr)
# UpStream: 46,349,032, i.e. the expected number of unique barcodes. Good.

# some format fixes, then save
bcCounts <- data.frame(bcCountsDplyr)
bcCounts[,4] <- gsub("[^0-9,]", '', bcCounts[,4])
bcCounts[,5] <- gsub("[^ACGT,]", '', bcCounts[,5])
#save(bcCounts, file=paste0("R_bcCounts_", libraryID[lane], "_161203.RData", sep=""))

# examine the top few by eye:
# UpStream:
# the first barcodes are all GGGGGGG with some other bases strewn in here and there. These have terrible "fraction best oligo" metrics


# plot
pdf(paste0("barcodeVsOligoCounts_", libraryID[lane], ".pdf", sep=""), width=15, height=5)
par(mfrow=c(1,3))
# instead of random, do top 1k
plotSample = c(1:5000)
plot(bcCounts$count[plotSample], bcCounts$oligoCount[plotSample], cex=.5, col="#00000022", xlab="barcode count", ylab="number oligos tagged", main="most frequent barcodes")
abline(0, 1, col="grey")
plotSample = sample(1:length(bcCounts$count), 5e3)
plot(bcCounts$count[plotSample], bcCounts$oligoCount[plotSample], cex=.5, col="#00000022", xlab="barcode count", ylab="number oligos tagged", main="random barcodes")
abline(0, 1, col="grey")
smoothScatter(bcCounts$count, bcCounts$oligoCount, cex=.5, col="#00000022", xlab="barcode count", ylab="number oligos tagged", main="all barcodes")
abline(0, 1, col="grey")
dev.off()

fractionCountsPerOligo = t(sapply(bcCounts$oligoDistribution, function(x){
    ret = rep(0, 11)
    thisRes = as.numeric(strsplit(x, ",")[[1]])
    thisRes = thisRes/sum(thisRes)
    thisLen = length(thisRes)
    ret[1:min(c(10, thisLen))] <- thisRes[1:min(c(10, thisLen))]
    if (thisLen > 10){
        ret[11] <- sum(thisRes[11:length(thisRes)])
    }
    ret
}, USE.NAMES=FALSE))
#save(fractionCountsPerOligo, file=paste0("R_fractionCountsPerOligo_", libraryID[lane], ".RData"))

# what is the distribution of mass in the top oligo?
pdf(paste0("oligoCountsPerBarcodeHist_", libraryID[lane], ".pdf"), width=12, height=4)
par(mfrow=c(1,3))
hist(fractionCountsPerOligo[,1], breaks=50, main="fraction counts that are top oligo")
hist(fractionCountsPerOligo[,2], breaks=50, xlim=c(0,1), main="fraction counts that are 2nd oligo")
hist(rowSums(fractionCountsPerOligo[,3:11]), breaks=50, xlim=c(0,1), main="fraction counts that are 3rd and higher oligo")
dev.off()

pdf(paste0("oligoCountsPerBarcodeHist_atLeast10Counts_", libraryID[lane], ".pdf"), width=12, height=4)
par(mfrow=c(1,3))
hist(fractionCountsPerOligo[bcCounts$count > 9, 1], breaks=50, main="fraction counts that are top oligo")
hist(fractionCountsPerOligo[bcCounts$count > 9, 2], breaks=50, xlim=c(0,1), main="fraction counts that are 2nd oligo")
hist(rowSums(fractionCountsPerOligo[bcCounts$count > 9, 3:11]), breaks=50, xlim=c(0,1), main="fraction counts that are 3rd and higher oligo")
dev.off()
# overall, tagging of multiple oligos by a barcode is not a big problem



######################################
# map to design

readOligosRaw = read.table("jointDesignFinal27k_150612_ReverseComplementWhenALargerT.txt", head=TRUE, stringsAsFactors=FALSE)
libraryDesignNames = c("firstTileAUG", "singleVariantTiles", "firstTileUTR")

rownames(readOligosRaw) = readOligosRaw$oligoID
readOligosRaw = readOligosRaw[readOligosRaw$designType==libraryDesignNames[lane],]
# note that we have revcomped some of the oligos for the array design (to avoid A stretches)
# these do NOT need to revcomped now - they are already sense compared to the annotation run
# revcomp the others
readOligos = readOligosRaw
readOligos$sequence[!readOligos$reverseComplementFlag] = as.character(reverseComplement(DNAStringSet(readOligos$sequence[!readOligos$reverseComplementFlag])))

# adjust the designed oligos: cut off PS2 but NOT the AscI site
PS2s <- c("GCCGTGTGAAGCTGG", "GTTCCGCAGCCACAC", "GGTTTAGCCGGCGTG")

readOligos$sequence <- substr(readOligos$sequence, regexpr(PS2s[lane], readOligos$sequence) + nchar(PS2s[lane]), regexpr("GGCCGTAATGGCC", readOligos$sequence)-1)

oligoDesign = readOligos$sequence
names(oligoDesign) = readOligos$oligoID

# recode the allele strings
readOligosTemp = readOligos
ubo <- unique(readOligosTemp$oligoBlock)
for(i in 1:length(ubo)){
    print(i/length(ubo))
    x = ubo[i]
    aSIn = as.numeric(readOligosTemp[readOligosTemp$oligoBlock == x,]$alleleString)
    aSOut = with_options(c(scipen = 999), str_pad(aSIn, max(nchar(aSIn)), pad=0))
    # order is important!
    aSOut = gsub("1", "2", aSOut)
    aSOut = gsub("0", "1", aSOut)
    #print(aSOut)
    readOligosTemp[readOligosTemp$oligoBlock == x,]$alleleString <- aSOut
}
readOligosTemp$alleleString[readOligosTemp$alleleString == "2222221111111111"] <- "2222222222222222"
readOligosTemp$alleleString[readOligosTemp$alleleString == "22222211111111"] <- "22222222222222"
unique(readOligosTemp$alleleString)

readOligos <- readOligosTemp
#save(readOligos, file=paste0("R_readOligos_", libraryID[lane],".RData", sep=""))


# process all the barcodes, no culling
# pull out the respective first oligo sequence
temp = sapply(bcCounts[,5], function(x){strsplit(x, ",")[[1]][1]})
names(temp) = NULL
bcCountsTopOligoOnly = bcCounts
bcCountsTopOligoOnly[,5] = temp
colnames(bcCountsTopOligoOnly) = c("bcSeq", "count", "oligoCount", "oligoCountDistribution", "oligoSeq", "fractionTopOligo")
save(bcCountsTopOligoOnly, file=paste0("R_bcCountsTopOligoOnly_161203_", libraryID[lane],".RData", sep=""))

# for GEO: write to file
#write.table(bcCountsTopOligoOnly, file=paste0("R_bcCountsTopOligoOnly_161203_", libraryID[lane],".txt", sep=""), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


#############
# now find the oligo from the design which matches the observed sequences
bcCountsAssigned = merge(bcCountsTopOligoOnly, readOligos, by.x = "oligoSeq", by.y = "sequence")
save(bcCountsAssigned, file = paste0("R_bcCountsAssigned_", libraryID[lane],".RData", sep=""))


# some stats
nrow(bcCountsAssigned)
# upstream:
# 20,008,691

length(unique(bcCountsAssigned$bcSeq))
# UpStream 19,983,696

length(unique(bcCountsAssigned$oligoSeq))
# UpStream
# 9636 oligos, 97.5% of 9,882 designed

# count for each oligo how many barcodes it has, and how often those barcodes were seen
# using plyr (need to add the zeros after this)
assignedOligosBarcodeCounts = ddply(bcCountsAssigned, .(bcCountsAssigned$oligoID), .fun = function(x) {
    nrow(x)
})
assignedOligosBarcodeSums = ddply(bcCountsAssigned, .(bcCountsAssigned$oligoID), .fun = function(x) {
    sum(x[,"count"])
})
rownames(assignedOligosBarcodeSums) = assignedOligosBarcodeSums[,1]
rownames(assignedOligosBarcodeCounts) = assignedOligosBarcodeCounts[,1]

# now drop the ones with counts into those with none (fill)
assignedOligosBarcodeCountsWithZeros = rep(0, length(oligoDesign))
names(assignedOligosBarcodeCountsWithZeros) = names(oligoDesign)
assignedOligosBarcodeSumsWithZeros = assignedOligosBarcodeCountsWithZeros

assignedOligosBarcodeCountsWithZeros[names(assignedOligosBarcodeCountsWithZeros) %in% rownames(assignedOligosBarcodeCounts)] <- assignedOligosBarcodeCounts[names(assignedOligosBarcodeCountsWithZeros)[names(assignedOligosBarcodeCountsWithZeros) %in% rownames(assignedOligosBarcodeCounts)], 2]
assignedOligosBarcodeSumsWithZeros[names(assignedOligosBarcodeSumsWithZeros) %in% rownames(assignedOligosBarcodeCounts)] <- assignedOligosBarcodeSums[names(assignedOligosBarcodeSumsWithZeros)[names(assignedOligosBarcodeSumsWithZeros) %in% rownames(assignedOligosBarcodeCounts)], 2]

pdf(paste0("oligoCountsAssigned_", libraryID[lane], ".pdf", sep=""), width=12, height=8)
par(mfrow=c(2,3))

hist(assignedOligosBarcodeCountsWithZeros, breaks=50, xlab="number of barcodes per oligo", main="all oligos")
hist(assignedOligosBarcodeCountsWithZeros[assignedOligosBarcodeCountsWithZeros <= 5000], breaks=50, xlab="number of barcodes per oligo", main="oligos with barcode number < 5000")
hist(log10(assignedOligosBarcodeCountsWithZeros + 1), breaks=50, xlab="log10(number of barcodes per oligo)", main="all oligos")

hist(assignedOligosBarcodeCounts[,2], breaks=50, xlab="number of barcodes per oligo", main="oligos with at least 1 barcode")
hist(assignedOligosBarcodeCounts[,2][assignedOligosBarcodeCounts[,2] <= 5000], breaks=50, xlab="number of barcodes per oligo", main="oligos with barcode number 0 < x < 5000")
hist(log10(assignedOligosBarcodeCounts[,2] + 1), breaks=50, xlab="log10(number of barcodes per oligo)", main="oligos with at least 1 barcode")

dev.off()

summary(assignedOligosBarcodeCountsWithZeros)
# upstream
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0     407     972    2025    2215   60000

# drop the oligos that aren't seen:
summary(assignedOligosBarcodeCounts[,2])
# upStream
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#1     442    1008    2074    2267   60000



# USED IN PAPER
# distribution of just the assigned barcodes:
pdf("barcodeHistogram_bcCountsTopOligoOnly.pdf", width=9, height=9)
par(mfrow=c(2,2))
hist(bcCountsTopOligoOnly$count, breaks=100, xlab="count of barcode/oligo combination", main="all")
hist(bcCountsTopOligoOnly$count[bcCountsTopOligoOnly$count < 50], breaks=100, xlab="count of barcode/oligo combination", main="< 50")
hist(bcCountsTopOligoOnly$count[bcCountsTopOligoOnly$count < 50 & bcCountsTopOligoOnly$count > 1], breaks=100, xlab="count of barcode/oligo combination", main="1 < x < 50")
hist(bcCountsTopOligoOnly$count[bcCountsTopOligoOnly$count < 50 & bcCountsTopOligoOnly$count > 10], breaks=100, xlab="count of barcode/oligo combination", main="10 < x < 50")
dev.off()

summary(bcCountsTopOligoOnly$count)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#1.000    1.000    2.000    2.687    3.000 2587.000



###################################
# examine oligo dropouts again

designedInsets = readOligosRaw$sequence
# flip them into the same space so that they can be parsed the same way
designedInsets[!readOligosRaw$reverseComplementFlag] = as.character(reverseComplement(DNAStringSet(designedInsets[!readOligosRaw$reverseComplementFlag])))
# cut out the inset from the invariant seqs
designedInsets <- substr(designedInsets, start = nchar("GGTTTAGCCGGCGTGGGCGCGCC")+1, stop = regexpr("GGCCGTAATGGCCCACCAAGGGCGATCG", designedInsets) - 1)
# flip them all back into design space
# DON'T DO FOR POSITIONAL ANALYSES (which are not included in this round so far)
#designedInsets[!readOligosRaw$reverseComplementFlag] = as.character(reverseComplement(DNAStringSet(designedInsets[!readOligosRaw$reverseComplementFlag])))
names(designedInsets) = rownames(readOligosRaw)

# now feed these (NOT the whole oligos; these are anyway invariant - or only come in two versions at least; they cannot explain the dropouts)
zeroOligos = DNAStringSet(designedInsets[names(assignedOligosBarcodeCountsWithZeros[assignedOligosBarcodeCountsWithZeros == 0])])
nonZeroOligos = DNAStringSet(designedInsets[names(assignedOligosBarcodeCountsWithZeros[assignedOligosBarcodeCountsWithZeros > 0])])
allOligos = DNAStringSet(designedInsets)

# plot dinucleotide frequencies at the first two bases vs oligo count
# make sure the names are aligned, although they probably are going in
oligoCountsAndFirstTwoBases = data.frame(substr(allOligos, 1, 2), assignedOligosBarcodeCountsWithZeros[names(allOligos)], readOligosRaw[names(allOligos),]$reverseComplementFlag)
names(oligoCountsAndFirstTwoBases) = c("dinucleotide", "count", "rcFlag")
# ugly hack to force the plot to show the empty factor levels (two DNs are only seen in forward or reverse)
tempPadder = data.frame(c("CG", "CG", "TG", "TG"), rep(NA, 4), c(TRUE, FALSE, TRUE, FALSE))
names(tempPadder) = names(oligoCountsAndFirstTwoBases)
oligoCountsAndFirstTwoBases = rbind(oligoCountsAndFirstTwoBases, tempPadder)
names(oligoCountsAndFirstTwoBases) = c("dinucleotide", "count", "rcFlag")


pdf("dinucleotidesVsCount.pdf", width=8, height=5)
p <- ggplot(oligoCountsAndFirstTwoBases[oligoCountsAndFirstTwoBases$rcFlag,], aes(factor(dinucleotide), count))
p + geom_boxplot() + geom_jitter(position = position_jitter(width = .4, height=0), alpha=0.5, size=0.5) + ggtitle("reverse synthesis")

p <- ggplot(oligoCountsAndFirstTwoBases[!oligoCountsAndFirstTwoBases$rcFlag,], aes(factor(dinucleotide), count))
p + geom_boxplot() + geom_jitter(position = position_jitter(width = .4, height=0), alpha=0.5, size=0.5) + ggtitle("forward synthesis")

dev.off()

# dropout due to G(G) at the first dinlucleotide
