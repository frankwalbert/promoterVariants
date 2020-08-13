
# Prior to running this R code:
# run PEAR on paired end data
# concatenate the output from the two lanes for this TSS annotation

library(ShortRead)
library(plyr)
library(dplyr)
library(doSNOW) # for parallel ddply
library(gplots)
library(RColorBrewer)
library(Biostrings)
library(ggplot2)
library(devtools) # for with_options to reformat allele strings in the design
library(stringr) # for with_options to reformat allele strings in the design

# would need to generate this from GEO; too large to be shared here
# PEAR the reads, then cat the two TSS annotation lanes
allMergedReads = FastqStreamer("allAssembledReads.fastq.gz")

# in here, 'yield' gives 1M reads (by default) in order of the file
# can then process them in chunks
# also note that the streamer works on gz files
outFile = "allOligosBarcodes.txt"
#write.table("", file=outFile, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
file.create(outFile)
repeat {
    thisChunk <- sread(yield(allMergedReads))
    if (length(thisChunk) == 0) { break }

    invariantStretch1 = vmatchPattern("CCTGCAGGGGTTTAGCCGGCGTG", thisChunk, max.mismatch=1)
    invariantStretch2 = vmatchPattern("GGCCGTAATGGCC", thisChunk, max.mismatch=1)
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

barcodeOligoCounts = barcodeOligoCounts[order(barcodeOligoCounts[,"freq"], decreasing=TRUE),]
names(barcodeOligoCounts) = c("barcodes", "oligos", "freq")

# FIX THE FACTORS
barcodeOligoCounts[,1] <- as.character(barcodeOligoCounts[,1])
barcodeOligoCounts[,2] <- as.character(barcodeOligoCounts[,2])
#save(barcodeOligoCounts, file="R_barcodeOligoCounts_160627")

# at this point, we have made a table that simply lists a barcode sequence, the oligo sequence it tags, and how often this exact pair was seen
# it does not yet consider the same barcode tagging different unique sequences


pdf("barcodeHistogramRaw.pdf", width=9, height=9)
par(mfrow=c(2,2))
hist(barcodeOligoCounts[,3], breaks=100, xlab="count of barcode/oligo combination", main="all")
hist(barcodeOligoCounts[,3][barcodeOligoCounts[,3] < 100], breaks=100, xlab="count of barcode/oligo combination", main="< 100")
hist(barcodeOligoCounts[,3][barcodeOligoCounts[,3] < 100 & barcodeOligoCounts[,3] > 1], breaks=100, xlab="count of barcode/oligo combination", main="1 < x < 100")
hist(barcodeOligoCounts[,3][barcodeOligoCounts[,3] < 100 & barcodeOligoCounts[,3] > 10], breaks=100, xlab="count of barcode/oligo combination", main="10 < x < 100")
dev.off()



################################
# for every barcode, count how many oligos, their distribution, and their seqs
# this step gets very intense
# for this TSS annotation, it had been run on MSI using "02_TSS_oligoCounterSingleThreadDplyr.R"
# see that file for the code
# that code will produce the file we now load:

load("R_bcCounts_singleThread_160627")
# resulting object has 23,648,568 rows
# i.e. the expected number of unique barcodes
# good
# the oligo count and sequence strings got messed up:
# they look like "c(x, x, x, x)"
# get rid of the c(), and spaces:

bcCounts <- data.frame(bcCounts)
bcCounts[,4] <- gsub("[^0-9,]", '', bcCounts[,4])
bcCounts[,5] <- gsub("[^ACGT,]", '', bcCounts[,5])
#save(bcCounts, file="R_bcCounts_singleThread_160627_correctFormat")

# plot in R
pdf("barcodeVsOligoCounts.pdf", width=15, height=5)
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

# examine distribution of barcode/oligo counts
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
#save(fractionCountsPerOligo, file="R_fractionCountsPerOligo")


# what is the distribution of mass in the top oligo?
pdf("oligoCountsPerBarcodeHist.pdf", width=12, height=4)
par(mfrow=c(1,3))
hist(fractionCountsPerOligo[,1], breaks=50, main="fraction counts that are top oligo")
hist(fractionCountsPerOligo[,2], breaks=50, xlim=c(0,1), main="fraction counts that are 2nd oligo")
hist(rowSums(fractionCountsPerOligo[,3:11]), breaks=50, xlim=c(0,1), main="fraction counts that are 3rd and higher oligo")
dev.off()
# this looks great
# tagging of multiple oligos by a barcode is not a problem



######################################
# map to design

# I'm using the file that used to be this one from kserver:
# moved from its parallel home on my dropbox code folder
readOligosRaw = read.table("jointDesignFinal27k_150612_ReverseComplementWhenALargerT.txt", head=TRUE, stringsAsFactors=FALSE)

rownames(readOligosRaw) = readOligosRaw$oligoID
readOligosRaw = readOligosRaw[readOligosRaw$designType=="firstTileUTR",]
# note that we have revcomped some of the oligos for the array design (to avoid A stretches)
# these do NOT need to revcomped now - they are already sense compared to the annotation run
# revcomp the others
readOligos = readOligosRaw
readOligos$sequence[!readOligos$reverseComplementFlag] = as.character(reverseComplement(DNAStringSet(readOligos$sequence[!readOligos$reverseComplementFlag])))

# adjust the designed oligos: cut off PS2 but NOT the AscI site
# PS2(GGTTTAGCCGGCGTG for the UTR library; revcomp) - AscI(GGCGCGCC) - oligo - SfiIA(GGCCGTAATGGCC in revcomp)
readOligos$sequence <- substr(readOligos$sequence, regexpr("GGTTTAGCCGGCGTG", readOligos$sequence) + nchar("GGTTTAGCCGGCGTG"), regexpr("GGCCGTAATGGCC", readOligos$sequence)-1)

oligoDesign = readOligos$sequence
names(oligoDesign) = readOligos$oligoID

# recode the allele strings!!!!
# (allele strings had been coded as 0/1 in the original design, and were mangled by excel at some point to cut of leading zeros – fix here, and replace with 1/2 encoding
# ("1" = BY, "2" = RM)
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
#save(readOligos, file="R_readOligos")
# this is one of the files we'll need to work with the actual samples later on


# process all the barcodes, no culling
# pull out the respective first oligo sequence
temp = sapply(bcCounts[,5], function(x){strsplit(x, ",")[[1]][1]})
names(temp) = NULL
bcCountsTopOligoOnly = bcCounts
bcCountsTopOligoOnly[,5] = temp
colnames(bcCountsTopOligoOnly) = c("bcSeq", "count", "oligoCount", "oligoCountDistribution", "oligoSeq", "fractionTopOligo")
#save(bcCountsTopOligoOnly, file="R_bcCountsTopOligoOnly_160629")
# another output from the annotation runs that will be used later

# for GEO: write to file
#write.table(bcCountsTopOligoOnly, file="TSS_bcCountsTopOligoOnly_160629.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


# now find the oligo from the design which matches the observed sequences
bcCountsAssigned = merge(bcCountsTopOligoOnly, readOligos, by.x = "oligoSeq", by.y = "sequence")

nrow(bcCountsAssigned)
# 9,153,139

length(unique(bcCountsAssigned$oligoSeq))
# 6,564 oligos covered; 91%

#save(bcCountsAssigned, file = "R_bcCountsAssigned_160629")
# the final of three output files from the annotation run that will be used later



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

pdf("oligoCountsAssigned.pdf", width=12, height=8)
par(mfrow=c(2,3))

hist(assignedOligosBarcodeCountsWithZeros, breaks=50, xlab="number of barcodes per oligo", main="all oligos")
hist(assignedOligosBarcodeCountsWithZeros[assignedOligosBarcodeCountsWithZeros <= 5000], breaks=50, xlab="number of barcodes per oligo", main="oligos with barcode number < 5000")
hist(log10(assignedOligosBarcodeCountsWithZeros + 1), breaks=50, xlab="log10(number of barcodes per oligo)", main="all oligos")

hist(assignedOligosBarcodeCounts[,2], breaks=50, xlab="number of barcodes per oligo", main="oligos with at least 1 barcode")
hist(assignedOligosBarcodeCounts[,2][assignedOligosBarcodeCounts[,2] <= 5000], breaks=50, xlab="number of barcodes per oligo", main="oligos with barcode number 0 < x < 5000")
hist(log10(assignedOligosBarcodeCounts[,2] + 1), breaks=50, xlab="log10(number of barcodes per oligo)", main="oligos with at least 1 barcode")

dev.off()

summary(assignedOligosBarcodeCountsWithZeros)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0      28     436    1269    1662   18750

summary(assignedOligosBarcodeCounts[,2])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#1      67     590    1394    1849   18750


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
#1.000    2.000    4.000    7.202    9.000 4303




###################################
# examine oligo dropouts

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

# dropout due to G(G) at the first dinlucleotides
