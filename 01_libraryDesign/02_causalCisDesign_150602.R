
# note that read.table by default suppresses lines starting with '#'
#allVariants = read.table("/scratch/kruglyak/raid1/home/falbert/BY_RM_sequences/data/BAM_sacCer3/RM.var.flt.MoreIndels.vcf", stringsAsFactors=FALSE)

# as of 05/22/2015, switch over to Josh's GATK based version:
# can directly read in the version where the chromomse ends are already blocked out
allVariants = read.table("A_forFrank051815.sort.blocked.vcf", stringsAsFactors=FALSE)
rownames(allVariants) = paste(allVariants[,1], allVariants[,2], allVariants[,4], allVariants[,5], sep="_")
# there are 2,520 INDELS in the BY RM data
colnames(allVariants) = c("chr", "pos", "JoshName", "ref", "alt", "qual", "dot", "fields", "GT", "info1", "info2")


genes = read.table("SGD_Protein_coding_Genes_and_RNAGenes_ncOverlapsRemoved.gff3", stringsAsFactors=FALSE)
genes = genes[genes[,3] %in% c("gene", "ncRNA", "tRNA", "snoRNA", "snRNA"),]
rownames(genes) = sapply(as.character(genes[,9]), function(x){strsplit(strsplit(x, ";")[[1]][1], "=")[[1]][2]})
# remove 2micron & mitchondrion
genes = genes[genes[,1] != "2-micron",]
genes = genes[genes[,1] != "chrM",]
colnames(genes) = c("chr", "source", "type", "start", "end", "dot1", "strand", "dot2", "info")

# TSS we made from Pelechano et al
load("/scratch/kruglyak/raid1/home/falbert/causalCisVariants/design/R_TSSs_140925")


# note that this version is NOT smart enough to know about gene start/stop with respect to strand
# it also outputs the TSSs and the strand of the gene used to pull the variants
# offset is to shift the TSS by some distance (useful for tiling)
getVariantsPerTSS <- function(theseGenes, theseVariants, thisDistance, offset=0){
	ret1 = lapply(rownames(theseGenes), function(y){
					x = theseGenes[y,]
					ret = NA
					thisTSS = NA
								
					if (x[,"strand"] == "+"){
						thisTSS = x[,"TSS"] - offset
						ret = theseVariants[theseVariants[,"chr"] == x[,"chr"] & theseVariants[,"pos"] > (thisTSS - thisDistance) & theseVariants[,"pos"] < thisTSS,]
					}
					if (x[,"strand"] == "-"){
						thisTSS = x[,"TSS"] + offset		
						ret = theseVariants[theseVariants[,"chr"] == x[,"chr"] & theseVariants[,"pos"] > thisTSS & theseVariants[,"pos"] < (thisTSS + thisDistance),]
					}
				ret = cbind(ret, rep(thisTSS, nrow(ret)), rep(x[,"strand"], nrow(ret)), rep(thisDistance, nrow(ret)), rep(offset, nrow(ret)))
				colnames(ret)[(ncol(ret)-3):ncol(ret)] = c("TSS_pos", "TSS_strand", "oligoLength", "tilingOffset")
# catch indels that are too close to the end of the oligo:
				ret = ret[(ret[,"pos"] + nchar(ret[,"ref"]) < ret[,"TSS_pos"] & ret[,"TSS_strand"] == "+")
                        |
                        (ret[,"pos"] + nchar(ret[,"ref"]) < ret[,"TSS_pos"] + ret[,"oligoLength"] & ret[,"TSS_strand"] == "-")
                    ,]
										# for some reason, strand becomes a factor. change it back:
				ret[,"TSS_strand"] <- as.character(ret[,"TSS_strand"])
				ret
		})
	names(ret1) = rownames(theseGenes)
	ret1
}



# we need a function to turn single variants into something the oligo maker can use
# for each variant, center on it and adjust the TSS site etc accordingly
# the output of this function will be different from the one above:
# the one above makes lists of the same length for each run, irrespective of whether there are variants in that window or not
# this one will add one block per variant to a growing output
# the idea is to use it to pull all variants further away than "offset" from the gene start (AUG, presumably)
# and that are within either some "thisMaxDist" and/or in the entire distance until the next gene
# in here, ""offset" is the distance from the gene start within which no variants are pulled
# 05/22/2015
# not sure if this ever happens, but:
# a long deletion, that is farther away than "offset", but long enough to reach the gene start
# would be kept by this design (unless it is too long to fit in the oligo)
# added code in the variant selection to catch this - such deletions are now excluded

getSingleVariantsPerTSS <- function(theseGenes, theseVariants, thisOligoLength, offset=0, maxDist = Inf, limitToUpstreamRegion=TRUE, thisMaxOligoLength=149){
    # make sure the genes are sorted by coordinate
    theseGenes = theseGenes[order(theseGenes[,"chr"], theseGenes[,"start"], theseGenes[,"end"]),]
    #    print(theseGenes[,1:7])
    ret1 = list()
    for(y in 1:nrow(theseGenes)){
        print(rownames(theseGenes)[y])
        thisGene = theseGenes[y,]
        subGenes = theseGenes[theseGenes[,"chr"] == thisGene[,"chr"],]
        thisMaxDist = maxDist
        # this is only for chromsome ends:
        maxDistHardLimit = 1000
        retVars = NA
        if (thisGene[,"strand"] == "+"){
            if (limitToUpstreamRegion){
                if(y==1){ thisMaxDist = maxDistHardLimit }
                else {
                    if (theseGenes[y-1,"chr"]!=thisGene[,"chr"]) { thisMaxDist = maxDistHardLimit }
                    # in fact, simply looking at the neighboring element does not work becuase of overlapping genes
                    # need to find the closest gene end instead
                    # else {thisMaxDist = min(c(maxDist, thisGene[,"start"] - theseGenes[y-1,"end"]))}
                     else {   thisMaxDist = min(c(maxDist, thisGene[,"start"] - sort(subGenes[,"end"][subGenes[,"end"] < thisGene[,"start"]], decreasing=TRUE)[1]))
                    }
                }
            }
            #print(c(rownames(theseGenes)[y], thisMaxDist))
            # 05/22/2015
            # edit this to take into account long deletions that, when centering oligos on them, could make the oligo reach into the gene
            # this problem does not arise on the minus strand because the annotation is always on the plus
            #retVars = theseVariants[theseVariants[,"chr"] == thisGene[,"chr"] & theseVariants[,"pos"] > (thisGene[,"start"] - thisMaxDist) & theseVariants[,"pos"] < thisGene[,"start"] - offset,]
            retVars = theseVariants[theseVariants[,"chr"] == thisGene[,"chr"] &
                theseVariants[,"pos"] > (thisGene[,"start"] - thisMaxDist) &
                ((theseVariants[,"pos"] + nchar(theseVariants[,"ref"]) < thisGene[,"start"] - offset & nchar(theseVariants[,"ref"]) > nchar(theseVariants[,"alt"])) |
                    (theseVariants[,"pos"] < thisGene[,"start"] - offset & nchar(theseVariants[,"ref"]) <= nchar(theseVariants[,"alt"])))
                ,]
        }
        if (thisGene[,"strand"] == "-"){
            if (limitToUpstreamRegion){
                if(y==nrow(theseGenes)){ thisMaxDist = maxDistHardLimit }
                else {
                    if (theseGenes[y+1,"chr"]!=thisGene[,"chr"]) { thisMaxDist = maxDistHardLimit }
                    #else {thisMaxDist = min(c(maxDist, theseGenes[y+1,"start"] - thisGene[,"end"]))}
                    else{
                        thisMaxDist = min(c(maxDist, sort(subGenes[,"start"][subGenes[,"start"] > thisGene[,"end"]], decreasing=FALSE)[1] - thisGene[,"end"]))
                    }
                }
            }
            #print(c(rownames(theseGenes)[y], thisMaxDist))
            retVars = theseVariants[theseVariants[,"chr"] == thisGene[,"chr"] & theseVariants[,"pos"] > thisGene[,"end"] + offset & theseVariants[,"pos"] < (thisGene[,"end"] + thisMaxDist),]
        }
        # now, step through the list of variants and add them to the result list one at a time
        if (nrow(retVars) > 0){
            for (j in 1:nrow(retVars)){
                tOL = thisOligoLength
                if (nchar(retVars[j,"alt"]) > nchar(retVars[j,"ref"]) & thisOligoLength + nchar(retVars[j,"alt"]) - nchar(retVars[j,"ref"]) > thisMaxOligoLength){
                    tOL = thisMaxOligoLength - nchar(retVars[j,"alt"])
                }
                if (thisGene[,"strand"] == "+"){
                    # 05/22/2015 add if clauses to keep oligos symmetric when there is a long deletion
                    if (nchar(retVars[j,"alt"]) >= nchar(retVars[j,"ref"])){ # SNP or insertion
                        ret = cbind(retVars[j,], rep(retVars[j,2] + ceiling(tOL/2), 1), rep(thisGene[,"strand"], 1), rep(tOL, 1), rep(offset, 1))
                    }
                    else{ # deletion
                        ret = cbind(retVars[j,], rep(retVars[j,2] + ceiling(tOL/2) + ceiling(abs(nchar(retVars[j,"alt"]) - nchar(retVars[j,"ref"]))/2), 1), rep(thisGene[,"strand"], 1), rep(tOL, 1), rep(offset, 1))
                    }
                }
                if (thisGene[,"strand"] == "-"){
                    # 05/22/2015 add if clause
                    if (nchar(retVars[j,"alt"]) >= nchar(retVars[j,"ref"])){ # SNP or insertion
                        ret = cbind(retVars[j,], rep(retVars[j,2] - ceiling(tOL/2), 1), rep(thisGene[,"strand"], 1), rep(tOL, 1), rep(offset, 1))
                    }
                    else { # deletion
                        ret = cbind(retVars[j,], rep(retVars[j,2] - ceiling(tOL/2) + ceiling(abs(nchar(retVars[j,"alt"]) - nchar(retVars[j,"ref"]))/2), 1), rep(thisGene[,"strand"], 1), rep(tOL, 1), rep(offset, 1))
                    }
                }
                colnames(ret)[(ncol(ret)-3):ncol(ret)] = c("TSS_pos", "TSS_strand", "oligoLength", "tilingOffset")
                # catch indels that are too close to the end of the oligo (i.e., deletions that would be too big):
                ret = ret[(ret[,"pos"] + nchar(ret[,"ref"]) < ret[,"TSS_pos"] & ret[,"TSS_strand"] == "+")
                    |
                    (ret[,"pos"] + nchar(ret[,"ref"]) < ret[,"TSS_pos"] + ret[,"oligoLength"] & ret[,"TSS_strand"] == "-")
                ,]
                if (nrow(ret) > 0){
                    # for some reason, strand becomes a factor. change it back:
                    ret[,"TSS_strand"] <- as.character(ret[,"TSS_strand"])
                    rownames(ret) <- paste(retVars[j,1], retVars[j,2], retVars[j,4], retVars[j,5], sep="_")
                    ret1 = c(ret1, list(ret))
                    names(ret1)[length(ret1)] <- paste(rownames(theseGenes)[y], retVars[j,1], retVars[j,2], retVars[j,4], retVars[j,5], sep="_")
                }
            }
        }
    }
    ret1
}

# tests:
#test = getSingleVariantsPerTSS(genes[14:20,], allVariants, thisOligoLength=137, offset=100)
#makeOligos(test[[25]], indelPadder=0)
#makeOligos(test[[26]], indelPadder=0)
#makeOligos(getSingleVariantsPerTSS(genes[14:20,], allVariants, thisOligoLength=144, offset=100)[[26]], indelPadder=0)
#makeOligos(upstreamVariantsAUGCentered1kMaxDist[["tY(GUA)J2_chrX_543602_GAGTGACAAAGTATCAAAGAAGGAATAAAAAAAAGACATTATTATGTTGAACTAGTGACACCCTTTTGTGGATTCCTACATACCCAGAGAGAACTCCT_G"]], indelPadder=0)
#makeOligos(getSingleVariantsPerTSS(genes["tY(GUA)J2",], allVariants, thisOligoLength=oligoLength, offset=round(oligoLength/2), maxDist = 1000)[[6]], indelPadder=0)
#makeOligos(getSingleVariantsPerTSS(genes["tS(AGA)H",], allVariants, thisOligoLength=oligoLength, offset=round(oligoLength/2), maxDist = 1000)[[2]], indelPadder=0)


##########################
# 05/22/2015
# now we're ready to design

# want three types:
# 1. native upstream AUG (with and without minimal promoter) <= this library had been designed and ordered, but was never run
# 2. native upstream UTR (with and without minimal promoter)
# 3. tiling of whole 1kb upstream of AUG (with minimal promoter only, one oligo pair per variant)

# the re-written getVariantsPerTSS() function should be all we need
# go through the gene/TSS list, and cycle through the offset

# create the offset steps
# reduce oligoLength to allow Rocky to add the extra restriction site upstream of the promoter
# the cloning sites are now 2*15 (priming sites) + 8 + 13 = 51 bp
# if we allow maybe 5bp for the longest indels
# that leaves 144 bp
# also confirmed that from the earlier design, only 92/45,622 oligos exceed that length, or 0.2%

oligoLength = 144
overlapLength = 10
maxDist = 600

# like this, the last step can stick out over the max length, so that the interval is fully covered
#offsetSteps = seq(from = 0, to = maxDist, by = oligoLength-overlapLength)

# like this, the last step fits inside of the interval, and the last part of the stretch is not covered
# in the worst case, the missing distance can be as long as the oligo length
offsetSteps = seq(from = 0, to = maxDist - oligoLength, by = oligoLength-overlapLength)


# 05/22/2015
# note in 2020: this is the TSS design
# note that we're only doing the first tiling step
# the more upstream variants will be individually targeted by an oligo pair
upstreamVariantsPerTSSBYRMTiledMulti = lapply(offsetSteps[1], function(thisOffset){print(thisOffset); getVariantsPerTSS(TSSs, allVariants, oligoLength, offset=thisOffset)})

# get upstream variants for all the genes
# note in 2020: this is the Upstream design
#upstreamVariantsAUGCentered = getSingleVariantsPerTSS(genes, allVariants, thisOligoLength=oligoLength, offset=round(oligoLength/2))
upstreamVariantsAUGCentered1kMaxDist = getSingleVariantsPerTSS(genes, allVariants, thisOligoLength=oligoLength, offset=round(oligoLength/2), maxDist = 1000)

#save(upstreamVariantsPerTSSBYRMTiledMulti, file="R_upstreamVariantsPerTSSBYRMTiledMulti_150522")
#save(upstreamVariantsAUGCentered1kMaxDist, upstreamVariantsAUGCentered, file="R_upstreamVariantsAUGCentered1kMaxDist_150522")

# when using all the intergenic space, this is 22,206 variants (from allVariants, 44% of total)
# when limiting to max 1k of intergenic space, this is 17,566 variants

# an example with a multiallele variant is upstreamVariantsPerTSSBYRMTiledMulti[[1]][2108]
# count how many MULTI variants there are in a design:
countMultis <- function(thisDesign){
        ret = sapply(thisDesign, function(x){
            sapply(x[,"alt"], function(theseAlts){
                length(strsplit(theseAlts, ",")[[1]]) - 1
            })
        })
        return(ret)
}
test = countMultis(upstreamVariantsAUGCentered1kMaxDist)
length(test[test>0])
# 15
test1 = c(); for (i in 1:length(upstreamVariantsPerTSSBYRMTiledMulti)){test1=c(test1, unlist(countMultis(upstreamVariantsPerTSSBYRMTiledMulti[[i]])))}
length(test1[test1>0])
# 10


# how many sequences is this?
# how many with > 0?
# make a matrix for the variant counts (the different distances are in columns)
upstreamVariantsCountPerTSSBYRMTiled = sapply(upstreamVariantsPerTSSBYRMTiledMulti, function(x){sapply(x, nrow)})

# the grand total of oligo sequences will be:
2*length(upstreamVariantsCountPerTSSBYRMTiled[upstreamVariantsCountPerTSSBYRMTiled > 0]) + 
sum(upstreamVariantsCountPerTSSBYRMTiled[upstreamVariantsCountPerTSSBYRMTiled > 1])
# 144 bp first tile only: 16751


# how many variants does this assay, and how often is each assayed?
# (consider different TSS, strand, and tiling)
countVariantsAssayedByDesign = function(thisDesign, theseVariants){
    numberVariantAssayed = data.frame(rep(0, nrow(theseVariants)), rep(0, nrow(theseVariants)))
    rownames(numberVariantAssayed) = rownames(theseVariants)
    colnames(numberVariantAssayed) = c("+", "-")

    if (class(thisDesign[[1]]) != "list"){thisDesign = list(thisDesign)}
    for (j in 1:length(thisDesign)){
        for(x in names(thisDesign[[j]])){
            print(c(j, x))
            thisGene = strsplit(x, "_")[[1]][1]
            thisStrand = genes[thisGene, 7]
            rownames(thisDesign[[j]][[x]])
            numberVariantAssayed[rownames(thisDesign[[j]][[x]]), thisStrand] = numberVariantAssayed[rownames(thisDesign[[j]][[x]]), thisStrand] + 1
        }
    }

    countRange = 0:max(numberVariantAssayed)
    # this matrix will have the + strand counts along the horizontal and the - counts along the vertical
    variantsAssayedCount = sapply(countRange, function(x){ sapply(countRange, function(y){nrow(numberVariantAssayed[numberVariantAssayed[,1] == x & numberVariantAssayed[,2] == y,])})})
    colnames(variantsAssayedCount) = countRange
    rownames(variantsAssayedCount) = countRange
    variantsAssayedCount
}

vAC = countVariantsAssayedByDesign(upstreamVariantsPerTSSBYRMTiledMulti, allVariants)
vAC1 = countVariantsAssayedByDesign(upstreamVariantsAUGCentered1kMaxDist, allVariants)
vACCombined = countVariantsAssayedByDesign(c(upstreamVariantsPerTSSBYRMTiledMulti, list(upstreamVariantsAUGCentered1kMaxDist)), allVariants)

# how many are assayed at least once?
sum(vAC) - vAC[1,1]
# 6789 for TSS
sum(vAC1) - vAC1[1,1]
# 13720
sum(vACCombined) - vACCombined[1,1]
# 15519
# note that these will be downsampled later to arrive at the design in the paper

##############
# create sequences and write
# for each sequence, store:

# gene/TSS this sequence is "for"
# chr
# 5prime position in genome
# strand
# sequence
# number of variants
# IDs of variants (format = chr_pos_ref_alt)
# position of variant in this sequence (from 5prime)
# allelic state at variants in this oligo (ref (0) or alt (1))



# extract sequences
library("Biostrings")
BYGenome = readDNAStringSet("sacCer3.fa")


# 03/12/2015
# updated to fish out first alt allele when there are multiple alts
# also changed the variant ID separator from comma to semicolon (otherwise conflict with the multi-alleles)

# 05/22/2015
# with the GATK file from Josh, we have an issue with overlapping ranges
# eg YBL029W_AUG: two variants chrII_166025_CTT_CT and chrII_166027_T_C
# a one-bp deletion immediately followed by a SNP
# in reality, two T substituted by one C, just coded weirdly
# notice how the second variant is contained in the ref string of the first
# this causes replaceAt to break
# the one at-a-time flips should be OK
# how do I create the entire alt sequence?
# should be able to do this by stepping through the variants one at a time - build up the all-alt sequence sequentially

# this function pulls out sequence from the genome and adds variants that were collected above
makeOligos <- function(theseVariants, geneName="dummy", indelPadder = 0){
# indelPadder serves to pull more sequence within which we can inject indels
# then, cut back to the oligoLength for the output
		
	if (nrow(theseVariants) < 2){
# 0 if none, 2 if one variant
		ret = data.frame(matrix(NA, ncol=10, nrow=(2*nrow(theseVariants))), stringsAsFactors=FALSE)
	}
	if (nrow(theseVariants) >= 2){
# two for the complete flips, one more for each variant
		ret = data.frame(matrix(NA, ncol=10, nrow=(2 + nrow(theseVariants))), stringsAsFactors=FALSE)
	}
	colnames(ret) = c("gene", "chr", "TSS_pos", "strand", "sequence", "numberVariants", "variantIDs", "variantPositions", "alleleString", "tilingOffset")
	ret[,"TSS_pos"] <- as.integer(ret[,"TSS_pos"])
	ret[,"numberVariants"] <- as.integer(ret[,"numberVariants"])	
	if (nrow(theseVariants) == 0){return(ret)}
	
	ret[,"gene"] <- rep(geneName, nrow(ret))
	ret[,"chr"] <- rep(theseVariants[1,"chr"], nrow(ret))
	ret[,"TSS_pos"] <- rep(theseVariants[1,"TSS_pos"], nrow(ret))
	ret[,"strand"] <- rep(theseVariants[1,"TSS_strand"], nrow(ret))
	ret[,"numberVariants"] <- rep(nrow(theseVariants), nrow(ret))
	ret[,"variantIDs"] <- rep(paste(rownames(theseVariants), collapse=";"), nrow(ret))
	ret[,"tilingOffset"] <- rep(theseVariants[1,"tilingOffset"], nrow(ret))
#	ret[,"variantPositions"] <- rep(paste(theseVariants[,"TSS_pos"] - theseVariants[,"pos"], collapse=";"), nrow(ret))

# pick the first allele from a multiallelic variant
    theseVariants[,"alt"] <- sapply(theseVariants[,"alt"], function(x){strsplit(x, ",")[[1]][1]})

# make a matrix that shows the allelic states
# 05/22/2015 change so that all-alt is at the end of the matrix
# note: here we had made the poor design choice of printing the allele states as 0 & 1, which got eaten by excel later on when it trimmed leading zeros. this is why there are fixes to allele strings in later code bits
	alleleMatrix = matrix(0, nrow=nrow(ret), ncol=nrow(theseVariants))
	if (nrow(theseVariants) > 1){
		for (i in 1:nrow(theseVariants)){
			alleleMatrix[1+i, i] <- 1
		}
	}
    alleleMatrix[nrow(alleleMatrix),] = rep(1, nrow(theseVariants))
    #    print(alleleMatrix)
	ret[,"alleleString"] <- apply(alleleMatrix, 1, function(x){paste(x, collapse="")})
	
# for pulling the reference sequence, note the 1bp correction necessitated by the fact that TSS_start is the first base of the gene!
	if (theseVariants[1,"TSS_strand"] == "+"){
		refSequence = toString(subseq(BYGenome[theseVariants[1,"chr"]], start = theseVariants[1, "TSS_pos"] - theseVariants[1, "oligoLength"] - indelPadder, end = theseVariants[1, "TSS_pos"]-1))
		replacePositions = indelPadder + theseVariants[1, "oligoLength"] - (theseVariants[,"TSS_pos"] - theseVariants[,"pos"]) + 1
# fill in the allelePositions in the read
		ret[,"variantPositions"] <- rep(paste(replacePositions - indelPadder, collapse=";"), nrow(ret))
	}
	
	if (theseVariants[1,"TSS_strand"] == "-"){
		refSequence = toString(subseq(BYGenome[theseVariants[1,"chr"]], start = theseVariants[1, "TSS_pos"] + 1, width = theseVariants[1, "oligoLength"] + indelPadder))
		replacePositions = theseVariants[,"pos"] - theseVariants[,"TSS_pos"]
#		"revcomp" the positions within the reads:
		ret[,"variantPositions"] <- rep(paste(theseVariants[1, "oligoLength"] - replacePositions + 1, collapse=";"), nrow(ret))
	}

# fill in ret with all refSequence
	ret[,"sequence"] <- rep(refSequence, nrow(ret))

# loop through the alleles and replace one allele at a time - directly change them in the ret data.frame
# note that replaceLetterAt() only does substitutions, not indels; use the more generic replaceAt() instead
# 05/22/2015:
# as we step through the individual swaps, also edit the all-alt sequence, one at a time
# need to keep track of position shifts in the oligo when we insert an indel!
# if there is only one variant, can do as before:
    if (nrow(alleleMatrix) == 2){
        ret[2,"sequence"] = as.character(replaceAt(DNAString(ret[2, "sequence"]), at = IRanges(start = replacePositions[1], width = nchar(theseVariants[,"ref"][1])), value = theseVariants[,"alt"][1]))
    }
# otherwise:
    else{
        allAltReplacePositions <- replacePositions
        for (i in 2:(nrow(alleleMatrix)-1)){
            ret[i,"sequence"] = as.character(replaceAt(DNAString(ret[i, "sequence"]), at = IRanges(start = replacePositions[which(alleleMatrix[i,] == 1)], width = sapply(theseVariants[,"ref"], nchar)[which(alleleMatrix[i,] == 1)]), value = theseVariants[,"alt"][which(alleleMatrix[i,] == 1)]))
            # 05/22/2015: this is the second call, now altering the last (all-alt) sequence
            # before we call, need to adjust the positions of the remaining variants in the all-alt
            # the adjustment is always caluclated, but does nothing (adds zero) unless the current variant is an indel
            # UPDATE 1/22/2020: it looks like this broke the all-RM string in cases where a big RM insertion shortened the base oligo. we then have a second allele string of "0" but no "11111". not sure how, but be aware.
            # single variant flips look fine though
            if (i < (nrow(alleleMatrix)-1)){
                allAltReplacePositions[(which(alleleMatrix[i,] == 1)+1):length(allAltReplacePositions)] <- allAltReplacePositions[(which(alleleMatrix[i,] == 1)+1):length(allAltReplacePositions)] + (nchar(theseVariants[,"alt"][which(alleleMatrix[i,] == 1)]) - nchar(theseVariants[,"ref"][which(alleleMatrix[i,] == 1)]))
            }
            ret[nrow(alleleMatrix),"sequence"] = as.character(replaceAt(DNAString(ret[nrow(alleleMatrix), "sequence"]), at = IRanges(start = allAltReplacePositions[which(alleleMatrix[i,] == 1)], width = sapply(theseVariants[,"ref"], nchar)[which(alleleMatrix[i,] == 1)]), value = theseVariants[,"alt"][which(alleleMatrix[i,] == 1)]))
        }
    }

# revcomp the sequences on the minus strand:
	if (theseVariants[1,"TSS_strand"] == "-"){
		ret[,"sequence"] <- sapply(ret[,"sequence"], function(x){toString(reverseComplement(DNAString(x)))})
		
	}
	
# cut back to the oligoLength (ie remove the indelPadder)
	if (indelPadder > 0){
		ret[,"sequence"] <- sapply(ret[,"sequence"], function(x){substr(x, start = nchar(x) - theseVariants[1, "oligoLength"] + 1, stop=nchar(x))})
	}
	# if indelPadder == 0, nothing was added
	# in that case, simply report the sequence, irrespective of whether it grew/shrank
	
	ret
}

# tests
#makeOligos(upstreamVariantsPerTSSBYRMTiledMulti[[1]][[2108]])
#makeOligos(upstreamVariantsPerTSSBYRMTiledMulti[[1]][["YBL029W_AUG"]], indelPadder=0)
#makeOligos(upstreamVariantsPerTSSBYRMTiledMulti[[1]][["YAL024C_AUG"]], indelPadder=0)


# let's join the single variants with the tiles

jointDesign = c(upstreamVariantsPerTSSBYRMTiledMulti, list(upstreamVariantsAUGCentered1kMaxDist))


# now loop through the variant lists and make the actual oligos
jointDesignOligos = lapply(jointDesign, function(x){
	ret = lapply(names(x), function(y){
							 print(y)
		makeOligos(x[[y]], y, indelPadder = 0)
	})
	names(ret) = names(x)
	ret
})


# some oligos are too long
# when this happens, shrink all (!) oligos in a block by the length it takes to shrink the largest to be exactly within range
# this might sometimes push out variants - need to adjust accordingly
# (in practice, variants are almost never pushed out)

# 149 is max length (51 bp go for cloning, see above)
maxOligoLength = 149
removeCounterShorten = 0
removeCounterShortenRemove = 0

# 05/22/2015:
# note that this leads to asymmetric oligos since the chew-back happens only from one side
# can't change this to bidirectional chewing here
# because the firstTile oligos start at the TSS - can't chew back from here
# the chewing always happens from the gene-distal end
# just leave as is for now. An example for asymmetry is tS(AGA)H_chrVIII_133342_T_TTTATGTTCTTTCCTTTCATATGTTTCG_72
# UPDATE 1/22/2020: I wonder if this is where the long RM insertions broke the all-RM string allele string assignment?
jointDesignOligosLengthFiltered = lapply(jointDesignOligos, function(x){
	lapply(x, function(y){
				 if(nrow(y) == 0){
						return(y)
				 }
				 else{
						ret = y
						theseSeqs = y[,"sequence"]
						if (sum(nchar(theseSeqs) > maxOligoLength) > 0){
							nBaseToCut = max(nchar(theseSeqs)) - maxOligoLength
							cutSeqs = substr(theseSeqs, nBaseToCut + 1, 999999999)
							theseVarPoss = as.numeric(strsplit(y[1, "variantPositions"], ";")[[1]])
							theseVarPoss = theseVarPoss - nBaseToCut
							varsToDelete = which(theseVarPoss < 0)
							
							ret[,"sequence"] = cutSeqs
							ret[,"variantPositions"] <- rep(paste(theseVarPoss, collapse=";"), nrow(ret))
							
							if (length(varsToDelete) > 0){
							print(y)
								ret[,"variantPositions"] <- rep(paste(theseVarPoss[-varsToDelete], collapse=";"), nrow(ret))
								ret[,"variantIDs"] <- rep(paste(as.character(strsplit(y[1, "variantIDs"], ";")[[1]])[-varsToDelete], collapse=";"), nrow(ret))
								ret[,"numberVariants"] <- ret[,"numberVariants"] - length(varsToDelete)
								ret[,"alleleString"] <- sapply(y[, "alleleString"], function(z){paste(as.character(strsplit(z, "")[[1]])[-varsToDelete], collapse="")})
							
								if (nrow(y) >= 2){
											ret = ret[-(varsToDelete + 2),]
								}
								if (nrow(y) == 2){
											ret = ret[-2,]
								}
								removeCounterShortenRemove <<- removeCounterShortenRemove + 1
							}
							
							if (nrow(ret) == 1){
								# hack to remove all lines from ret
								ret = ret[rep(FALSE, nrow(ret)) == TRUE,]
							}

							removeCounterShorten <<- removeCounterShorten + 1
						}
						ret
				 }
	})
})
# 05/22/2015
# we have a few longer variants that need to get axed:
# 8 removed
# number that needed shortening - probably because of the change in the single variant code that pre-shortens the ref sequence
# 67 shortened


# remove oligo sets where one of the oligos matches one of the restriction sites used for cloning:
# note that the recognition sequenes are already palindromic, so there is no need to revcomp them
avoidPatterns = c("GGCGCGCC", "CCTGCAGG", "GGCC[ACGT]{5}GGCC", "^CTGCAGG", "CCTGCA$", "CCTGCAG$", "GGCCATTACGGCC", "GGCCGTAATGGCC", "^TGCAGG", "N")
removeCounterByPattern = rep(0, length(avoidPatterns))
names(removeCounterByPattern) = avoidPatterns
removeCounter = 0

jointDesignOligosLengthRestrictionFiltered = lapply(jointDesignOligosLengthFiltered, function(x){
	lapply(x, function(y){
				 if(nrow(y) == 0){
						return(y)
				 }
				 else{
						ret = y
						theseSeqs = y[,"sequence"]
						decider = FALSE
						for(pat in avoidPatterns){
								if(length(grep(pat, theseSeqs)) > 0){
										print(pat)
										print(y)
# NOTE: <<- operator to modify variable in outer scope!!!
										removeCounterByPattern[pat] <<- removeCounterByPattern[pat] + 1
										decider = TRUE
								}
						}
						if (decider){
								# hack to remove all lines from ret
								ret = ret[rep(FALSE, nrow(ret)) == TRUE,]
								removeCounter <<- removeCounter + 1
						}
						ret
				 }
	})
})

# but fewer are removed here (including no more seqs with 'N'??):
t(t(removeCounterByPattern))
#GGCGCGCC            26
#CCTGCAGG            23
#GGCC[ACGT]{5}GGCC    8
#^CTGCAGG             1
#CCTGCA$              4
#CCTGCAG$             1
#GGCCATTACGGCC        0
#GGCCGTAATGGCC        0
#^TGCAGG              5
#N                    0
removeCounter
# 68

# Only looses a handful of sets


# attach the cloning sites
# skip adding priming (!) sites for now
# just attach the cloning sites and write out the genomic sequences and add priming sites later
# other wise: headache with primers
leftPrimerToBeAdded = ""
rightPrimerToBeAdded = ""
leftCloneSeqToBeAdded = "GGCCATTACGGCC"
rightCloneSeqToBeAdded = "GGCGCGCC"

jointDesignOligosLengthRestrictionFilteredSeqsAdded = lapply(jointDesignOligosLengthRestrictionFiltered, function(x){
	lapply(x, function(y){
			if(nrow(y) == 0){
				 return(y)
			}
			else{
				 ret = y
				 ret[,"sequence"] = paste(leftPrimerToBeAdded, leftCloneSeqToBeAdded, ret[,"sequence"], rightCloneSeqToBeAdded, rightPrimerToBeAdded, sep="")
			}
			ret
			})
})

#save(jointDesign, jointDesignOligos, jointDesignOligosLengthRestrictionFiltered, jointDesignOligosLengthRestrictionFilteredSeqsAdded, offsetSteps, file="R_jointDesign_150529")


# the design is a little bigger than we might like
# 05/22/2015: 51k and some change
# cut it down as follows:
# use all oligos immediately upstream of the AUG
# for the genes that have a distinct UTR, use the first 5UTR tile
# do not use any tiles further way from the TSS/AUG
# instead, use all the upstream single oligo guys
# to do this now, make a big dataframe with all the oligos

jointDesignAsTable = data.frame(matrix(NA, ncol=ncol(jointDesignOligosLengthRestrictionFilteredSeqsAdded[[1]][[1]]),
	nrow=sum(unlist(lapply(jointDesignOligosLengthRestrictionFilteredSeqsAdded, function(x){sapply(x, nrow)})))))
colnames(jointDesignAsTable) = colnames(jointDesignOligosLengthRestrictionFilteredSeqsAdded[[1]][[1]])
tempCount = 1
	
for (i in 1:length(jointDesignOligosLengthRestrictionFilteredSeqsAdded)){
print(i)
	for (j in 1: length(jointDesignOligosLengthRestrictionFilteredSeqsAdded[[i]])){
	print(j)
		thisThing = jointDesignOligosLengthRestrictionFilteredSeqsAdded[[i]][[j]]
		if (nrow(thisThing) > 0){
			jointDesignAsTable[tempCount:(tempCount + nrow(thisThing) - 1),] = thisThing
			tempCount = tempCount + nrow(thisThing)
		}
	}
}

# and filter from here, using the TSSs object for information on whether each gene has an annotated 5UTR

# stick on the flag for whether there is a UTR for a given AUG
temp = table(TSSs[,"gene"])
AUGUTRCount = sapply(jointDesignAsTable[,1], function(x){
	temp[strsplit(x, split="_")[[1]][1]]
})
# in this, 1 => only AUG; 2 => AUG and UTR

# this is now "AUG", "5UTR" or a chromosome name (for the single variant design)
AUGOrUTR = sapply(jointDesignAsTable[,1], function(x){
	strsplit(x, split="_")[[1]][2]
})
AUGOrUTR[!(AUGOrUTR %in% c("AUG", "5UTR"))] <- "singleVariant"

jointDesignAsTableWithFlags = cbind(paste(jointDesignAsTable[,"gene"], jointDesignAsTable[,"tilingOffset"], sep="_"), jointDesignAsTable, AUGUTRCount, AUGOrUTR)
colnames(jointDesignAsTableWithFlags)[1] <- "oligoBlock"

# are some oligos are longer than 200bp?
table(nchar(jointDesignAsTableWithFlags[,"sequence"]))
# all <= 170, i.e. total length <= 200 once priming sites are added

firstTileAUG = jointDesignAsTableWithFlags[,"AUGOrUTR"] == "AUG" & jointDesignAsTableWithFlags[,"tilingOffset"] == 0
# 9699 oligos
firstTileUTR = (jointDesignAsTableWithFlags[,"AUGOrUTR"] == "5UTR" & jointDesignAsTableWithFlags[,"tilingOffset"] == 0)
# 7003 oligos (0 shared with the first AUG)

upstreamTiles = jointDesignAsTableWithFlags[,"tilingOffset"] > 0 & (jointDesignAsTableWithFlags[,"AUGOrUTR"] == "AUG" | jointDesignAsTableWithFlags[,"AUGOrUTR"] == "5UTR")
# none! (didn't design any)

uniqueUpstreamTiles = (jointDesignAsTableWithFlags[,"AUGOrUTR"] == "5UTR" & jointDesignAsTableWithFlags[,"tilingOffset"] > 0) | (jointDesignAsTableWithFlags[,"AUGOrUTR"] == "AUG" & jointDesignAsTableWithFlags[,"tilingOffset"] > 0 & jointDesignAsTableWithFlags[,"AUGUTRCount"] == 1)
# none!

singleVariantTiles = jointDesignAsTableWithFlags[,"AUGOrUTR"] == "singleVariant"
# 35018 oligos

# single + first AUG + first 5UTR = 51,720 oligos


jointDesignAsTableWithFlags = cbind(jointDesignAsTableWithFlags, firstTileAUG, firstTileUTR, upstreamTiles, uniqueUpstreamTiles, singleVariantTiles)
jointDesignAsTableReduced = jointDesignAsTableWithFlags[firstTileAUG | firstTileUTR | singleVariantTiles,]
# collapse the last three columns
designType = rep(NA, nrow(jointDesignAsTableReduced))
designType[jointDesignAsTableReduced[,c("firstTileAUG")]] <- "firstTileAUG"
designType[jointDesignAsTableReduced[,c("firstTileUTR")]] <- "firstTileUTR"
designType[jointDesignAsTableReduced[,c("singleVariantTiles")]] <- "singleVariantTiles"
jointDesignAsTableReduced = cbind(jointDesignAsTableReduced[,1:11], designType)
# 51,720 oligos
# 05/22/2015 notice that in this design run, the reduction didn't so anything!

# now, how many variants does *this* assay?
numberVariantAssayed = data.frame(rep(0, nrow(allVariants)), rep(0, nrow(allVariants)))
rownames(numberVariantAssayed) = rownames(allVariants)
colnames(numberVariantAssayed) = c("+", "-")

for(x in unique(jointDesignAsTableReduced[,"oligoBlock"])){
	print(x)
	thisGene = strsplit(x, "_")[[1]][1]
	thisStrand = genes[thisGene, 7]
	theseVars = strsplit(jointDesignAsTableReduced[jointDesignAsTableReduced[,"oligoBlock"] == x, "variantIDs"][1], split=";")[[1]]
	numberVariantAssayed[theseVars, thisStrand] = numberVariantAssayed[theseVars, thisStrand] + 1
}

countRange = 0:max(numberVariantAssayed)
# this matrix will have the + strand counts along the horizontal and the - counts along the vertical
variantsAssayedCount = sapply(countRange, function(x){ sapply(countRange, function(y){nrow(numberVariantAssayed[numberVariantAssayed[,1] == x & numberVariantAssayed[,2] == y,])})})
colnames(variantsAssayedCount) = countRange
rownames(variantsAssayedCount) = countRange
variantsAssayedCount
#0    1    2   3 4
#0 30063 3947 1025 352 0
#1  3889 2924  661 206 0
#2   980  616  121  69 0
#3   362  221   62  43 0
#4     1    1    0   0 0

sum(variantsAssayedCount) - variantsAssayedCount[1,1]
# 15480


# now, attach control oligos to each design
# the four we designed for Hunter Fraser's indel and SNP:
# these had been designed separately (not described in paper)
# BY-ref, RM-double flip, indel only fliP, SNP only flip
# and two pairs in "UpStream" configuration, i.e. centered on the variant
hunterSeqs = c( "CTTCTTCCGGACTTGCCCATGATTAACCTAATCTTATACGAACTGAATTAAACTTTACGGTATTACCGATAGGAAACTTCTATTTTATGATTTTTTCGTTCGGGGACGGAACGAACAGGAAACAAAAAAAAAGGTACGATCCATTGT",
                "CTTCTTCCGGACTTGTCCATGATTAACCTAATCTTATACGAACTGAATTAAACTTTACGGTATTACCGATAGGAAACTTCTATTTTATGATTTTTTCGTTCGGGGACGGAACGAACAGGAAACAAAAAAAAAAAGGTACGATCCATTGT",
                "CTTCTTCCGGACTTGCCCATGATTAACCTAATCTTATACGAACTGAATTAAACTTTACGGTATTACCGATAGGAAACTTCTATTTTATGATTTTTTCGTTCGGGGACGGAACGAACAGGAAACAAAAAAAAAAAGGTACGATCCATTGT",
                "CTTCTTCCGGACTTGTCCATGATTAACCTAATCTTATACGAACTGAATTAAACTTTACGGTATTACCGATAGGAAACTTCTATTTTATGATTTTTTCGTTCGGGGACGGAACGAACAGGAAACAAAAAAAAAGGTACGATCCATTGT",
                "TACGGTATTACCGATAGGAAACTTCTATTTTATGATTTTTTCGTTCGGGGACGGAACGAACAGGAAACAAAAAAAAAGGTACGATCCATTGTATTCTCTACCCCCGTATATAAAACTAAGCTGAACAAGCCTGTTGTTTTGCTTTAC",
                "TACGGTATTACCGATAGGAAACTTCTATTTTATGATTTTTTCGTTCGGGGACGGAACGAACAGGAAACAAAAAAAAAAAGGTACGATCCATTGTATTCTCTACCCCCGTATATAAAACTAAGCTGAACAAGCCTGTTGTTTTGCTTTAC",
                "CTTCAGTTGAAAATTACAGTGAACACAACATCTTCCCCAACAGACCTACATTAAAACGCTTCTTCCGGACTTGCCCATGATTAACCTAATCTTATACGAACTGAATTAAACTTTACGGTATTACCGATAGGAAACTTCTATTTTATG",
                "CTTCAGTTGAAAATTACAGTGAACACAACATCTTCCCCAACAGACCTACATTAAAACGCTTCTTCCGGACTTGTCCATGATTAACCTAATCTTATACGAACTGAATTAAACTTTACGGTATTACCGATAGGAAACTTCTATTTTATG"
)
hunterOligoBlock = data.frame(
    c(rep("HunterControl_both", 4), rep("HunterControl_Indel", 2), rep("HunterControl_SNP", 2)),
    rep("YER044C", 8),
    rep("chrV", 8),
    rep(238182, 8),
    rep("-", 8),
    hunterSeqs,
    c(2,2,2,2,1,1,1,1),
    c(rep("hunterIndel_CTTTTTTTTT_CTTTTTTTTTTT,CTTTTTTTTTTTT;hunterSNP_G_A", 4), rep("hunterIndel_CTTTTTTTTT_CTTTTTTTTTTT,CTTTTTTTTTTTT", 2), rep("hunterSNP_G_A", 2)),
    c(rep("133;16", 4), "78", "78", "74", "74"),
    c("00", "11", "10", "01", "0", "1", "0","1"),
    rep(0, 8),
    rep("control", 8),
    stringsAsFactors = FALSE
)
colnames(hunterOligoBlock) = colnames(jointDesignAsTableReduced)
# attach adapters
hunterOligoBlock[,"sequence"] <- paste("GGCCATTACGGCC", hunterOligoBlock[,"sequence"], "GGCGCGCC", sep="")

# and one set of 200 oligos from Eran's paper
eranData = read.table("Sharon_2012_nbt_S2_mod.txt", sep="\t", stringsAsFactors=FALSE, head=TRUE)
eranData = eranData[complete.cases(eranData),]
eranData[,1] = paste("sharon2012", eranData[,1], sep="_")

# need to take off their cloning sites:
# note that this contains the 10bp barcode
trimFromLeft = "GGGGACCAGGTGCCGTAAG.........."
trimFromRight = "GCGATCCTAGGGCGATCA"
eranDataTrimmed = eranData
eranDataTrimmed[,"OligoSequence"] <- substr(eranDataTrimmed[,"OligoSequence"], start=nchar(trimFromLeft)+1, stop = nchar(trimFromLeft)+103)


avoidPatterns = c("GGCGCGCC", "CCTGCAGG", "GGCC[ACGT]{5}GGCC", "^CTGCAGG", "CCTGCA$", "^TGCAGG", "CCTGCAG$")

avoidOligo = sapply(avoidPatterns, function(pat){
    grepl(pat, eranDataTrimmed[,"OligoSequence"])
})
avoidOligo = rowSums(avoidOligo) > 0
# this costs at 10 of their oligos
# 3658 are left when killing T stretches

eranDataTrimmedFiltered = eranDataTrimmed[!avoidOligo,]
rownames(eranDataTrimmedFiltered) = paste("ID", eranDataTrimmedFiltered[,"LibraryID"], sep="_")

# let us just sample 200 oligos from their experiment:
nEranControl=200
positiveControls = sample(rownames(eranDataTrimmedFiltered), nEranControl)

# make sure they cover the space:
pdf("eranExpressionValues.pdf")
plot(log(eranDataTrimmedFiltered[,"ExpressionReplicate2"]), log(eranDataTrimmedFiltered[,"ExpressionReplicate1"]), pch=19, col="#00000022", xlab="expression replicate 1", ylab="expression replicate 2")
points(log(eranDataTrimmedFiltered[positiveControls,"ExpressionReplicate2"]), log(eranDataTrimmedFiltered[positiveControls,"ExpressionReplicate1"]), col="red")
dev.off()

eranOligoBlock = data.frame(
    rep("EranControl", nrow(eranDataTrimmedFiltered)),
    rep("control", nrow(eranDataTrimmedFiltered)),
    rep("control", nrow(eranDataTrimmedFiltered)),
    rep(NA, nrow(eranDataTrimmedFiltered)),
    rep("+", nrow(eranDataTrimmedFiltered)),
    eranDataTrimmedFiltered[,"OligoSequence"],
    rep(1, nrow(eranDataTrimmedFiltered)),
    rownames(eranDataTrimmedFiltered),
    rep("1", nrow(eranDataTrimmedFiltered)),
    rep("0", nrow(eranDataTrimmedFiltered)),
    rep(0, nrow(eranDataTrimmedFiltered)),
    rep("EranNBT12", nrow(eranDataTrimmedFiltered)),
    stringsAsFactors = FALSE
)
colnames(eranOligoBlock) = colnames(jointDesignAsTableReduced)
rownames(eranOligoBlock) = rownames(eranDataTrimmedFiltered)

# need to attach the padding sequence to Eran's oligos, as well as adapters:
eranOligoBlock[,"sequence"] <- paste("GGCCATTACGGCC", eranOligoBlock[,"sequence"], "GGCGCGCC", "TATAGAACGGAATCACCTCTGACAAGTAGCGTCAAATCGGT", sep="")
# now they're all 165 bp
# i.e., once the priming sites are added they will be the correct length, just like my design

eranOligoControls = eranOligoBlock[positiveControls,]
eranOligoControls[,"designType"] <-"control"

# assemble the final design table
# 05/22/2015 leave out Eran's whole set now!
# this will free up oligos for the single variants
jointDesignFinal = rbind(jointDesignAsTableReduced[jointDesignAsTableReduced[,"designType"] == "firstTileAUG",], hunterOligoBlock, eranOligoControls)
jointDesignFinal[,"designType"][jointDesignFinal[,"designType"] == "control"] <- "firstTileAUG"
jointDesignFinal = rbind(jointDesignFinal, jointDesignAsTableReduced[jointDesignAsTableReduced[,"designType"] == "firstTileUTR",], hunterOligoBlock, eranOligoControls)
jointDesignFinal[,"designType"][jointDesignFinal[,"designType"] == "control"] <- "firstTileUTR"
jointDesignFinal = rbind(jointDesignFinal, jointDesignAsTableReduced[jointDesignAsTableReduced[,"designType"] == "singleVariantTiles",], hunterOligoBlock, eranOligoControls)
jointDesignFinal[,"designType"][jointDesignFinal[,"designType"] == "control"] <- "singleVariantTiles"
#jointDesignFinal = rbind(jointDesignFinal, eranOligoBlock, hunterOligoBlock)
#jointDesignFinal[,"designType"][jointDesignFinal[,"designType"] == "control"] <- "EranNBT12"

# expect this to have 9699+7003+35018+3*8+3*200 = 52,344 oligos
# that is the case


# add rownames
rownames(jointDesignFinal) = c(
    paste("firstTileAUG", formatC(1:nrow(jointDesignFinal[jointDesignFinal[,"designType"] == "firstTileAUG",]), width=5, format="d", flag="0"), sep="_"),
    paste("firstTileUTR", formatC(1:nrow(jointDesignFinal[jointDesignFinal[,"designType"] == "firstTileUTR",]), width=5, format="d", flag="0"), sep="_"),
    paste("singleVariantTiles", formatC(1:nrow(jointDesignFinal[jointDesignFinal[,"designType"] == "singleVariantTiles",]), width=5, format="d", flag="0"), sep="_")
#    paste("EranNBT12", formatC(1:nrow(jointDesignFinal[jointDesignFinal[,"designType"] == "EranNBT12",]), width=5, format="d", flag="0"), sep="_")
)


# now, attach the priming sites for the designs:
jointDesignFinal[,"sequence"][jointDesignFinal[,"designType"] == "firstTileAUG"] <- paste("CGCGTCGAGTAGGGT", jointDesignFinal[,"sequence"][jointDesignFinal[,"designType"] == "firstTileAUG"], "CCAGCTTCACACGGC", sep="")
jointDesignFinal[,"sequence"][jointDesignFinal[,"designType"] == "firstTileUTR"] <- paste("CGATCGCCCTTGGTG", jointDesignFinal[,"sequence"][jointDesignFinal[,"designType"] == "firstTileUTR"], "CACGCCGGCTAAACC", sep="")
jointDesignFinal[,"sequence"][jointDesignFinal[,"designType"] == "singleVariantTiles"] <- paste("GGGTCACGCGTAGGA", jointDesignFinal[,"sequence"][jointDesignFinal[,"designType"] == "singleVariantTiles"], "GTGTGGCTGCGGAAC", sep="")
#jointDesignFinal[,"sequence"][jointDesignFinal[,"designType"] == "EranNBT12"] <- paste("GGTCGAGCCGGAACT", jointDesignFinal[,"sequence"][jointDesignFinal[,"designType"] == "EranNBT12"], "TCTGGGTGCGCATCC", sep="")

table(nchar(jointDesignFinal[,"sequence"]))
# looks good!


#save(jointDesignFinal, jointDesignAsTableWithFlags, jointDesignAsTableReduced, file="R_jointDesignAsTable_150524")

# write:
#write.table(jointDesignFinal, file="jointDesignFinal_150522.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)


########################################
# bring the total oligo count to >= 27k
# do this by keeping AUG, TSS, Eran
# and fill up with single variants
# focus on genes that have the most (and some with the least) evidence for ASE

# we need to fill out 27000 - (9699+7003+3*8+3*200) = 9674 oligos

# load the previous ASE analyses (had been done separately; combined ASE from Torabi et al ("Noori") and Albert 2014
load("R_jointASEData.RData")

# for each gene, how many variants (oligos) are in the single variant design?
toCutFrom = as.character(jointDesignFinal[jointDesignFinal[,"designType"] == "singleVariantTiles","oligoBlock"])
toCutFromGene = sapply(toCutFrom, function(x){strsplit(x, "_")[[1]][1]})
oligoCountsPerGene = as.numeric(table(toCutFromGene))
names(oligoCountsPerGene) = names(table(toCutFromGene))

# restrict ASE data to genes that actually have variants (i.e. oligos)
jointASEData = jointASEData[rownames(jointASEData) %in% names(oligoCountsPerGene),]
bonfCut = 0.05/nrow(jointASEData)
nrow(jointASEData[jointASEData[,"NooriB"] < bonfCut & jointASEData[,"AlbertB"] < bonfCut,])
# 117
bigASEGenes = rownames(jointASEData[jointASEData[,"NooriB"] < 0.05 & jointASEData[,"AlbertB"] < 0.05,])
sum(oligoCountsPerGene[bigASEGenes])
# if do not require same direction (which is interesting as well...): 4560 oligos, 451 genes => 5114 oligos to go
# let's do 1k from "no ASE" genes and fill up with random oligos

notSigASE = jointASEData[jointASEData[,"NooriB"] > 0.2 & jointASEData[,"AlbertB"] > 0.2,]
# notice that it matters where we place the abs() here:
# if we take the abs of the mean, we get anticorrelated measures between Albert & Noori
# with their individual abs not necessarily being super small
# instead, need to abs each individually
notSigASE = notSigASE[order(abs(notSigASE[,"NooriFC"]) + abs(notSigASE[,"AlbertFC"])/2),]
runningSumNotSigASE = cumsum(oligoCountsPerGene[rownames(notSigASE)])
smallASEGenes = rownames(notSigASE[runningSumNotSigASE <= 1000,])
sum(oligoCountsPerGene[smallASEGenes])
# 990
# 4124 left => can select 2062 random variants


pdf("ASEforSelectedSingleTiles_150522.pdf")
plot(jointASEData[bigASEGenes, c("NooriFC", "AlbertFC")], col="#FF000044", pch=19)
points(jointASEData[smallASEGenes, c("NooriFC", "AlbertFC")], col="#0000FF44", pch=19)
abline(h=0, lty=2, col="grey", lwd=2)
abline(v=0, lty=2, col="grey", lwd=2)
abline(0, 1, lty=2, col="grey", lwd=2)
dev.off()

# looks good, now find the corresponding oligos
oligosToPickFrom = jointDesignFinal[jointDesignFinal[,"designType"] == "singleVariantTiles",]
singleVariantOligosFor27k = oligosToPickFrom[toCutFromGene %in% bigASEGenes | toCutFromGene %in% smallASEGenes,]
singleVariantControlOligos = jointDesignFinal[jointDesignFinal[,"designType"] == "singleVariantTiles" & jointDesignFinal[,"oligoBlock"] %in% c("HunterControl_both", "HunterControl_Indel", "HunterControl_SNP", "EranControl"),]

# a random smattering of oligo pairs
# note that in spite of sampling, this code preserves the order of the oligos - it just skips those it doesn't want
oligosToPickFromForRandom = oligosToPickFrom[!(rownames(oligosToPickFrom) %in% c(rownames(singleVariantControlOligos), rownames(singleVariantOligosFor27k))),]
randomSingleOligoPairs = oligosToPickFromForRandom[oligosToPickFromForRandom[,"oligoBlock"] %in% sample(unique(oligosToPickFromForRandom[,"oligoBlock"]), 2062),]

#jointDesignFinal27k = jointDesignFinal[jointDesignFinal[,"designType"] %in% c("firstTileAUG", "firstTileUTR", "EranNBT12"),]
jointDesignFinal27k = jointDesignFinal[jointDesignFinal[,"designType"] %in% c("firstTileAUG", "firstTileUTR"),]

#jointDesignFinal27k = rbind(jointDesignFinal27k, singleVariantOligosFor27k, hunterOligoBlock, eranOligoControls)
# CAREFUL! now the controls are not primed yet!
# need to prime them before changing their library type
#jointDesignFinal27k[,"sequence"][jointDesignFinal27k[,"designType"] == "control"] <- paste("GGGTCACGCGTAGGA", jointDesignFinal27k[,"sequence"][jointDesignFinal27k[,"designType"] == "control"], "GTGTGGCTGCGGAAC", sep="")
#rownames(jointDesignFinal27k[,"designType"][jointDesignFinal27k[,"designType"] == "control"]) <- rownames()
#jointDesignFinal27k[,"designType"][jointDesignFinal27k[,"designType"] == "control"] <- "singleVariantTiles"

# this is safer:
jointDesignFinal27k = rbind(jointDesignFinal27k, singleVariantOligosFor27k, randomSingleOligoPairs, singleVariantControlOligos)

#save(jointDesignFinal27k, file="R_jointDesignFinal27k_150529")

# write:
#write.table(jointDesignFinal27k, file="jointDesignFinal27k_150529.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)





# now, how many variants does *this* assay?
numberVariantAssayed = data.frame(rep(0, nrow(allVariants)), rep(0, nrow(allVariants)))
rownames(numberVariantAssayed) = rownames(allVariants)
colnames(numberVariantAssayed) = c("+", "-")

designToCount = jointDesignFinal27k[!jointDesignFinal27k[,"oligoBlock"] %in% c("EranControl", "HunterControl_both", "HunterControl_Indel", "HunterControl_SNP"),]

for(x in unique(designToCount[,"oligoBlock"])){
    print(x)
    thisGene = strsplit(x, "_")[[1]][1]
    thisStrand = genes[thisGene, 7]
    theseVars = strsplit(designToCount[designToCount[,"oligoBlock"] == x, "variantIDs"][1], split=";")[[1]]
    numberVariantAssayed[theseVars, thisStrand] = numberVariantAssayed[theseVars, thisStrand] + 1
}

countRange = 0:max(numberVariantAssayed)
# this matrix will have the + strand counts along the horizontal and the - counts along the vertical
variantsAssayedCount = sapply(countRange, function(x){ sapply(countRange, function(y){nrow(numberVariantAssayed[numberVariantAssayed[,1] == x & numberVariantAssayed[,2] == y,])})})
colnames(variantsAssayedCount) = countRange
rownames(variantsAssayedCount) = countRange
variantsAssayedCount
# note that the numbers of a rerun may not exactly be the same as below, due to sampling in creating a smaller library
#0    1   2   3
#0 35938 3178 901 138
#1  3004  656 223  42
#2   981  176  59  13
#3   187   30  12   5

sum(variantsAssayedCount) - variantsAssayedCount[1,1]
# 9,505


