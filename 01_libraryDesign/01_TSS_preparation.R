# wrangle the data from Pelechano et al to get the most used TSS for each gene

genes = read.table("SGD_Protein_coding_Genes_and_RNAGenes_ncOverlapsRemoved.gff3", stringsAsFactors=FALSE)
genes = genes[genes[,3] %in% c("gene", "ncRNA", "tRNA", "snoRNA", "snRNA"),]
rownames(genes) = sapply(as.character(genes[,9]), function(x){strsplit(strsplit(x, ";")[[1]][1], "=")[[1]][2]})

# remove 2micron & mitchondrion
genes = genes[genes[,1] != "2-micron",]
genes = genes[genes[,1] != "chrM",]



# read in the Pelechano 2013 transcripts
# their original file uses spaces to separate columns AND ONE OF THE COLUMNS IS FREE TEXT WITH SPACES
# needed to operate on in excel before we can use it
pelechano1 = read.table("S2_tcd_mTIFAnno_mod.txt", stringsAsFactors=FALSE, sep="\t", head=TRUE)
pelechanoGeneNamesAndTypes = t(sapply(pelechano1[,1], function(x){
    splitted = strsplit(x, " ")[[1]]
    splitted[c(length(splitted), 7)]
}))

pelechano2 = read.table("S2_tcd_mTIFAnno_mod_1st6columns.txt", stringsAsFactors=FALSE, sep="\t", head=TRUE)
pelechano2 = cbind(pelechano2, pelechanoGeneNamesAndTypes)

# restrict to those transcripts that span a given gene completely
# i.e. avoid truncated forms
pelechano2 = pelechano2[pelechano2[,7] %in% rownames(genes) & pelechano2[,8] == "Covering",]

# now, for each gene, find the most frequent 5prime end
pelechano5PrimeEnds = sapply(as.character(unique(pelechano2[,7])), function(x){
    print(x)
    these5Primes = pelechano2[pelechano2[,7] == x, c(3, 5)]
#    print(these5Primes)
# the 5prime ends recur (due to different 3prime ends)
# add them up for 5prime
    theseSummed5Primes = sapply(unique(these5Primes[,1]), function(y){sum(these5Primes[these5Primes[,1] == y,2])})
#    print(theseSummed5Primes)
    unique(these5Primes[,1])[theseSummed5Primes == max(theseSummed5Primes)][1]
})

# for use with the countVariantsPerBin function, slide these UTRs into gene table
genesWithPelechano = genes
genesWithPelechano[rownames(genesWithPelechano) %in% names(pelechano5PrimeEnds) & genesWithPelechano[,7]=="+", 4] <- pelechano5PrimeEnds[rownames(genesWithPelechano[rownames(genesWithPelechano) %in% names(pelechano5PrimeEnds) & genesWithPelechano[,7]=="+", ])]
genesWithPelechano[rownames(genesWithPelechano) %in% names(pelechano5PrimeEnds) & genesWithPelechano[,7]=="-", 5] <- pelechano5PrimeEnds[rownames(genesWithPelechano[rownames(genesWithPelechano) %in% names(pelechano5PrimeEnds) & genesWithPelechano[,7]=="-", ])]

# save the genesWithPelechano
save(genesWithPelechano, file="R_genesWithPelechano")




##########################
# unique TSSs
# the object will be a data frame that holds the gene, the strand, the coordinate, and the type (AUG or UTRs)
# slow and ugly code
# note that the TSS_start is the first base of the ORF - ie this base should *not* be part of the oligo, but the bases immediately upstream of it!

TSSs = t(sapply(rownames(genes), function(x){
    if (genes[x,7]== "+"){return(c(x, as.character(genes[x, c(1, 4, 7)]), "AUG"))}
    if (genes[x,7]== "-"){return(c(x, as.character(genes[x, c(1, 5, 7)]), "AUG"))}
}))
TSSs = rbind(TSSs, t(sapply(rownames(genesWithPelechano), function(x){
    if (genesWithPelechano[x,7]== "+"){return(c(x, as.character(genesWithPelechano[x, c(1, 4, 7)]), "5UTR"))}
    if (genesWithPelechano[x,7]== "-"){return(c(x, as.character(genesWithPelechano[x, c(1, 5, 7)]), "5UTR"))}
})))
rownames(TSSs) = paste(TSSs[,1], TSSs[,5], sep="_")
# throw out redundant UTR-based definitions (i.e. when there is no UTR)
countOccurence = table(paste(TSSs[,1], TSSs[,2], TSSs[,3], TSSs[,4], sep="_"))
TSSs = TSSs[TSSs[,5] == "AUG" | (TSSs[,5] == "5UTR" & countOccurence[paste(TSSs[,1], TSSs[,2], TSSs[,3], TSSs[,4], sep="_")] == 1),]

TSSs = data.frame(TSSs[,1:2], as.integer(TSSs[,3]), TSSs[,4:5], stringsAsFactors=FALSE)
colnames(TSSs) = c("gene", "chr", "TSS", "strand", "TSSType")

dim(TSSs[TSSs[,5]=="AUG",])
# 6947 for AUG, 4881 for 5UTR
# total is 11,828 TSSs

# decorate TSSs with a column indictating whether each AUG also has a 5UTR annotation
also5UTR = apply(TSSs, 1, function(x){
print(x)
    if (x["TSSType"] == "5UTR"){return(NA)}
    else{
        x["gene"] %in% TSSs[TSSs[,"TSSType"] == "5UTR",][,1]
    }
})
TSSs = cbind(TSSs, also5UTR)

save(TSSs, file="R_TSSs_140925")
