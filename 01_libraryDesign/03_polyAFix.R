# read the design, reverse complement the oligo if it has more As than Ts

a = read.table("jointDesignFinal27k_150529.txt", stringsAsFactors=FALSE)
AContent = sapply(a[,"sequence"], function(x){b=strsplit(x,"")[[1]]; length(b[b=="A"])/nchar(x)})
TContent = sapply(a[,"sequence"], function(x){b=strsplit(x,"")[[1]]; length(b[b=="T"])/nchar(x)})

pdf("AContent.pdf")
hist(AContent, breaks=100)
dev.off()

# no obvious breaks in A content, perhaps cut at 40%
reverseComplementFlag = AContent > 0.4
# for 40%, this is 291 oligos

# or, can RC all with more As than Ts
reverseComplementFlag = AContent > TContent
# this is TRUE for 15059 oligos, 55.8%
# as one might expect - promoters are A rich, apparently more so than T rich


library("Biostrings")

aRC = a
aRC[,"sequence"][reverseComplementFlag] <- sapply(aRC[,"sequence"][reverseComplementFlag], function(x){toString(reverseComplement(DNAString(x)))})
aRC = cbind(aRC, reverseComplementFlag)

# test if the A>T criterium worked
didItWork = sapply(aRC[,"sequence"], function(x){b=strsplit(x,"")[[1]]; length(b[b=="A"])/nchar(x)}) > sapply(aRC[,"sequence"], function(x){b=strsplit(x,"")[[1]]; length(b[b=="T"])/nchar(x)})
summary(didItWork)
# these are all false => good

#pdf("AContentAfterRC_40PercentA.pdf")
#pdf("AContentAfterRC_ALargerT.pdf")
hist(sapply(aRC[,"sequence"], function(x){b=strsplit(x,"")[[1]]; length(b[b=="A"])/nchar(x)}), breaks=100, main="AContent after reverse complement", xlab="AContent")
dev.off()


#write.table(aRC, file="jointDesignFinal27k_150611_REVERSECOMPLEMENT.txt", quote=FALSE, sep="\t")
#write.table(aRC, file="jointDesignFinal27k_150612_ReverseComplementWhenALargerT.txt", quote=FALSE, sep="\t")
