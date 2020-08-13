# this was helper code to throw this step onto the cluster before we had our own private compute node on which to do unlimited-length computation

library(dplyr)

setwd("YOUR_PATH")
load("R_barcodeOligoCounts_160627")


# for single thread:
doFunc <- function(x){
    ret = data.frame(count=sum(x[,3]), oligoCount=nrow(x), oligoDistribution=paste(x[,3], collapse=","), oligoSeqs=paste(x[,2], collapse=","), fractionCountsPerTopOligo=max(x[,3])/sum(x[,3]), stringsAsFactors=FALSE)
    ret
}

bcCounts = barcodeOligoCounts %>% group_by(barcodes) %>% do(doFunc(.))
bcCounts <- bcCounts[order(bcCounts$count, decreasing=TRUE),]

save(bcCounts, file="R_bcCounts_singleThread_160627")
