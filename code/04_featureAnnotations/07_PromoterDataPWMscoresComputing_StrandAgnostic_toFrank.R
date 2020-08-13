library(GenomicRanges)
library(Biostrings)
library(magrittr)
library(stringr)
library(dplyr)
library(readr)
library(filesstrings)

dir = "/Users/renga011/Desktop/3UTR/"

#genome -----
genome = readDNAStringSet(paste0(dir,"data/sacCer3.fa"))

#VariantInfo --------
  #load promoterData to get variant info
load(paste(dir, "data/resultsAndAnnotations_200501.rdata", sep = ""))

  #extract variant info
variantsDF <- data.frame(t(sapply(resultsAndAnnotations$variantID, function(x){strsplit(x, "_")[[1]]})), stringsAsFactors=FALSE)
colnames(variantsDF) <- c("chr", "pos", "ref", "alt", "strand")
variantsDF$SNP <- nchar(variantsDF[,3]) == 1 & nchar(variantsDF[,4]) == 1
  # in cases with commas in the alt allele, we only made the first
variantsDF$alt[str_detect(variantsDF$alt, ",")] <- sapply(variantsDF$alt[str_detect(variantsDF$alt, ",")], function(x){str_split(x, ",")[[1]][1]})
variantsGRanges <- GRanges(seqnames = variantsDF$chr, ranges = IRanges(start=as.numeric(variantsDF$pos), end=as.numeric(variantsDF$pos) + nchar(variantsDF$ref) - 1), strand = variantsDF$strand, ref=variantsDF$ref, alt=variantsDF$alt)

#PWMs ----------
  #load the cutoffs matrix
optimal_matrix_cutoffs = read_delim(paste0(dir,"data/optimal matrix cutoffs.cutoffs"), 
                                     "\t", escape_double = FALSE, col_names = FALSE, 
                                     trim_ws = TRUE)
colnames(optimal_matrix_cutoffs) = c("TF", "originalDataset", "TF_2(?)", "Experiment(?)", "cutoff", "(?)")
optimal_matrix_cutoffs = select(optimal_matrix_cutoffs, c("TF", "originalDataset", "cutoff")) #keeping only the columns that are not (?)

  #some TFs have isoforms and thus more than one cutoff.
  #we will just average across the isoforms to get an average cutoff
TFMatrix = as.data.frame(unique(optimal_matrix_cutoffs$TF))
colnames(TFMatrix) = "TF"
TFMatrix$cutoff = 0

for (transcriptionFactor in unique(optimal_matrix_cutoffs$TF)) {
  subsetted_matrix = filter(optimal_matrix_cutoffs, TF == transcriptionFactor)
  avgCutoff = mean(subsetted_matrix$cutoff)
  TFMatrix[TFMatrix$TF == transcriptionFactor, "cutoff"] = avgCutoff
}

  #load the PWMs 
  #all PWMs were put into a single file using the command cat /Users/Desktop/3UTR/PromoterPWMs/PWM/* > allPWMs.txt
allPWMs =  read_delim(paste0(dir,"data/allPWMs.txt"), 
                      "|", escape_double = FALSE, col_names = FALSE, 
                      trim_ws = TRUE)
allPWMs$PositionSpecificScores = lapply(extract_numbers(allPWMs$X2, decimals = TRUE, negs = TRUE), as.numeric)

colnames(allPWMs) = c("nucleotide", "PWM_strings", "PositionSpecificScores")

TFMatrix$PWM_A = filter(allPWMs[,-2], nucleotide == "A")[,"PositionSpecificScores"]
TFMatrix$PWM_C = filter(allPWMs[,-2], nucleotide == "C")[,"PositionSpecificScores"]
TFMatrix$PWM_G = filter(allPWMs[,-2], nucleotide == "G")[,"PositionSpecificScores"]
TFMatrix$PWM_T = filter(allPWMs[,-2], nucleotide == "T")[,"PositionSpecificScores"]

for(i in 1:nrow(TFMatrix)){
  TFMatrix$PWMSize[i] = length(TFMatrix$PWM_A[[1]][[i]])
}

#Window analyses to compute sequence scores for all variants for the 196 TFs --------------

  variant_WithoutCutOff_plus = as.data.frame(matrix(data = 0, nrow = nrow(variantsDF), ncol = 2*nrow(TFMatrix)))
  colnames(variant_WithoutCutOff_plus) = paste(TFMatrix$TF, rep(c("differenceOfMax", "differenceOfMeans"), 
                                                           each = length(TFMatrix$TF)), sep = "_")
  variant_WithCutOff_plus =  as.data.frame(matrix(data = 0, nrow = nrow(variantsDF), ncol = 3*nrow(TFMatrix)))
  colnames(variant_WithCutOff_plus) = paste(TFMatrix$TF, 
                                          rep(c("differenceOfMax", "differenceOfMeans", "deltaSignificantBindingSites"), 
                                              each = length(TFMatrix$TF)), 
                                          sep = "_")
  
  variant_WithoutCutOff_plus = as.data.frame(matrix(data = 0, nrow = nrow(variantsDF), ncol = 2*nrow(TFMatrix)))
  colnames(variant_WithoutCutOff_plus) = paste(TFMatrix$TF, rep(c("differenceOfMax", "differenceOfMeans"), 
                                                                each = length(TFMatrix$TF)), sep = "_")
  variant_WithCutOff_plus =  as.data.frame(matrix(data = 0, nrow = nrow(variantsDF), ncol = 3*nrow(TFMatrix)))
  colnames(variant_WithCutOff_plus) = paste(TFMatrix$TF, 
                                            rep(c("differenceOfMax", "differenceOfMeans", "deltaSignificantBindingSites"), 
                                                each = length(TFMatrix$TF)), 
                                            sep = "_")
  
  variant_WithoutCutOff_minus = as.data.frame(matrix(data = 0, nrow = nrow(variantsDF), ncol = 2*nrow(TFMatrix)))
  colnames(variant_WithoutCutOff_minus) = paste(TFMatrix$TF, rep(c("differenceOfMax", "differenceOfMeans"), 
                                                                each = length(TFMatrix$TF)), sep = "_")
  variant_WithCutOff_minus =  as.data.frame(matrix(data = 0, nrow = nrow(variantsDF), ncol = 3*nrow(TFMatrix)))
  colnames(variant_WithCutOff_minus) = paste(TFMatrix$TF, 
                                            rep(c("differenceOfMax", "differenceOfMeans", "deltaSignificantBindingSites"), 
                                                each = length(TFMatrix$TF)), 
                                            sep = "_")
  
  variant_WithoutCutOff_agnostic = as.data.frame(matrix(data = 0, nrow = nrow(variantsDF), ncol = 2*nrow(TFMatrix)))
  colnames(variant_WithoutCutOff_agnostic) = paste(TFMatrix$TF, rep(c("differenceOfMax", "differenceOfMeans"), 
                                                                 each = length(TFMatrix$TF)), sep = "_")
  variant_WithCutOff_agnostic =  as.data.frame(matrix(data = 0, nrow = nrow(variantsDF), ncol = 3*nrow(TFMatrix)))
  colnames(variant_WithCutOff_agnostic) = paste(TFMatrix$TF, 
                                             rep(c("differenceOfMax", "differenceOfMeans", "deltaSignificantBindingSites"), 
                                                 each = length(TFMatrix$TF)), 
                                             sep = "_")
  
  
  TFScoresForVariants_plus = as.data.frame(TFMatrix$TF)
  colnames(TFScoresForVariants_plus) = "TF"
  
  TFScoresForVariants_minus = as.data.frame(TFMatrix$TF)
  colnames(TFScoresForVariants_minus) = "TF"
  
  TFScoresForVariants_agnostic = as.data.frame(TFMatrix$TF)
  colnames(TFScoresForVariants_agnostic) = "TF"
  
  
  for(thisVariant in 1:nrow(variantsDF)){
allele1 = variantsDF$ref[thisVariant]
allele2 = variantsDF$alt[thisVariant]

  TFScoresForThisVariant_plus = as.data.frame(TFMatrix$TF) # temporary matrix which contains the scores for all PWMs for this Variant
  rownames(TFScoresForThisVariant_plus) = TFMatrix$TF
  TFScoresForThisVariant_plus$SequenceScores_refSeq = 0
  TFScoresForThisVariant_plus$SequenceScores_altSeq = 0
  
  TFScoresForThisVariant_minus = as.data.frame(TFMatrix$TF) # temporary matrix which contains the scores for all PWMs for this Variant
  rownames(TFScoresForThisVariant_minus) = TFMatrix$TF
  TFScoresForThisVariant_minus$SequenceScores_refSeq = 0
  TFScoresForThisVariant_minus$SequenceScores_altSeq = 0
  
  TFScoresForThisVariant_agnostic = as.data.frame(TFMatrix$TF) # temporary matrix which contains the scores for all PWMs for this Variant
  rownames(TFScoresForThisVariant_agnostic) = TFMatrix$TF
  TFScoresForThisVariant_agnostic$SequenceScores_refSeq = 0
  TFScoresForThisVariant_agnostic$SequenceScores_altSeq = 0
  
  
  # for loop for PWMs for this variant begins here
  for(thisPWM in 1:nrow(TFMatrix)){

extractedRefSeqLeft <- as.character(subseq(genome[seqnames(variantsGRanges)[thisVariant]], 
                                           start=start(variantsGRanges)[thisVariant] - TFMatrix$PWMSize[thisPWM] + 1, 
                                           end=start(variantsGRanges)[thisVariant] - 1))
extractedRefSeqRight <- as.character(subseq(genome[seqnames(variantsGRanges)[thisVariant]], 
                                            start=start(variantsGRanges)[thisVariant] + nchar(allele1), 
                                            end=start(variantsGRanges)[thisVariant] + 
                                              nchar(allele1) + TFMatrix$PWMSize[thisPWM] - 2))

  #paste the sequences together with the alleles to product the ref and alt seq for this PWM
refSeq = paste(extractedRefSeqLeft, allele1, extractedRefSeqRight, sep="")
altSeq = paste(extractedRefSeqLeft, allele2, extractedRefSeqRight, sep="")

  # variant annotation is always + strand, so we can revComp the plus in case the strand is negative
if(as.character(strand(variantsGRanges[thisVariant])) == "-"){
  refSeq <- as.character(reverseComplement(DNAString(refSeq)))
  altSeq <- as.character(reverseComplement(DNAString(altSeq)))
}

#PWM Matrix for this TF
PWMMatrix = matrix(TFMatrix$PWM_A[[1]][[thisPWM]])
PWMMatrix = as.data.frame(cbind(PWMMatrix, 
                  TFMatrix$PWM_C[[1]][[thisPWM]], 
                  TFMatrix$PWM_G[[1]][[thisPWM]], 
                  TFMatrix$PWM_T[[1]][[thisPWM]]))
colnames(PWMMatrix) = c("A", "C", "G", "T")

#cutoff for this TF
cutoff = TFMatrix$cutoff[thisPWM]

#get sequence scores for all possible windows for this PWM

  for(seq in c(refSeq, altSeq)){
  maxStartIndexForThisPWM = nchar(seq)-TFMatrix$PWMSize[thisPWM] +1
  SequenceScoresForThisPWM_plus = c()
  SequenceScoresForThisPWM_minus = c()
  SequenceScoresForThisPWM = c()
  SequenceScoresForThisPWM_WithCutoff = c()

    #plus strand
  for(startIndex in 1:maxStartIndexForThisPWM){
    sequenceScore = 0
    startPosition = startIndex
    endPosition = startPosition + TFMatrix$PWMSize[thisPWM] - 1
    subSeq = substring(text = seq,first = startPosition, last = endPosition)
        for(position in 1:nchar(subSeq)){
          nucleotideAtPosition = substring(text = subSeq, first = position, last = position)
          sequenceScore = sequenceScore + PWMMatrix[position, nucleotideAtPosition]
        }
    SequenceScoresForThisPWM_plus[subSeq] = sequenceScore
  }
  
    #minus strand
  seq_minus = as.character(reverseComplement(DNAString(seq)))
  for(startIndex in 1:maxStartIndexForThisPWM){
    sequenceScore = 0
    startPosition = startIndex
    endPosition = startPosition + TFMatrix$PWMSize[thisPWM] - 1
    subSeq = substring(text = seq_minus,first = startPosition, last = endPosition)
    for(position in 1:nchar(subSeq)){
      nucleotideAtPosition = substring(text = subSeq, first = position, last = position)
      sequenceScore = sequenceScore + PWMMatrix[position, nucleotideAtPosition]
    }
    SequenceScoresForThisPWM_minus[subSeq] = sequenceScore
  }
  
  #if refSeq assign the list of sequence scores to the refSeq column 
  computeTFScoreForThisVariant = function(seq, TFScoresForThisVariant, 
                                          SequenceScoresForThisPWM, SequenceScoresForThisPWM_WithCutoff){
    if(seq == refSeq){
      TFScoresForThisVariant$SequenceScores_refSeq[thisPWM] = list(SequenceScoresForThisPWM)
      TFScoresForThisVariant$SequenceScores_refSeq_WithCutoff[thisPWM] = list(SequenceScoresForThisPWM_WithCutoff)
      TFScoresForThisVariant$numberOfSignificantBindingSites_refSeq[thisPWM] = length(SequenceScoresForThisPWM_WithCutoff[!(is.na(SequenceScoresForThisPWM_WithCutoff))])
    }
    #if altSeq assign the list of sequence scores to the altSeq column
    if(seq == altSeq){
      TFScoresForThisVariant$SequenceScores_altSeq[thisPWM] = list(SequenceScoresForThisPWM)
      TFScoresForThisVariant$SequenceScores_altSeq_WithCutoff[thisPWM] = list(SequenceScoresForThisPWM_WithCutoff)
      TFScoresForThisVariant$numberOfSignificantBindingSites_altSeq[thisPWM] = length(SequenceScoresForThisPWM_WithCutoff[!(is.na(SequenceScoresForThisPWM_WithCutoff))])
    }
    
    return(TFScoresForThisVariant)
  }
  
    #agnostic -----
  SequenceScoresForThisPWM_agnostic = c(SequenceScoresForThisPWM_plus, SequenceScoresForThisPWM_minus)
  SequenceScoresForThisPWM_WithCutoff_agnostic = ifelse(SequenceScoresForThisPWM_agnostic >= cutoff, SequenceScoresForThisPWM_agnostic, NA) #strong binding only 
  SequenceScoresForThisPWM_agnostic = ifelse(SequenceScoresForThisPWM_agnostic < cutoff, SequenceScoresForThisPWM_agnostic, NA) #weak binding only
  
  TFScoresForThisVariant_agnostic =  computeTFScoreForThisVariant(seq, 
                                                                  TFScoresForThisVariant_agnostic, 
                                                                  SequenceScoresForThisPWM_agnostic, 
                                                                  SequenceScoresForThisPWM_WithCutoff_agnostic) 
                               
  
    #plus ------
  SequenceScoresForThisPWM_WithCutoff_plus = ifelse(SequenceScoresForThisPWM_plus >= cutoff, SequenceScoresForThisPWM_plus, NA) #strong binding only 
  SequenceScoresForThisPWM_plus = ifelse(SequenceScoresForThisPWM_plus < cutoff, SequenceScoresForThisPWM_plus, NA) #weak binding only
  
  TFScoresForThisVariant_plus = computeTFScoreForThisVariant(seq, 
                                                             TFScoresForThisVariant_plus, 
                                                            SequenceScoresForThisPWM_plus, 
                                                            SequenceScoresForThisPWM_WithCutoff_plus)

  
    #minus -----
  SequenceScoresForThisPWM_WithCutoff_minus = ifelse(SequenceScoresForThisPWM_minus >= cutoff, SequenceScoresForThisPWM_minus, NA) #strong binding only 
  SequenceScoresForThisPWM_minus = ifelse(SequenceScoresForThisPWM_minus < cutoff, SequenceScoresForThisPWM_minus, NA) #weak binding only
  
  TFScoresForThisVariant_minus = computeTFScoreForThisVariant(seq, 
                                                              TFScoresForThisVariant_minus,
                                                              SequenceScoresForThisPWM_minus, 
                                                              SequenceScoresForThisPWM_WithCutoff_minus)
  
                               
  
  } #for loop for refSeq, altSeq ends here
  } #for loop for PWMs for this variant ends here
  
  #compute different summary statistics for each TF for this PWM ----
    #maxNonZero = function(x){return(max(x[which(x!=0)]))}
  computeDifferentialSummaryStats = function(TFScoresForThisVariant){
    #maximum sequence scores for each PWM
  TFScoresForThisVariant$maxRefScore_noCutoff = unlist(lapply(TFScoresForThisVariant$SequenceScores_refSeq, max,na.rm = TRUE))
  TFScoresForThisVariant$maxAltScore_noCutoff = unlist(lapply(TFScoresForThisVariant$SequenceScores_altSeq, max,na.rm = TRUE))
    #average of sequence scores for each PWM
  TFScoresForThisVariant$meanRefScore_noCutoff = unlist(lapply(TFScoresForThisVariant$SequenceScores_refSeq, mean,na.rm = TRUE))
  TFScoresForThisVariant$meanAltScore_noCutoff = unlist(lapply(TFScoresForThisVariant$SequenceScores_altSeq, mean,na.rm = TRUE))
  
    #change the NA in max and mean alt scores to 0
  TFScoresForThisVariant$maxRefScore_noCutoff[which(is.infinite(TFScoresForThisVariant$maxRefScore_noCutoff))] = 0
  TFScoresForThisVariant$maxAltScore_noCutoff[which(is.infinite(TFScoresForThisVariant$maxAltScore_noCutoff))] = 0
  
  TFScoresForThisVariant$meanRefScore_noCutoff[which(is.na(TFScoresForThisVariant$meanRefScore_noCutoff))] = 0
  TFScoresForThisVariant$meanAltScore_noCutoff[which(is.na(TFScoresForThisVariant$meanAltScore_noCutoff))] = 0
  
    #Ref-Alt for each PWM - uncorrected
  TFScoresForThisVariant$differenceOfMax_noCutoff = TFScoresForThisVariant$maxAltScore_noCutoff - TFScoresForThisVariant$maxRefScore_noCutoff
  TFScoresForThisVariant$differenceOfMeans_noCutoff = TFScoresForThisVariant$meanAltScore_noCutoff - TFScoresForThisVariant$meanRefScore_noCutoff
  
    #maximum sequence scores for each PWM - with Cutoff
  TFScoresForThisVariant$maxRefScore_withCutoff = unlist(lapply(TFScoresForThisVariant$SequenceScores_refSeq_WithCutoff, max,na.rm = TRUE))
  TFScoresForThisVariant$maxAltScore_withCutoff = unlist(lapply(TFScoresForThisVariant$SequenceScores_altSeq_WithCutoff, max,na.rm = TRUE))
    #average sequence scores for each PWM - with Cutoff
  TFScoresForThisVariant$meanRefScore_withCutoff = unlist(lapply(TFScoresForThisVariant$SequenceScores_refSeq_WithCutoff, mean,na.rm = TRUE))
  TFScoresForThisVariant$meanAltScore_withCutoff = unlist(lapply(TFScoresForThisVariant$SequenceScores_altSeq_WithCutoff, mean,na.rm = TRUE))
    #change the NA in max and mean alt scores to 0
  TFScoresForThisVariant$maxRefScore_withCutoff[which(is.infinite(TFScoresForThisVariant$maxRefScore_withCutoff))] = 0
  TFScoresForThisVariant$maxAltScore_withCutoff[which(is.infinite(TFScoresForThisVariant$maxAltScore_withCutoff))] = 0
  
  TFScoresForThisVariant$meanRefScore_withCutoff[which(is.na(TFScoresForThisVariant$meanRefScore_withCutoff))] = 0
  TFScoresForThisVariant$meanAltScore_withCutoff[which(is.na(TFScoresForThisVariant$meanAltScore_withCutoff))] = 0
  
    #Ref - alt for each PWM - with cut off
  TFScoresForThisVariant$differenceOfMax_withCutoff = TFScoresForThisVariant$maxAltScore_withCutoff - TFScoresForThisVariant$maxRefScore_withCutoff
  TFScoresForThisVariant$differenceOfMeans_withCutoff = TFScoresForThisVariant$meanAltScore_withCutoff - TFScoresForThisVariant$meanRefScore_withCutoff
  TFScoresForThisVariant$deltaSignificantBindingSites = TFScoresForThisVariant$numberOfSignificantBindingSites_altSeq - TFScoresForThisVariant$numberOfSignificantBindingSites_refSeq
  
  return(TFScoresForThisVariant)
  }
  
  TFScoresForThisVariant_agnostic = computeDifferentialSummaryStats(TFScoresForThisVariant_agnostic)
  TFScoresForThisVariant_minus = computeDifferentialSummaryStats(TFScoresForThisVariant_minus)
  TFScoresForThisVariant_plus = computeDifferentialSummaryStats(TFScoresForThisVariant_plus)
  
  
  ##bind data to the variant row-----
    #plus ---
  variant_WithoutCutOff_plus[thisVariant,] = cbind(t(TFScoresForThisVariant_plus$differenceOfMax_noCutoff), 
                                              t(TFScoresForThisVariant_plus$differenceOfMeans_noCutoff))
  variant_WithCutOff_plus[thisVariant,] = cbind(t(TFScoresForThisVariant_plus$differenceOfMax_withCutoff), 
                                           t(TFScoresForThisVariant_plus$differenceOfMeans_withCutoff), 
                                           t(TFScoresForThisVariant_plus$deltaSignificantBindingSites))
  
    #minus ---
  variant_WithoutCutOff_minus[thisVariant,] = cbind(t(TFScoresForThisVariant_minus$differenceOfMax_noCutoff), 
                                                   t(TFScoresForThisVariant_minus$differenceOfMeans_noCutoff))
  variant_WithCutOff_minus[thisVariant,] = cbind(t(TFScoresForThisVariant_minus$differenceOfMax_withCutoff), 
                                                t(TFScoresForThisVariant_minus$differenceOfMeans_withCutoff), 
                                                t(TFScoresForThisVariant_minus$deltaSignificantBindingSites))
    #agnostic ---
  variant_WithoutCutOff_agnostic[thisVariant,] = cbind(t(TFScoresForThisVariant_agnostic$differenceOfMax_noCutoff), 
                                                   t(TFScoresForThisVariant_agnostic$differenceOfMeans_noCutoff))
  variant_WithCutOff_agnostic[thisVariant,] = cbind(t(TFScoresForThisVariant_agnostic$differenceOfMax_withCutoff), 
                                                t(TFScoresForThisVariant_agnostic$differenceOfMeans_withCutoff), 
                                                t(TFScoresForThisVariant_agnostic$deltaSignificantBindingSites))
  
  
  ##separate data frame with just the sequence scores for each TF for each variant-----
  consolidateTFSequenceScoresForAllVariants = function(TFScoresForThisVariant, TFScoresForVariants){
  sequenceScores = select(TFScoresForThisVariant, contains("SequenceScores"))
  colnames(sequenceScores) = paste(colnames(sequenceScores),rownames(variantsDF)[thisVariant], sep = "." )
  TFScoresForVariants = cbind(TFScoresForVariants, sequenceScores)
  return(TFScoresForVariants)
  }
  
  TFScoresForVariants_agnostic = consolidateTFSequenceScoresForAllVariants(TFScoresForThisVariant_agnostic, TFScoresForVariants_agnostic)
  TFScoresForVariants_plus = consolidateTFSequenceScoresForAllVariants(TFScoresForThisVariant_plus, TFScoresForVariants_plus)
  TFScoresForVariants_minus = consolidateTFSequenceScoresForAllVariants(TFScoresForThisVariant_minus, TFScoresForVariants_minus)
  } #for loop for all variants ends here
  
  #get number of TFs whose TFBS is affected by variant
  otherInfo = resultsAndAnnotations[,1:44]
  otherInfo$phastConsScore = resultsAndAnnotations$phastConScore
  otherInfo$GCContentOfOligo = resultsAndAnnotations$GCContentOfOligo
  otherInfo$GCFoldChange = resultsAndAnnotations$GCFoldChange
  

  deltaSignificant = select(variant_WithCutOff_plus, contains("deltaSignificant")) #define a dataframe with only the 196 columns of 'deltaSignificantBindingSites'
  variant_WithCutOff_plus$numberOfTFBSGainOrLost = rowSums(deltaSignificant!=0)
  
  variantData_WithoutCutOff_plus = cbind(resultsAndAnnotations[,c("variantID", "logFC")], variant_WithoutCutOff_plus)
  variantData_WithCutOff_plus = cbind(resultsAndAnnotations[, c("variantID", "logFC")], variant_WithCutOff_plus)
  
  deltaSignificant = select(variant_WithCutOff_minus, contains("deltaSignificant")) #define a dataframe with only the 196 columns of 'deltaSignificantBindingSites'
  variant_WithCutOff_minus$numberOfTFBSGainOrLost = rowSums(deltaSignificant!=0)
  
  variantData_WithoutCutOff_minus = cbind(resultsAndAnnotations[,c("variantID", "logFC")], variant_WithoutCutOff_minus)
  variantData_WithCutOff_minus = cbind(resultsAndAnnotations[, c("variantID", "logFC")], variant_WithCutOff_minus)
  
  deltaSignificant = select(variant_WithCutOff_agnostic, contains("deltaSignificant")) #define a dataframe with only the 196 columns of 'deltaSignificantBindingSites'
  variant_WithCutOff_agnostic$numberOfTFBSGainOrLost = rowSums(deltaSignificant!=0)
  
  variantData_WithoutCutOff_agnostic = cbind(resultsAndAnnotations[,c("variantID", "logFC")], variant_WithoutCutOff_agnostic)
  variantData_WithCutOff_agnostic = cbind(resultsAndAnnotations[, c("variantID", "logFC")], variant_WithCutOff_agnostic)
  
  
  #save --------
  
  save(variantData_WithoutCutOff_plus, variantData_WithCutOff_plus, otherInfo, 
       file = paste0(dir, "data/PromoterDataWithNewTFScores_PlusStrand_ScerTF.RData"))
  
  #individual sequence scores for all variants for the 196 TF in ScerTF db - with and without cutoff
  save(TFScoresForVariants_plus, 
       file = paste0(dir, "results/TFSequenceScoresForAllVariants_PlusStrand_ScerTF.rda"))
  
  
  save(variantData_WithoutCutOff_minus, variantData_WithCutOff_minus, otherInfo, 
       file = paste0(dir, "results/PromoterDataWithNewTFScores_MinusStrand_ScerTF.RData"))
  
  #individual sequence scores for all variants for the 196 TF in ScerTF db - with and without cutoff
  save(TFScoresForVariants_minus, 
       file = paste0(dir, "results/TFSequenceScoresForAllVariants_MinusStrand_ScerTF.rda"))
  
  
  save(variantData_WithoutCutOff_agnostic, variantData_WithCutOff_agnostic, otherInfo, 
       file = paste0(dir, "results/PromoterDataWithNewTFScores_StrandAgnostic_ScerTF.RData"))
  
  #individual sequence scores for all variants for the 196 TF in ScerTF db - with and without cutoff
  save(TFScoresForVariants_agnostic, 
       file = paste0(dir, "results/TFSequenceScoresForAllVariants_StrandAgnostic_ScerTF.rda"))
  

  
#join withCutOff and without CutOff dfs with the other variantFeatures and save
  
  #all variant data
  
  #TF-PWMs from ScerTF db
    #load the ScerTF gene list with the systematic names - sourced from http://yeastgenome.org
    ScerTF_geneList = read.csv(paste0(dir,"data/ScerTF_geneList.csv"), sep="", stringsAsFactors=FALSE)
  
    #add systematic name to the TF-PWM matrix
    for(i in 1:nrow(TFMatrix)){
      TFName = TFMatrix$TF[i]
      TFSystematicName = ScerTF_geneList$secondaryIdentifier[which(ScerTF_geneList$input == TFName)]
      TFMatrix$systematicName[i] = TFSystematicName
    }
  
    #save the file
    save(TFMatrix, file = paste0(dir, "results/TF_PWMs_ScerTF.rda"))
    
 
 
  
  