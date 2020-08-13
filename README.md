# promoterVariants
### Code, objects and annotations for Cheung, Renganaath, et al., 2020

The code is broken into the six major areas listed below.

All code written by Frank Albert, except for code in area 5 and TFBS annotation code in area 4, which was written by Kaushik Renganaath.

For each area, there is one folder with code, and one folder with various objects, annotations, etc, that get used by the code. In some cases, objects get used by several analysis areas (such as the oligo design file, gene annotations, barcode counts). Be sure to look through all object folders and adjust your paths as necessary. Some of the annotation files are gzipped – need to unzip before use with the R code.

Throughout, the code files are numbered in the order in which they should be used. Sometimes a later file won't work if it is run before an earlier file because it needs intermediate objects saved by some of the earlier files. We provided the most important of these intermediate objects, but some minor files may need to be created using the code here. Some huge files are also missing from this code repository (e.g. barcode level counts, not just those to designed oligos). These are available at the GEO repository listed in the paper or at https://figshare.com/s/1c23d927e17fc203ac3b

A general note:
Many "save" or "write.table" commands are commented out. This can prevent accidental overwriting of files. Check for these commented lines – you want to actually run these, once you've confirmed that no important files will be overwritten.

## The analysis areas are:

### 1. Oligo library design
a) TSS_preparations.R: Make TSS annotations from Pelechano 2013 data

b) The actual design code

c) reverse-complement oligos with more A than T nucleotides

d) analyze the design to get descriptives for the paper


### 2. Annotation runs
Map barcodes to oligos, count the combinations, and map them to the MPRA design. There is one code file for TSS and one for Upstream. For each library, the code produces three files that will be needed in downstream analyses. Four of these six files are big (TSS library: R_bcCountsAssigned_160629.RData and R_bcCountsTopOligoOnly_160629.RData; Upstream library: R_bcCountsAssigned_UpStream.RData and R_bcCountsTopOligoOnly_161203_UpStream.RData). They are available at https://figshare.com/s/1c23d927e17fc203ac3b.


### 3. Test for causal variants and do basic analyses on them
"01_countAndcombineAllSamples": does just that. For every replicate sample, count the barcodes, aggregate them into a big table. This uses a lot of RAM and runs for hours. The key output of this is "R_countsMappedToDesignedOligos_180616.RData", which is available at https://figshare.com/s/1c23d927e17fc203ac3b.

The three files starting "0x_mpralm..." contain the statistical analyses of single variants and epistasis.

"06_aggregateVariantResults.R" combines the results across the TSS and Upstream library, does reproducibility of single variants, and makes the big volcano plot.

"07_variantsPerGene_and_eQTLs.R" compares the single variant effects to local eQTLs.


### 4. variant annotation
Make and gather various features. While all code needed to recapitulate how these features were assembled is here, it is messy due to having been written over multiple years, with multiple updates to various intermediate objects. It may be more efficient to just work with the final product ("R_resultsAndAnnotations_200719.RData"). This file has all non-TFBS annotations and variant results, which then got fed into feature association testing.


### 5. feature analyses & prediction
Single features and the various multiple regression models

### 6. Application of the deBoer 2020 model to our data
The model is described here: https://www.nature.com/articles/s41587-019-0315-8

All required code and data are in a dedicated folder. We made minor modifications to the code from the paper above.
