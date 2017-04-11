## Load packages
library("DESeq2")
library("tximport")

## Prepare vectors with sample information
sampleNames = sort(c("ERR458500","ERR458514","ERR458528","ERR458885","ERR458899","ERR458493","ERR458507","ERR458521","ERR458878","ERR458892"))
sampleTypes = c("WT","SNF2","SNF2","SNF2","SNF2","SNF2","WT","WT","WT","WT")
sampleFiles = c(paste(sampleNames, '/abundance.tsv', sep=''))

## Read in files
txDat = tximport(sampleFiles, type="kallisto", txOut=TRUE)
coldata = data.frame(condition=sampleTypes)
rownames(coldata) = sampleNames

## DESeq2 object
dds = DESeqDataSetFromTximport(txDat, colData=coldata, design=~condition)
dex = DESeq(dds)
res = results(dex)

## Plot
plotMA(res)
sizeFactorEst = estimateSizeFactors(dex)
dispEst = estimateDispersions(dex)
plotDispEsts(dex, 1)

## Get FDR
fdr = res$padj
numGeneBelow5Percent = sum(fdr<0.05, na.rm=TRUE)
numGeneBelow5Percent
