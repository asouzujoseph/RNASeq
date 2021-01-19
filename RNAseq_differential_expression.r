# load necessary packages
library("DESeq2")
library(dplyr)
library(scales)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)


getwd()

# include path to the gene counts file
matrixFile = path.expand("/mnt/nfs/data/RNA-seq/MCF7/all.gene.counts")

# create a matrix of the gene counts file
countData <- as.matrix(read.csv(matrixFile, sep="\t", row.names="Geneid"))
head(countData)

# modify the column names and remove the .bam suffix
colnames(countData) <- c("NS1", "NS2", "S1", "S2")
head(countData)

# create the conditions which are stimulated (S) vs non-stimulated(NS)
condition <- factor(c("NS", "NS", "S", "S"))
print(condition)

# create a data frame to be used in generating the DESEQ object
colData <- data.frame(sampleName = colnames(countData),
                      condition = condition)
colData

# Create a DESEQ2 object
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ condition)
# colData = a data frame with columns as the variables known about the samples (such as conditions)and rownames as unique sample names
# countData = a matrix of the actual gene count values

# Set the column names to be same as the treatments. we need to do this to ensure the column names are correct.
colnames(dds) <- colnames(countData)

dds$condition

dim(dds)

# filter out genes with no reads
keep_genes <- rowSums(counts(dds)) > 0
dds <- dds[keep_genes,]
dim(dds)


# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)
options(repr.plot.width=10, repr.plot.height=2)
plotPCA(rld, intgroup="condition")

options(repr.plot.width=10, repr.plot.height=10)
plot(data.frame(log10(counts(dds))), cex=0.1, col=alpha("black", 0.4))


deseq2VST <- vst(dds)
sampleDists <- dist(t(assay(deseq2VST)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- c("grey", "pink")
pheatmap(sampleDistMatrix,clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         cluster_rows =T, cluster_cols = T,
         col = colors)

dE <- DESeq(dds)

# The MA plot is used to visualize the significantly and differentially expressed genes which are shown in red
options(repr.plot.width=5, repr.plot.height=4)
plotMA(dE,ylim=c(-3,3))

res <- results(dE)
res <- res[order(res$padj),]
res.sign <- res[(!is.na(res$padj)) & (res$padj < 0.05),] #select rows where padj is not NA and less than 0.05
nrow(res.sign) # number of significantly expressed genes


# Histogram showing the distribution of significant genes (p-adj < 0.05)
hist(res$padj, 
     col="black", border="white", xlab="", ylab="", main="frequencies of adj. p-values\n(all genes)")

par(mfrow=c(1,2))
options(repr.plot.width=10, repr.plot.height=4)
plotCounts(dds, gene="BBC3", normalized = TRUE) # UPREGULATED BY p53 in stimulated cells
plotCounts(dds, gene="TIMELESS", normalized = TRUE) # DOWNREGULATED BY p53 in Stimulated cells

setwd('/mnt/storage/r0605462/jupyternotebooks/RNASeq/')

write.table(res, "deseq.results.tsv", sep="\t", col.names=NA, quote=FALSE)

sessionInfo()
