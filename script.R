library(SummarizedExperiment)
library(DESeq2)
library(ggplot2)
library(pcaExplorer)
library(limma)
library(Glimma)
library(edgeR)
library(tidyverse)
library(goseq)
setwd("~/Desktop/Genomic Data Science/Capstone/Task 5/")

############
# Setup:  Download data, remove rows of zeros (they're uninformative), and create a SummarizedExperiment object
# Variables:  counts - Raw data from files, imported as a data frame
#             Counts - A matrix version of counts. This was needed because the SummarizedExperiment constructor expected one
############
counts <- read.table("geneCounts.txt", sep = "\t", header = TRUE)
rownames(counts) <- counts[,1]
counts <- counts[,2:ncol(counts)]
Counts <- matrix(as.numeric(as.character(unlist(counts))),nrow=nrow(counts))
rownames(Counts) <- rownames(counts)

############
# Pre-processing: We'll use CPM to normalize counts to account for sequencing depth, and then remove lowly expressed
#                 genes. 
# Variables:  (l)cpm_Counts - Tranformed Counts matrix.
#             a_df - An intermediate object used for visualizing
#             filteredCounts - The Counts matrix, with lowly expressed genes filtered out
############
# First, we'll visualize the distribution of counts in each sample to
# see if we have to do anything about sequencing depth.
a_df <- as.data.frame(log(Counts))
names(a_df) <- names(counts)
a_df %>% gather(Sample, logCount) %>% 
  ggplot(aes(Sample, logCount)) + 
  geom_boxplot() + 
  ggtitle("Expression Levels of Samples (Raw Data)") +
  labs(y = "Log of Gene Read Counts")

# So one of them looks under-represented. We'll transform with a CPM
# transformation.
lcpm_Counts <- cpm(Counts, log = T)
cpm_Counts <- cpm(Counts)
a_df <- as.data.frame(lcpm_Counts)
names(a_df) <- names(counts)
a_df %>% gather(Sample, logCount) %>% 
  ggplot(aes(Sample, logCount)) + 
  geom_boxplot() + 
  ggtitle("Expression Levels of Samples after CPM transformation") +
  labs(y = "Log of Gene Read Counts")
# This now looks like the counts have been normalized.

# Now we'll see the density of counts to see how many lowly-expressed
# genes there are
a_df <- as.data.frame(lcpm_Counts)
names(a_df) <- names(counts)
a_df %>% gather(Sample, logCount) %>% 
  ggplot(aes(x = logCount, color = Sample)) + 
  geom_density() + 
  ggtitle("Distribution of Read Counts per Sample") +
  labs(y = "Log of Gene Read Counts")
# It looks like there are a lot of lowly-expressed genes.

# We must remove these lowly expressed genes.
filteredCounts <- Counts[which(filterByExpr(Counts)),]
a_df <- as.data.frame(cpm(filteredCounts, log = T))
names(a_df) <- names(counts)
a_df %>% gather(Sample, logCount) %>% 
  ggplot(aes(x = logCount, color = Sample)) + 
  geom_density() + 
  ggtitle("Distribution of Read Counts per Sample (Filtered)") +
  labs(y = "Log of Gene Read Counts")
# With those genes removed, we're left with a much smaller data frame.
# Also, it looks like the distribution of read counts in each sample
# is roughly similar, so we don't need to account for external factors.

############
# Visualization: boxplot of expression levels + PCA plot of first two principal components
# Variables:  dds(2) - A DESeqDataSet object; intermediate set used to give input to the vst() function
#             vsd(2) - A vst-transformed version of dds used for purposes of input for the plotpca() function
#             p_data - Matrix of phenotype data from the experiment. To be used to create summarizedExperiment object
#             pData - A dataframe-version of p_data. Used because SummarizedExperiment requires a data frame
#             seo - The summarized experiment object. Consolidation of (filtered) count data and p_data
############
# First, create the seo object to feed into DESeq.
p_data <- read.table("pData.txt", header = TRUE)
rownames(p_data) = p_data[,1]
pData <- DataFrame(p_data)
seo <- SummarizedExperiment(assays = list(counts=filteredCounts), colData = pData)

# Create DESeq objects to visualize using PCA
dds <- DESeqDataSet(seo, design = ~ Group)
vsd <- vst(dds, blind = F)
pcaplot(vsd, intgroup = "Group")

dds2 <- DESeqDataSet(seo, design = ~ Sex)
vsd2 <- vst(dds, blind = F)
pcaplot(vsd2, intgroup = "Sex")
# By the first two principal components, the samples seem to cluster nicely by Group, but not by Sex.

############
# Differential Expression Analysis: We'll determine and visualize the differentially epxressed genes with a volcano plot
# Variables:  dds - A DESeqDataSet object, with differential expression analysis done
#             res - A data frame with the log fold changes and associated p-values for the null hypothesis that the gene isn't differentially epxressed.
#                   Also contains additional information
#             resSig - Filtered dataframe from res. Only includes those genes that have an adjusted p-value of < 0.05
#             a_df - Intermediate dataframe used to visualize the volcano plot
############
# We'll use DESeq to find differentially expressed genes.
dds <- DESeq(dds)
res <- results(dds)
resSig <- subset(res, padj < 0.05)
resSig <- resSig[order(resSig$padj),]

# Now we'll visualize the data in a volcano plot.
a_df <- as.data.frame(resSig[,c("log2FoldChange", "padj")])
a_df[,2] <- -log(a_df[,2])
a_df[,3] <- sapply(c(1:nrow(a_df)), function(x) if (x <= 0.05*nrow(a_df)) return("Very Significant") else return("Not Very Significant"))
names(a_df)[3] <- "Significant"
ggplot(a_df, aes(x = log2FoldChange, y = padj, color = Significant)) + 
  geom_point() + 
  ggtitle("Fold Change of Significant Differentially Expressed Genes") +
  labs(y = "Adjusted p-value")


############
# GO term annotation: Now that we have a list of significantly, differentially-expressed genes, we'll annnotate them with GO terms.
# Variables:  genes - Names of all the significant genet
#             go - A list, where each item is a list of GO terms associated with a gene in gnens
#             goList - The union of all the GO terms in go
############
genes <- rownames(resSig)[c(1:as.integer(0.05 * nrow(resSig)))]
go <- getgo(genes, "hg38", "knownGene")
goList <- unique(unlist(go))
write(goList, file = "goTerms.txt")
# Now that we have a list of GO terms, we can perform an enrichment analysis on the list of significantly differentially epxressed genes.
