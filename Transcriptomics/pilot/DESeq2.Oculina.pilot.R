# necessary libraries
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(pheatmap)
library(RColorBrewer)
library(limma)
library(ggforce)
library(ggVennDiagram)
library(NbClust)
library(ComplexHeatmap)
library(edgeR)
library(DEGreport)
library(fgsea)
library(purrr)
library(tibble)
library(reshape2)
library(patchwork)
library(stringr)
library(ggpattern)


# setting working directory 
setwd("/home/gospozha/haifa/O.patagonica/rna/analysis/de")

#### Preparing necessary files ####

# list files with gene counts for each sample
dir = "/home/gospozha/haifa/O.patagonica/rna/analysis/de/"
files = list.files(paste0(dir, "count.tables"), "*ReadsPerGene.out.tab", full.names = T)
countData = data.frame(fread(files[1]))[c(1,2)]

# looping and reading the 2nd column from the remaining files
for(i in 2:length(files)) {
  countData = cbind(countData, data.frame(fread(files[i]))[2])
}

# skipping the first 4 lines, since count data starts on the 5th line
countData = countData[c(5:nrow(countData)),]

# renaming columns as sample names
colnames(countData) = c("GeneID", gsub(paste0(dir,"count.tables/"), "", files))
colnames(countData) = gsub(".fastq.gzReadsPerGene.out.tab", "", colnames(countData))
rownames(countData) = countData$GeneID
countData = countData[,c(2:ncol(countData))]
names <- colnames(countData)

# writing count matrix to a file
write.csv(countData, file="CountMatrix.csv")  

# reading count matrix from a file
countData  <- read.csv2('CountMatrix.csv', header=TRUE, row.names=1, sep=',', check.names = F)
# reading metadata file
MetaData <- read.csv2('Metadata', header=TRUE, sep="\t")

MetaData$condition <- as.factor(MetaData$condition)
MetaData$origin <- as.factor(MetaData$sample)


#### Initial quality check ####

# Convert counts to DGEList
dge <- DGEList(counts = countData)
# remove low counts
smallestGroupSize <- 2
keep <- rowSums(dge$counts >= 10) >= smallestGroupSize
dge <- dge[keep,]

# Calculate FPM (Fragments Per Million)
fpm_values <- cpm(dge, normalized.lib.sizes = TRUE)  # edgeR's CPM is equivalent to FPM

# Convert to long format for plotting
fpm_df <- as.data.frame(fpm_values) %>%
  tibble::rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "FPM") %>%
  left_join(MetaData, by = c("Sample" = "id"))  # Merge with metadata

ggplot(fpm_df, aes(x = FPM, color = condition)) +
  geom_density(alpha = 0.3) +
  scale_x_log10() +
  theme_minimal() +
  labs(title="Density Plot of FPM Values per Condition",
       x="FPM (log10 scaled)")

# statistical comparison
anova_res <- aov(FPM ~ condition, data = fpm_df)
summary(anova_res)
TukeyHSD(anova_res)
# they are the same

#### DESeq2 model ####
# creating DESeq2 object 
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = MetaData,
                              design = ~ condition)

# pre filtering
smallestGroupSize <- 2
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dim(dds) # 2702 genes have left

# running a model
dds <- DESeq(dds)
res <- results(dds)

# Plotting histograms of p-values
hist(res$pvalue, breaks=50, col="skyblue", main="~ condition",
     xlab="p-value", xlim=c(0,1), ylim=c(0, max(table(cut(res$pvalue, breaks=50)))))

# saving a DESeq2 model to an R object
saveRDS(dds, file = "dds_site_rrna_rin_condition.rds")
#dds <- readRDS(file = "../dds_site_rrna_rin_condition.rds")

#### PCA and sample distances using rlog ####

# estimating size factors to determine if it's better to use rlog
# to transform our data. rlog is more robust if size factors differ a lot
SF <- estimateSizeFactors(dds) 
print(sizeFactors(SF))

# the same using rlog transformation
rlog <- rlog(dds)
saveRDS(rlog, file = "rlog_site_rrna_rin_condition.rds")
#rlog <- readRDS(file = "rlog_site_rrna_rin_condition.rds")
# since vst does not remove variation that can be associated with covariates, 
# we manually remove the effect of covariates to be able to visualize it on PCA
mat <- assay(rlog)
mm <- model.matrix(~condition, colData(rlog))
treatment.design <- mm[,1:6]
batch.design <- mm[,-(1:6)]
mat <- limma::removeBatchEffect(mat, covariates=batch.design, design=treatment.design)
assay(rlog) <- mat
# PCA plot
pcaData <- plotPCA(rlog, intgroup=c("condition"), ntop = 500, returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
#pdf("PCA.full.pdf",width=7)
pca<-ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  ggtitle("PCA of gene counts") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_bw()
pca
ggsave("pca.jpg", pca, width = 6.5, height = 6)
norm.counts <- assay(rlog)
write.csv(norm.counts, file="./rlog.counts.csv")

#dev.off()

# sample distances
sampleDists <- dist(t(assay(rlog)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rlog$condition)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
#pdf("Dist.all.pdf",width=7)
dist <- pheatmap(sampleDistMatrix,
                 clustering_distance_rows=sampleDists,
                 clustering_distance_cols=sampleDists,
                 col=colors)
dist
ggsave("dist.jpg", dist, width = 6, height = 6)
#dev.off()


res <- results(dds, contrast=c("condition","10","30"), alpha = 0.05)
summary(res)
