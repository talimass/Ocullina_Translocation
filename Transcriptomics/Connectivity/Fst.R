library(adegenet)
library(StAMPP)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(SNPRelate)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(ggforce)
library(dplyr)
library(cowplot)

setwd("/home/gospozha/haifa/hiba/op_align_new/snp/pruned")


# PLINK to GDS
snpgdsBED2GDS("Ocupat_vcf_pruned_filt_0.3.bed", 
              "Ocupat_vcf_pruned_filt_0.3.fam", 
              "Ocupat_vcf_pruned_filt_0.3.bim", 
              "temp.gds", cvt.chr="char")

genofile <- snpgdsOpen("temp.gds")
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
geno <- snpgdsGetGeno(genofile) # Rows = Samples, Cols = SNPs
snpgdsClose(genofile)

# Load Metadata 
meta <- read.csv("pop.csv")

meta <- meta[match(sample.id, meta$id), ]

# Genotype conversion 
# convert 0,1,2 directly to AA, AB, BB
genotype_matrix <- matrix(NA, nrow=nrow(geno), ncol=ncol(geno))
genotype_matrix[geno == 0] <- "AA"
genotype_matrix[geno == 1] <- "AB"
genotype_matrix[geno == 2] <- "BB"

# StAMPP 'r' format dataframe
# StAMPP format 'r' requires: Sample, Pop, Ploidy, Format, SNP1, SNP2...
st_df <- data.frame(
  Sample = sample.id,
  Pop = meta$condition,  
  Ploidy = 2,
  Format = "BiA",
  stringsAsFactors = FALSE
)

# combine metadata with genotype Matrix
final_df <- cbind(st_df, genotype_matrix)

# filter out NAs in Pop if any
final_df <- final_df[!is.na(final_df$Pop), ]

# convert to StAMPP format
genotype.st <- stamppConvert(final_df, "r")

# run Fst
genotype.fst <- stamppFst(genotype.st, nboots = 1000, percent = 95, nclusters = 4)

print(genotype.fst$Fsts)
print(genotype.fst$Pvalues)


# Calculate Nei's Distance 
neis_dist <- stamppNeisD(genotype.st, pop = TRUE)


#### Fst plot ####

# prepare the data from stampp object
fst_mat <- genotype.fst[["Fsts"]]
p_mat <- genotype.fst[["Pvalues"]]

# ensure P-values are symmetric 
p_mat[upper.tri(p_mat)] <- t(p_mat)[upper.tri(p_mat)]

# combined matrix for display
combined_mat <- fst_mat
#  upper triangle with the P-values
combined_mat[upper.tri(combined_mat)] <- p_mat[upper.tri(p_mat)]

# label matrix
n <- nrow(combined_mat)
label_matrix <- matrix("", nrow = n, ncol = n)

for(i in 1:n) {
  for(j in 1:n) {
    if(i > j) {
      label_matrix[i, j] <- sprintf("%.3f", fst_mat[i, j]) 
    } else if (i < j) {
      label_matrix[i, j] <- paste0("p=", sprintf("%.3f", p_mat[i, j]))
    }
  }
}

# color 
col_fun <- colorRamp2(
  breaks = pretty(range(fst_mat, na.rm = TRUE), n = 5), 
  colors = colorRampPalette(c("#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695"))(length(breaks))
)


ht <- Heatmap(
  fst_mat, 
  name = "Fst",
  col = col_fun,
  na_col = "white",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  cell_fun = function(j, i, x, y, width, height, fill) {

    if(i < j) {
      grid.rect(x, y, width, height, gp = gpar(fill = "white", col = NA))
    }
    
    label <- label_matrix[i, j]
    if (label != "") {
      text_col <- ifelse(i < j, "grey40", "black")
      grid.text(label, x, y, gp = gpar(fontsize = 10, col = text_col))
    }
    # diagonal slash
    if (i == j) {
      grid.lines(x = unit.c(x - 0.5 * width, x + 0.5 * width),
                 y = unit.c(y + 0.5 * height, y - 0.5 * height),
                 gp = gpar(col = "grey80"))
    }
  }
)

draw(ht)

#### PCA plot ####
genotype.st  -> genotype.st2
geno.nas <- genotype.st2[, !(names(genotype.st2) %in% c("Sample", "Pop", "pop.num", "ploidy", "format"))]
pop.vector <- as.factor(genotype.st2[[2]])
convert_to_genind_format <- function(x) {
  ifelse(is.na(x), NA,
         ifelse(x == 1, "11",
                ifelse(x == 0.5, "12",
                       ifelse(x == 0, "22", NA))))
}

geno.char <- as.data.frame(lapply(geno.nas, convert_to_genind_format), stringsAsFactors = FALSE)
geno.obj <- df2genind(geno.char, pop = pop.vector, ploidy = 2, NA.char = NA, sep = "")
dist.mtx <- dist(geno.obj)
pca <- dudi.pca(tab(geno.obj, NA.method = "mean"), scannf = F, nf = 2)
s.class(pca$li, fac = pop.vector, col = rainbow(length(unique(pop.vector))))


# Extract individual scores (coordinates)
scores <- as.data.frame(pca$li)
scores$Pop <- pop.vector
scores$SampleID <- rownames(scores)

# Variance explained
eig <- 100 * pca$eig / sum(pca$eig)

# Calculate centroids per population
centroids <- scores %>%
  group_by(Pop) %>%
  dplyr::summarise(
    Axis1 = mean(Axis1),
    Axis2 = mean(Axis2)
  )

# Add centroid coords to scores
scores <- scores %>%
  left_join(centroids, by = "Pop", suffix = c("", ".centroid"))

# Plot
pca_plot2 <- ggplot(scores, aes(x = Axis1, y = Axis2, color = Pop)) +
  # segments from sample to centroid
  geom_segment(aes(xend = Axis1.centroid, yend = Axis2.centroid),
               alpha = 0.4, linewidth = 0.3) +
  # centroid labels (drawn first, so underneath points)
  geom_text(data = centroids, aes(x = Axis1, y = Axis2, label = Pop, color = Pop),
            , size = 6, show.legend = FALSE, alpha = 0.6) +
  # points for samples (drawn after, so on top)
  geom_point(size = 3) +
  # centroid markers (optional cross shape)
  #geom_point(data = centroids, aes(x = Axis1, y = Axis2, color = Pop),
  #           size = 1, shape = 4, stroke = 1.5) +
  labs(
    x = paste0("PC1 (", round(eig[1], 1), "%)"),
    y = paste0("PC2 (", round(eig[2], 1), "%)")
  ) +
  theme_minimal() +
  theme(legend.position = "none")

pca_plot2

#### save pics ####

ht_grob <- grid.grabExpr(draw(ht))

plotA <- ggdraw(ht_grob)

plotB <- pca_plot2 

final_combined_plot <- plot_grid(
  plotA, plotB, 
  ncol = 2, 
  labels = "AUTO",       
  label_size = 14,       
  rel_widths = c(1, 1) 
)

# Export the result
save_plot("Fst_PCA_Combined_Plot.pdf", final_combined_plot, base_width = 10, base_height = 4)
save_plot("Fst_PCA_Combined_Plot.png", final_combined_plot, base_width = 10, base_height = 4, dpi = 300)

plot(final_combined_plot)


save.image("050426.Rdata")