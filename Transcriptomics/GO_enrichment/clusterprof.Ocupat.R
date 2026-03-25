library(goseq)
library(tidyverse)
library(GSEABase)               #BiocManager::install("GSEABase")
library(data.table)
library(ggplot2)
library(cowplot)                #install.packages("cowplot")
library(patchwork)
library(dplyr)
library(tidyr)
library(clusterProfiler)
library(ontologyIndex)
library(GOSemSim)
library(simplifyEnrichment)
library(org.Hs.eg.db)
library(rlang)
library(readr)
library(GO.db)

#save.image(file="220126.GO.planulae.RData") 
#load("/home/gospozha/haifa/hiba/pl_align/clusterprofiler/220126.GO.planulae.RData")

# setting working directory 
setwd("/home/gospozha/haifa/hiba/op_align_new/clusterprofiler/")

term2gene <- read_tsv("term2gene.tsv")

term2name <- AnnotationDbi::select(
  GO.db,
  keys = unique(term2gene$term),
  columns = c("TERM"),
  keytype = "GOID"
) %>%
  rename(term = GOID, name = TERM) %>%
  distinct()

write_tsv(x = term2name, file = "term2name_GO.tsv")


term2name <- read_tsv("term2name_GO.tsv")


term2gene <- term2gene %>%
  filter(term %in% term2name$term)

#### reading count matrix from a file ####
countData  <- read.csv2('../CountMatrix.csv', header=TRUE, row.names=1, sep=',', check.names = F)
# reading metadata file
MetaData <- read.csv2('../Metadata1.csv', header=TRUE, sep=",")

# sample names in both objects
samples_meta  <- MetaData$id
samples_count <- colnames(countData)
# find common samples
common_samples <- intersect(samples_meta, samples_count)
# subset and reorder count matrix
countData<- countData[, common_samples]
# reorder metadata to match countData
MetaData <- MetaData[match(common_samples, MetaData$id), ]
# must be TRUE
all(colnames(countData) == MetaData$id)

MetaData$condition <- as.factor(MetaData$condition)
countData$geneID <- rownames(countData)
smallestGroupSize <- 4
keep <- rowSums(countData >= 10) >= smallestGroupSize
countData <- countData[keep,]

background_genes <- countData %>%
  dplyr::select("geneID") %>%
  unlist() %>%
  as.vector()

  
  #### Run enrichment separately ####
file_path <- ("../de_genes_30v1.csv")
interesting_set <- read_csv(file_path, show_col_types = FALSE) %>%               # if you want only significant genes
  pull(gene_id) %>%                        # <-- THIS is your LOC column
  unique() %>%
  na.omit()

interesting_set_30 <- read_csv(file_path, show_col_types = FALSE) %>%     
  filter(padj < 0.05, log2FoldChange > 0) %>%
  pull(gene_id) %>%                       
  unique() %>%
  na.omit()

interesting_set_1 <- read_csv(file_path, show_col_types = FALSE) %>%     
  filter(padj < 0.05, log2FoldChange < 0) %>%
  pull(gene_id) %>%                       
  unique() %>%
  na.omit()

enrichment_30 <- enricher(interesting_set_30,
                         TERM2GENE = term2gene,
                         TERM2NAME = term2name,
                         pvalueCutoff = 0.05,
                         universe = background_genes,
                         pAdjustMethod = "fdr",
                         qvalueCutoff = 0.2)
  
  # Save enrichment results
write_csv(enrichment_30@result%>%
            filter(!grepl("^GO:", Description)), "enrichment_results.30v1.csv")

enrichment_1 <- enricher(interesting_set_1,
                          TERM2GENE = term2gene,
                          TERM2NAME = term2name,
                          pvalueCutoff = 0.05,
                          universe = background_genes,
                          pAdjustMethod = "fdr",
                          qvalueCutoff = 0.2)

# Save enrichment results
write_csv(enrichment_1@result%>%
            filter(!grepl("^GO:", Description)), "enrichment_results.1v30.csv")
# not significant

p <- dotplot(enrichment_30,
                            x = "geneRatio",
                            color = "p.adjust",
                            orderBy = "x",
                            showCategory = 100,
                            font.size = 7) 

ggsave(paste0("enrichment_dotplot_30.pdf"), p, width = 4, height = 9)

p <- dotplot(enrichment_1,
             x = "geneRatio",
             color = "p.adjust",
             orderBy = "x",
             showCategory = 100,
             font.size = 7) 

ggsave(paste0("enrichment_dotplot_1.pdf"), p, width = 4, height = 9)


#### reading count matrix from a file 25vs10vs5 ####
countData  <- read.csv2('../CountMatrix.csv', header=TRUE, row.names=1, sep=',', check.names = F)
# reading metadata file
MetaData <- read.csv2('../Metadata2.csv', header=TRUE, sep=",")

# sample names in both objects
samples_meta  <- MetaData$id
samples_count <- colnames(countData)
# find common samples
common_samples <- intersect(samples_meta, samples_count)
# subset and reorder count matrix
countData<- countData[, common_samples]
# reorder metadata to match countData
MetaData <- MetaData[match(common_samples, MetaData$id), ]
# must be TRUE
all(colnames(countData) == MetaData$id)

MetaData$condition <- as.factor(MetaData$condition)
countData$geneID <- rownames(countData)
smallestGroupSize <- 4
keep <- rowSums(countData >= 10) >= smallestGroupSize
countData <- countData[keep,]

background_genes <- countData %>%
  dplyr::select("geneID") %>%
  unlist() %>%
  as.vector()

# enrichment
deg_files <- list(
  "25v5"  = "../de_genes_25v5.csv",
  "25v10" = "../de_genes_25v10.csv",
  "10v5"  = "../de_genes_10v5.csv"
)

run_one_enrichment <- function(gene_vec, out_name) {
  
  if (length(gene_vec) == 0) {
    message("  No genes for ", out_name)
    return(NULL)
  }
  
  enr <- tryCatch(
    enricher(
      gene_vec,
      TERM2GENE = term2gene,
      TERM2NAME = term2name,
      pvalueCutoff = 0.05,
      universe = background_genes,
      pAdjustMethod = "fdr",
      qvalueCutoff = 0.2
    ),
    error = function(e) {
      message("  Enrichment failed for ", out_name, ": ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(enr)) {
    return(NULL)
  }
  
  enr_df <- tryCatch(
    as.data.frame(enr),
    error = function(e) {
      message("  Could not extract results for ", out_name, ": ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(enr_df) || nrow(enr_df) == 0) {
    message("  No enriched terms for ", out_name)
    return(NULL)
  }
  
  write_csv(enr_df%>%
              filter(!grepl("^GO:", Description)), paste0("enrichment_results.", out_name, ".csv"))
  
  p <- tryCatch(
    dotplot(enr,
            x = "geneRatio",
            color = "p.adjust",
            orderBy = "x",
            showCategory = min(100, nrow(enr_df)),
            font.size = 7),
    error = function(e) {
      message("  Dotplot failed for ", out_name, ": ", e$message)
      return(NULL)
    }
  )
  
  if (!is.null(p)) {
    ggsave(paste0("enrichment_dotplot_", out_name, ".pdf"),
           p, width = 4, height = 9)
  }
  
  return(enr_df)
}

for (nm in names(deg_files)) {
  
  message("Processing: ", nm)
  
  df <- read_csv(deg_files[[nm]], show_col_types = FALSE)
  
  interesting_up <- df %>%
    filter(!is.na(padj), padj < 0.05, log2FoldChange > 0) %>%
    pull(gene_id) %>%
    unique() %>%
    na.omit()
  
  interesting_down <- df %>%
    filter(!is.na(padj), padj < 0.05, log2FoldChange < 0) %>%
    pull(gene_id) %>%
    unique() %>%
    na.omit()
  
  message("  Up genes: ", length(interesting_up))
  message("  Down genes: ", length(interesting_down))
  
  parts <- strsplit(nm, "v")[[1]]
  reverse_nm <- paste0(parts[2], "v", parts[1])
  
  run_one_enrichment(interesting_up, nm)
  run_one_enrichment(interesting_down, reverse_nm)
}



#### bar plot ####

deg_results <- list(
  "30v1"  = read_csv("../de_genes_30v1.csv"),
  "25v5"  = read_csv("../de_genes_25v5.csv"),
  "25v10" = read_csv("../de_genes_25v10.csv"),
  "10v5"  = read_csv("../de_genes_10v5.csv")
)

enrich_files <- c(
  "./enrichment_results.10v25.csv",
  "./enrichment_results.10v5.csv",
  "./enrichment_results.25v10.csv",
  "./enrichment_results.25v5.csv",
  "./enrichment_results.5v10.csv",
  "./enrichment_results.30v1.csv",
  "./enrichment_results.1v30.csv"
)

enrich_results <- map(enrich_files, read_csv, show_col_types = FALSE)

names(enrich_results) <- enrich_files %>%
  basename() %>%
  str_remove("^enrichment_results\\.") %>%
  str_remove("\\.csv$")

compute_mean_logfc <- function(enrich_df, deg_df) {
  enrich_df %>%
    select(ID, Description, geneID, p.adjust) %>%
    separate_rows(geneID, sep = "/") %>%
    left_join(deg_df, by = c("geneID" = "gene_id")) %>%
    group_by(ID, Description, p.adjust) %>%
    summarise(
      mean_logFC = mean(log2FoldChange, na.rm = TRUE),
      .groups = "drop"
    )
}

mean_logfc_named <- imap_dfr(enrich_results, function(enrich_df, nm) {
  
  # nm examples: "10v25", "25v10", "30v1", "1v30"
  
  contrast_key <- case_when(
    nm %in% c("30v1", "1v30")   ~ "30v1",
    nm %in% c("25v5", "5v25")   ~ "25v5",
    nm %in% c("25v10", "10v25") ~ "25v10",
    nm %in% c("10v5", "5v10")   ~ "10v5",
    TRUE ~ NA_character_
  )
  
  if (is.na(contrast_key)) {
    warning("Skipping unexpected contrast: ", nm)
    return(tibble())
  }
  
  deg_df <- deg_results[[contrast_key]]
  if (is.null(deg_df)) {
    warning("No DEG table found for contrast_key = ", contrast_key, " (from file name: ", nm, ")")
    return(tibble())
  }
  
  out <- compute_mean_logfc(enrich_df, deg_df)
  
  out %>% mutate(contrast = contrast_key)
})

top_terms <- mean_logfc_named %>%
  filter(p.adjust < 0.05) %>%
  mutate(
    z_logFC = scale(mean_logFC)[,1],
    Description_ordered = factor(
      Description,
      levels = rev(unique(Description))
    )
  )

contrast_levels <- c("30v1", "25v5", "25v10", "10v5")
top_terms$contrast <- factor(top_terms$contrast, levels = contrast_levels)

contrast_colors <- c(
  "30v1"  = "#1f78b4",
  "25v5"  = "#33a02c",
  "25v10" = "#e31a1c",
  "10v5"  = "#ff7f00"
)

contrast_shapes <- c(
  "30v1"  = 21,
  "25v5"  = 22,
  "25v10" = 23,
  "10v5"  = 24
)

clust.barplot <- ggplot(
  top_terms,
  aes(x = Description_ordered,
      y = mean_logFC,
      fill = contrast,
      shape = contrast)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_col(position = position_dodge(width = 0.8), width = 0.6, alpha = 0.6) +
  geom_point(size = 2, position = position_dodge(width = 0.8), color = "black") +
  coord_flip() +
  scale_y_continuous(
    name = "Mean log2FC",
    breaks = scales::breaks_width(1)
  ) +
  scale_x_discrete(
    labels = function(x) stringr::str_trunc(x, width = 50),
    name = "GO term"
  ) +
  scale_fill_manual(values = contrast_colors) +
  scale_shape_manual(values = contrast_shapes) +
  theme_minimal() +
  guides(fill = guide_legend(reverse = TRUE), shape = guide_legend(reverse = TRUE)) +
  labs(title = "GO term enrichment") +
  theme(
    plot.title = element_text(size = 13, margin = margin(t = 9, b = 6)),
    plot.margin = margin(10, 10, 8, 8),
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 8),
    legend.position = "right"
  )

clust.barplot

ggsave("go_barplot.pdf", clust.barplot, width = 6.5, height = 10)

#### go slim function ####

# 1. Load GO ontology + GO slim

go <- get_ontology("http://purl.obolibrary.org/obo/go.obo",
                   extract_tags = "everything")

goslim <- get_ontology("http://current.geneontology.org/ontology/subsets/goslim_generic.obo",
                       extract_tags = "everything")

# 2. Function to map GO IDs to GO slim

map_go_to_goslim <- function(go_ids, go, goslim) {
  slim_map <- lapply(go_ids, function(go_id) {
    if (is.na(go_id) || !go_id %in% names(go$ancestors)) return(NA_character_)
    
    ancestors <- unique(c(go_id, go$ancestors[[go_id]]))  # include self
    slim_hits <- intersect(ancestors, names(goslim$name))
    
    if (length(slim_hits) == 0) return(NA_character_)
    
    paste(slim_hits, collapse = ";")
  })
  
  unlist(slim_map)
}

# optional: convert GO slim IDs to names too
map_goslim_names <- function(goslim_ids, goslim) {
  sapply(goslim_ids, function(x) {
    if (is.na(x) || x == "") return(NA_character_)
    ids <- strsplit(x, ";")[[1]]
    ids <- ids[ids %in% names(goslim$name)]
    if (length(ids) == 0) return(NA_character_)
    paste(goslim$name[ids], collapse = ";")
  })
}


# 3. Function to process one enrichment file

process_enrichment_file <- function(file,
                                    go,
                                    goslim,
                                    p_cutoff = 0.05,
                                    save_dir = ".") {
  
  message("Processing: ", file)
  
  enrich_results <- read_csv(file, show_col_types = FALSE)
  
  # keep only significant rows if wanted
  enrich_results2 <- enrich_results %>%
    filter(!is.na(p.adjust), p.adjust < p_cutoff)
  
  enrich_results2$GOSLIM_ID <- map_go_to_goslim(enrich_results2$ID, go, goslim)
  enrich_results2$GOSLIM_Name <- map_goslim_names(enrich_results2$GOSLIM_ID, goslim)
  
  # extract contrast name from file name
  contrast_name <- file %>%
    basename() %>%
    str_remove("\\.csv$") %>%
    str_remove("^enrichment_results\\.?") %>%
    str_remove("^lfc[0-9.]+\\.")
  
  enrich_results2 <- enrich_results2 %>%
    mutate(contrast = contrast_name)
  
  if (nrow(enrich_results2) == 0) {
    enrich_results2$GOSLIM_ID <- character(0)
    enrich_results2$GOSLIM_Name <- character(0)
  } else {
    enrich_results2$GOSLIM_ID <- map_go_to_goslim(enrich_results2$ID, go, goslim)
    enrich_results2$GOSLIM_Name <- map_goslim_names(enrich_results2$GOSLIM_ID, goslim)
  }
  
  out_file <- file.path(
    save_dir,
    paste0(tools::file_path_sans_ext(basename(file)), "_goslim.csv")
  )
  
  write_csv(enrich_results2, out_file)
  
  return(enrich_results2)
}

# 4. Define files to process

files_to_process <- c(
  "./enrichment_results.10v25.csv",
  "./enrichment_results.10v5.csv",
  "./enrichment_results.25v10.csv",
  "./enrichment_results.25v5.csv",
  "./enrichment_results.5v10.csv",
  "./enrichment_results.30v1.csv",
  "./enrichment_results.1v30.csv"
)
# if some of these two extra files are already included by list.files(),
# unique() prevents duplication


# 5. Process all files

all_goslim_tables <- lapply(files_to_process, process_enrichment_file,
                            go = go,
                            goslim = goslim,
                            p_cutoff = 0.05,
                            save_dir = ".")

combined_goslim <- bind_rows(all_goslim_tables)

write_csv(combined_goslim, "all_enrichment_results_goslim_combined.csv")



# 6. Prepare combined table for plotting

plot_df <- combined_goslim %>%
  filter(!is.na(GOSLIM_Name), GOSLIM_Name != "") %>%
  separate_rows(GOSLIM_ID, GOSLIM_Name, sep = ";") %>%
  filter(!is.na(GOSLIM_Name), GOSLIM_Name != "")

# summarise per contrast x goslim
heatmap_df <- plot_df %>%
  group_by(contrast, GOSLIM_Name) %>%
  summarise(
    n_terms = n(),
    .groups = "drop"
  )

# optional: keep only the most frequent GO slim terms
top_goslim <- heatmap_df %>%
  group_by(GOSLIM_Name) %>%
  summarise(total_terms = sum(n_terms, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_terms)) %>%
  slice_head(n = 30) %>%
  pull(GOSLIM_Name)

heatmap_df_top <- heatmap_df %>%
  filter(GOSLIM_Name %in% top_goslim)

# order rows by total number of terms
goslim_order <- heatmap_df_top %>%
  group_by(GOSLIM_Name) %>%
  summarise(total_terms = sum(n_terms), .groups = "drop") %>%
  arrange(total_terms) %>%
  pull(GOSLIM_Name)

heatmap_df_top$GOSLIM_Name <- factor(
  heatmap_df_top$GOSLIM_Name,
  levels = goslim_order
)

# heatmap
p_heatmap <- ggplot(heatmap_df_top,
                    aes(x = contrast, y = GOSLIM_Name, fill = n_terms)) +
  geom_tile(color = "white") +
  scale_fill_gradient(name = "Number of GO terms") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  )

ggsave("goslim_heatmap_summary.pdf", p_heatmap, width = 6.5, height = 5)