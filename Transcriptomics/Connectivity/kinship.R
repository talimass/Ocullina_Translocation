library(igraph)
library(ggraph)
library(ggplot2)
library(stringr)



#save.image(file='kinship5.RData')
#load('kinship.RData')
#### all clones/no relatives ####
setwd("/home/gospozha/haifa/hiba/op_align_new/snp/pruned")
# Read KING table
king <- read.csv("idb.kin0", stringsAsFactors=FALSE, sep = "\t")
#king <- read.csv("idb_0.7.kin0", stringsAsFactors=FALSE, sep = "\t") # all
#king <- read.csv("./new/idb_w_rel.kin0", stringsAsFactors=FALSE, sep = "\t") # no_relatives
#king <- read.csv("idb.kin0.no_known_clones.tsv", stringsAsFactors=FALSE, sep = "\t") # no_relatives


# Define thresholds
clone_cutoff <- 0.354      # definite clones
firstdeg_cutoff <- 0.177   # first-degree relatives

# Create a new column for relationship type
king$rel_type <- ifelse(king$KINSHIP > clone_cutoff, "Clone",
                        ifelse(king$KINSHIP > firstdeg_cutoff, "1st-degree", NA))

# Keep only edges above 0.177
edges <- subset(king, !is.na(rel_type), select=c(ID1, ID2, KINSHIP, rel_type))

# Build graph
g <- graph_from_data_frame(edges, directed=FALSE)
# Short label: remove everything until first '-' (including '-')
V(g)$label <- sub("-[^-]*$", "", V(g)$name)

# Site from full name (keeps working even after shortening)
V(g)$site <- str_extract(V(g)$name, "(CC|MF)")

# Depth class from full name
V(g)$depth <- dplyr::case_when(
  grepl("(40|DS|DD)", V(g)$name) ~ "Deep",
  grepl("(10|SS|SD)", V(g)$name) ~ "Shallow",
  TRUE ~ NA_character_
)

# Plot
kinship <- ggraph(g, layout="fr") +
  geom_edge_link(aes(color=rel_type, width=KINSHIP), alpha=0.6) +
  
  # nodes: shape = site, color = depth
  geom_node_point(aes(shape=site, color=depth), size=4) +
  
  # labels: shortened
  geom_node_text(aes(label=label), repel=TRUE, size=3) +
  
  #ggtitle("Pairwise kinship relationships among samples") +
  
  # edge colors (pick anything you like here)
  scale_edge_color_manual(values=c("Clone"="#6F2DBD", "1st-degree"="#F4A261")) +
  
  # node colors: your requested palette
  scale_color_manual(values=c("Shallow"="#f75f55", "Deep"="#00A9FF"), na.value="grey70") +
  
  # node shapes by site
  scale_shape_manual(values=c("CC"=16, "MF"=17)) +
  
  scale_edge_width_continuous(name="Kinship", range=c(0.3, 2)) +
  theme_void() +
  labs(edge_color="Relationship", color="Origin", shape="Site")

kinship
#saveRDS(kinship, "fig_kinship.rds")
ggsave("kinship.jpg", kinship, width = 10, height = 10)

#### removing artificial clones ####

library(dplyr)
library(readr)
library(purrr)

# paths
meta_path <- "../Metadata.csv"
king_path <- "idb_0.7.kin0"

# 1) read metadata and clean sample IDs to match KING (metadata has "./" prefix)
meta <- read_csv(meta_path, show_col_types = FALSE) %>%
  mutate(id_clean = sub("^\\./", "", id)) %>%
  select(id_clean, origin)

# 2) build a set of "clone-pair keys" from metadata (all within-group pairs)
clone_pairs <- meta %>%
  filter(!is.na(origin)) %>%
  group_by(origin) %>%
  summarise(ids = list(unique(id_clean)), .groups = "drop") %>%
  mutate(pair_keys = map(ids, \(v) {
    v <- sort(v)
    if (length(v) < 2) return(character(0))
    cmb <- combn(v, 2)
    paste(cmb[1, ], cmb[2, ], sep = "__")
  })) %>%
  pull(pair_keys) %>%
  unlist(use.names = FALSE) %>%
  unique()

# 3) read KING table (first header is "#FID1", keep it)
king <- read_delim(king_path, delim = "\t", show_col_types = FALSE, comment = "") %>%
  rename_with(~ sub("^#", "", .x))  # turns "#FID1" -> "FID1"

# 4) make an order-invariant key for each pair in KING and remove clone pairs
king_filtered <- king %>%
  mutate(
    a = pmin(ID1, ID2),
    b = pmax(ID1, ID2),
    pair_key = paste(a, b, sep = "__")
  ) %>%
  filter(!pair_key %in% clone_pairs) %>%
  select(-a, -b, -pair_key)

# 5) optional: write out
write_tsv(king_filtered, "idb_0.7.kin0.no_known_clones.tsv")

# quick sanity check
cat("Original rows:", nrow(king), "\n")
cat("Filtered rows:", nrow(king_filtered), "\n")
cat("Removed rows:", nrow(king) - nrow(king_filtered), "\n")



