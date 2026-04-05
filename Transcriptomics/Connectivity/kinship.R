library(igraph)
library(ggraph)
library(ggplot2)
library(stringr)


setwd("/home/gospozha/haifa/hiba/op_align_new/snp/pruned")
# Read KING table
king <- read.csv("idb.kin0", stringsAsFactors=FALSE, sep = "\t")


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


