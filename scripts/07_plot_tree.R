library(ape)
library(ggplot2)
library(ggtree)

tr <- read.tree("results/partA_ml.treefile")
p <- ggtree(tr)

d <- p$data
d$support <- suppressWarnings(as.numeric(d$label))
d$species <- ifelse(d$isTip, sub("^([A-Za-z]+_[a-z]+)_.*$", "\\1", d$label), NA)
d$species <- gsub("_", " ", d$species)
d$lab <- ifelse(d$isTip, gsub("_", " ", d$label), NA)
p$data <- d

p <- p +
  geom_tiplab(aes(label = lab, color = species), size = 3, fontface = "italic") +
  geom_text2(aes(subset = !isTip & !is.na(support) & support >= 70, label = support),
             hjust = 1.15, vjust = -0.5, size = 3) +
  geom_treescale(x = 0, y = 0, width = 0.1) +
  xlim(0, max(d$x) * 2) +
  theme(legend.position = "none")

ggsave("figures/partA_tree.png", p, width = 9, height = 6, dpi = 300, bg = "white")
cat("Saved figures/partA_tree.png\n")
