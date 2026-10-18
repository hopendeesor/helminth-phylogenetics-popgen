library(ape)
library(ggplot2)
library(patchwork)

x <- read.dna("data/processed/partB_trimmed.fasta", format = "fasta")
meta <- read.csv("data/processed/partB_trimmed_metadata.csv", stringsAsFactors = FALSE)
country <- meta$country[match(rownames(x), meta$accession)]

# PCoA of pairwise genetic distances
d <- dist.dna(x, model = "raw")
pc <- cmdscale(d, k = 2, eig = TRUE)
ev <- pc$eig[pc$eig > 0]
pct <- round(100 * pc$eig[1:2] / sum(ev), 1)
df <- data.frame(PC1 = pc$points[, 1], PC2 = pc$points[, 2], country = country)
p1 <- ggplot(df, aes(PC1, PC2, color = country)) +
  geom_point(alpha = 0.6, size = 1.8) +
  labs(x = paste0("PC1 (", pct[1], "%)"), y = paste0("PC2 (", pct[2], "%)"),
       color = "Country", title = "PCoA of COI haplotypes") +
  theme_classic()

# Pairwise differences within each population
dn <- as.matrix(dist.dna(x, model = "N"))
pw <- do.call(rbind, lapply(sort(unique(country)), function(p) {
  m <- dn[country == p, country == p]
  data.frame(population = p, diff = m[upper.tri(m)])
}))
p2 <- ggplot(pw, aes(diff, fill = population)) +
  geom_density(alpha = 0.4) +
  labs(x = "Pairwise differences (out of 587 sites)", y = "Density",
       fill = "Country", title = "Differences within each population") +
  theme_classic()

ggsave("figures/partB_diversity.png", p1 + p2, width = 13, height = 5, dpi = 300, bg = "white")
cat("Saved figures/partB_diversity.png\n")
