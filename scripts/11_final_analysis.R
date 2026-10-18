library(ape)
library(pegas)
library(ggplot2)
library(patchwork)

x <- read.dna("data/processed/partB_trimmed.fasta", format = "fasta")
meta <- read.csv("data/processed/partB_trimmed_metadata.csv", stringsAsFactors = FALSE)
meta <- meta[match(rownames(x), meta$accession), ]
country <- meta$country

pc0 <- cmdscale(dist.dna(x, model = "raw"), k = 2)
excl <- country == "Nigeria" & pc0[, 1] > 0.05
cat("Excluded (H. placei-like):", sum(excl), "\n")
stopifnot(sum(excl) == 15)
write.csv(meta[excl, c("accession", "country", "host", "title")],
          "results/excluded_placei_like.csv", row.names = FALSE)

x2 <- x[!excl, ]
m2 <- meta[!excl, ]
c2 <- m2$country
cat("Final dataset:", nrow(x2), "sequences\n")
print(table(c2))

chars <- as.character(x2)
writeLines(unlist(lapply(seq_len(nrow(x2)), function(i)
  c(paste0(">", rownames(x2)[i]), paste(chars[i, ], collapse = "")))),
  "data/processed/partB_final.fasta")
write.csv(m2, "data/processed/partB_final_metadata.csv", row.names = FALSE)

pops <- c(sort(unique(c2)), "All")
res <- do.call(rbind, lapply(pops, function(p) {
  s <- if (p == "All") x2 else x2[c2 == p, ]
  tj <- tajima.test(s)
  data.frame(population = p, n = nrow(s),
             segregating_sites = length(seg.sites(s)),
             haplotypes = nrow(haplotype(s)),
             haplotype_diversity = round(hap.div(s), 3),
             nucleotide_diversity = round(nuc.div(s), 4),
             tajimas_D = round(tj$D, 3),
             D_pvalue = round(tj$Pval.normal, 4))
}))
cat("\nFINAL diversity:\n")
print(res, row.names = FALSE)
write.csv(res, "results/final_diversity.csv", row.names = FALSE)

prs <- combn(sort(unique(c2)), 2)
fst <- do.call(rbind, lapply(seq_len(ncol(prs)), function(i) {
  idx <- c2 %in% prs[, i]
  dd <- dist.dna(x2[idx, ], model = "N")
  g <- factor(c2[idx])
  a <- amova(dd ~ g, nperm = 999)
  data.frame(pop1 = prs[1, i], pop2 = prs[2, i],
             PhiST = round(a$varcomp[1, 1] / sum(a$varcomp[, 1]), 3),
             p_value = a$varcomp[1, 2])
}))
cat("\nFINAL pairwise PhiST:\n")
print(fst, row.names = FALSE)
write.csv(fst, "results/final_pairwise_phist.csv", row.names = FALSE)

pc <- cmdscale(dist.dna(x2, model = "raw"), k = 2, eig = TRUE)
ev <- pc$eig[pc$eig > 0]
pct <- round(100 * pc$eig[1:2] / sum(ev), 1)
df <- data.frame(PC1 = pc$points[, 1], PC2 = pc$points[, 2], country = c2)
p1 <- ggplot(df, aes(PC1, PC2, color = country)) +
  geom_point(alpha = 0.6, size = 1.8) +
  labs(x = paste0("PC1 (", pct[1], "%)"), y = paste0("PC2 (", pct[2], "%)"),
       color = "Country", title = "PCoA of H. contortus COI haplotypes") +
  theme_classic()

dn <- as.matrix(dist.dna(x2, model = "N"))
pw <- do.call(rbind, lapply(sort(unique(c2)), function(p) {
  m <- dn[c2 == p, c2 == p]
  data.frame(population = p, diff = m[upper.tri(m)])
}))
p2 <- ggplot(pw, aes(diff, fill = population)) +
  geom_density(alpha = 0.4) +
  labs(x = "Pairwise differences (out of 587 sites)", y = "Density",
       fill = "Country", title = "Differences within each population") +
  theme_classic()
ggsave("figures/partB_final.png", p1 + p2, width = 13, height = 5, dpi = 300, bg = "white")
cat("\nSaved figures/partB_final.png\nDONE\n")
