library(ape)
library(pegas)

x <- read.dna("data/processed/partB_trimmed.fasta", format = "fasta")
meta <- read.csv("data/processed/partB_trimmed_metadata.csv", stringsAsFactors = FALSE)
meta <- meta[match(rownames(x), meta$accession), ]
country <- meta$country

d <- dist.dna(x, model = "raw")
pc <- cmdscale(d, k = 2)
flag <- pc[, 1] > 0.01
cat("Flagged sequences:", sum(flag), "\n")
stopifnot(sum(flag) > 3, sum(flag) < 60)
print(table(country[flag]))

cat("\nFlagged sequences (accession, country, host, title):\n")
fl <- meta[flag, c("accession", "country", "host", "title")]
fl$PC1 <- round(pc[flag, 1], 3)
print(fl[order(fl$country, fl$PC1), ], row.names = FALSE)
write.csv(fl, "results/divergent_sequences.csv", row.names = FALSE)

dn <- as.matrix(dist.dna(x, model = "N"))
cat("\nMean differences (of 587 sites): flagged vs main group:",
    round(mean(dn[flag, !flag]), 1), "\n")
cat("Mean differences within main group:", round(mean(dn[!flag, !flag][upper.tri(dn[!flag, !flag])]), 1), "\n")

x2 <- x[!flag, ]
c2 <- country[!flag]
pops <- c(sort(unique(c2)), "All")
res <- do.call(rbind, lapply(pops, function(p) {
  s <- if (p == "All") x2 else x2[c2 == p, ]
  tj <- tajima.test(s)
  data.frame(population = p, n = nrow(s),
             haplotypes = nrow(haplotype(s)),
             nucleotide_diversity = round(nuc.div(s), 4),
             tajimas_D = round(tj$D, 3),
             D_pvalue = round(tj$Pval.normal, 4))
}))
cat("\nDiversity WITHOUT flagged sequences:\n")
print(res, row.names = FALSE)
write.csv(res, "results/diversity_without_divergent.csv", row.names = FALSE)

prs <- combn(sort(unique(c2)), 2)
fst <- do.call(rbind, lapply(seq_len(ncol(prs)), function(i) {
  idx <- c2 %in% prs[, i]
  dd <- dist.dna(x2[idx, ], model = "N")
  a <- amova(dd ~ factor(c2[idx]), nperm = 999)
  data.frame(pop1 = prs[1, i], pop2 = prs[2, i],
             PhiST = round(a$varcomp[1, 1] / sum(a$varcomp[, 1]), 3),
             p_value = a$varcomp[1, 2])
}))
cat("\nPairwise PhiST WITHOUT flagged sequences:\n")
print(fst, row.names = FALSE)
write.csv(fst, "results/pairwise_phist_without_divergent.csv", row.names = FALSE)
cat("\nDONE\n")
