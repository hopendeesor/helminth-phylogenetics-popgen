library(ape)
library(pegas)

x <- read.dna("data/processed/partB_trimmed.fasta", format = "fasta")
meta <- read.csv("data/processed/partB_trimmed_metadata.csv", stringsAsFactors = FALSE)
country <- meta$country[match(rownames(x), meta$accession)]
stopifnot(!any(is.na(country)))
cat("Sequences:", nrow(x), "| sites:", ncol(x), "\n\n")

pops <- c(sort(unique(country)), "All")
res <- do.call(rbind, lapply(pops, function(p) {
  s <- if (p == "All") x else x[country == p, ]
  tj <- tajima.test(s)
  data.frame(population = p, n = nrow(s),
             segregating_sites = length(seg.sites(s)),
             haplotypes = nrow(haplotype(s)),
             haplotype_diversity = round(hap.div(s), 3),
             nucleotide_diversity = round(nuc.div(s), 4),
             tajimas_D = round(tj$D, 3),
             D_pvalue = round(tj$Pval.normal, 4))
}))
print(res, row.names = FALSE)
write.csv(res, "results/diversity_by_population.csv", row.names = FALSE)

cat("\nPairwise PhiST (AMOVA, 999 permutations):\n")
prs <- combn(sort(unique(country)), 2)
fst <- do.call(rbind, lapply(seq_len(ncol(prs)), function(i) {
  idx <- country %in% prs[, i]
  d <- dist.dna(x[idx, ], model = "N")
  g <- factor(country[idx])
  a <- amova(d ~ g, nperm = 999)
  data.frame(pop1 = prs[1, i], pop2 = prs[2, i],
             PhiST = round(a$varcomp[1, 1] / sum(a$varcomp[, 1]), 3),
             p_value = a$varcomp[1, 2])
}))
print(fst, row.names = FALSE)
write.csv(fst, "results/pairwise_phist.csv", row.names = FALSE)
cat("\nDONE\n")
