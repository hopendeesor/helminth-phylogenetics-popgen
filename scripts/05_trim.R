read_aln <- function(f) {
  l <- readLines(f)
  l <- l[nzchar(trimws(l))]
  h <- grep("^>", l)
  ends <- c(h[-1] - 1, length(l))
  nm <- sub(" .*$", "", trimws(sub("^>", "", l[h])))
  seqs <- unlist(mapply(function(s, e) toupper(paste(trimws(l[(s + 1):e]), collapse = "")),
                        h, ends, USE.NAMES = FALSE))
  mat <- do.call(rbind, strsplit(seqs, ""))
  rownames(mat) <- nm
  mat
}
write_aln <- function(mat, f) {
  out <- unlist(lapply(seq_len(nrow(mat)), function(i)
    c(paste0(">", rownames(mat)[i]), paste(mat[i, ], collapse = ""))))
  writeLines(out, f)
}
is_base <- function(m) matrix(m %in% c("A", "C", "G", "T"), nrow = nrow(m))

# ---------------- Part A ----------------
a <- read_aln("data/processed/partA_aligned.fasta")
cat("PART A\n  sequences:", nrow(a), "| alignment length:", ncol(a), "\n")
keepA <- colMeans(is_base(a)) >= 0.8
a_t <- a[, keepA, drop = FALSE]
cat("  columns kept (>=80% of sequences have data):", ncol(a_t), "\n")
cat("  fraction of trimmed alignment covered, per sequence:\n")
print(round(setNames(rowMeans(is_base(a_t)), rownames(a_t)), 2))
write_aln(a_t, "data/processed/partA_trimmed.fasta")

# ---------------- Part B ----------------
b <- read_aln("data/processed/partB_aligned.fasta")
meta <- read.csv("data/processed/partB_3pop_metadata.csv", stringsAsFactors = FALSE)
country <- meta$country[match(rownames(b), meta$accession)]
cat("\nPART B\n  sequences:", nrow(b), "| alignment length:", ncol(b), "\n")
okb <- is_base(b)
core <- which(colMeans(okb) >= 0.9)
cat("  columns with data in >=90% of sequences:", length(core),
    "(span", min(core), "-", max(core), ")\n")
good <- rowSums(!okb[, core, drop = FALSE]) == 0
cat("  sequences complete across those columns:", sum(good), "of", nrow(b), "\n")
cat("  kept per country:\n"); print(table(country[good]))
cat("  dropped per country:\n"); print(table(country[!good]))
b_t <- b[good, core, drop = FALSE]
var_sites <- sum(apply(b_t, 2, function(x) length(unique(x)) > 1))
cat("  final: ", nrow(b_t), "sequences x", ncol(b_t), "sites;", var_sites, "variable sites\n")
write_aln(b_t, "data/processed/partB_trimmed.fasta")
write.csv(meta[meta$accession %in% rownames(b_t), ], "data/processed/partB_trimmed_metadata.csv",
          row.names = FALSE)
cat("\nDONE\n")
