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

a <- read_aln("data/processed/partA_aligned.fasta")
rownames(a) <- sub("^_R_", "", rownames(a))
drop <- c("Ascaris_suum_BM277881.1", "Ascaris_suum_BG733781.1")
a <- a[!rownames(a) %in% drop, , drop = FALSE]
a <- a[, colSums(a != "-") > 0, drop = FALSE]
cat("Sequences:", nrow(a), "| alignment length after removing empty columns:", ncol(a), "\n")

keep <- colMeans(is_base(a)) >= 0.8
a_t <- a[, keep, drop = FALSE]
cat("Columns kept (>=80% of sequences have data):", ncol(a_t), "\n")
cat("Coverage per sequence:\n")
print(round(setNames(rowMeans(is_base(a_t)), rownames(a_t)), 2))
write_aln(a_t, "data/processed/partA_trimmed.fasta")
cat("DONE\n")
