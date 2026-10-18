library(ape)

# ---------- Part B ----------
m <- read.csv("data/processed/Hcontortus_metadata_clean.csv", stringsAsFactors = FALSE)
x <- read.FASTA("data/raw/partB/Hcontortus_COI.fasta")
names(x) <- sub(" .*$", "", names(x))
m <- m[m$accession %in% names(x), ]

m$bin <- cut(m$length, c(0, 400, 600, 700, 2000),
             labels = c("<=400", "401-600", "601-700", ">700"))
big <- names(which(table(m$country[m$country != ""]) >= 10))
cat("Length bins by country (countries with 10+ records):\n")
print(table(m$country[m$country %in% big], m$bin[m$country %in% big]))

keep <- m[m$length >= 600 & m$country != "", ]
cat("\nPart B kept (length >= 600 and country known):", nrow(keep), "\n")
print(sort(table(keep$country), decreasing = TRUE))
write.FASTA(x[keep$accession], "data/processed/partB_filtered.fasta")
write.csv(keep, "data/processed/partB_filtered_metadata.csv", row.names = FALSE)

# ---------- Part A ----------
chosen_all <- list()
cat("\nPart A selected sequences:\n")
for (f in list.files("data/raw/partA", full.names = TRUE)) {
  sp <- sub(".fasta", "", basename(f))
  s <- read.FASTA(f)
  names(s) <- sub(" .*$", "", names(s))
  len <- sapply(s, length)
  ok <- which(len >= 600)
  if (length(ok) == 0) ok <- seq_along(s)
  pick <- head(ok[order(abs(len[ok] - 650))], 3)
  chosen <- s[pick]
  names(chosen) <- paste0(sp, "_", names(chosen))
  chosen_all[[sp]] <- chosen
  for (i in seq_along(chosen)) cat(sprintf("  %-45s %d bp\n", names(chosen)[i], length(chosen[[i]])))
}
write.FASTA(do.call(c, unname(chosen_all)), "data/processed/partA_selected.fasta")
cat("\nDONE\n")
