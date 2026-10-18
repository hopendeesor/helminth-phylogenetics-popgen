library(rentrez)

dir.create("data/raw/partA", recursive = TRUE, showWarnings = FALSE)
dir.create("data/raw/partB", recursive = TRUE, showWarnings = FALSE)

species <- c("Haemonchus contortus", "Teladorsagia circumcincta",
             "Trichostrongylus colubriformis", "Ostertagia ostertagi",
             "Cooperia oncophora", "Oesophagostomum dentatum",
             "Ascaris suum")

coi <- paste0('(COI[Gene] OR cox1[Gene] OR "cytochrome c oxidase subunit 1"[Title]',
              ' OR "cytochrome c oxidase subunit I"[Title])')

g <- function(v) if (is.null(v)) "" else paste(v, collapse = "")

# ---- Part A: COI for each species ----
cat("PART A\n")
for (sp in species) {
  term <- sprintf('"%s"[Organism] AND %s AND 500:1700[SLEN]', sp, coi)
  s <- entrez_search(db = "nuccore", term = term, retmax = 50)
  cat(sprintf("  %-35s %d records\n", sp, s$count))
  Sys.sleep(0.5)
  if (length(s$ids) > 0) {
    fa <- entrez_fetch(db = "nuccore", id = s$ids, rettype = "fasta", retmode = "text")
    writeLines(fa, file.path("data/raw/partA", paste0(gsub(" ", "_", sp), ".fasta")))
    Sys.sleep(0.5)
  }
}

# ---- Part B: all H. contortus COI + metadata ----
cat("PART B\n")
termB <- sprintf('"Haemonchus contortus"[Organism] AND %s AND 300:1700[SLEN]', coi)
sB <- entrez_search(db = "nuccore", term = termB, retmax = 2000)
cat("  H. contortus COI records:", sB$count, "\n")

fasta_all <- character(0)
rows_all <- list()
for (ch in split(sB$ids, ceiling(seq_along(sB$ids) / 100))) {
  fasta_all <- c(fasta_all, entrez_fetch(db = "nuccore", id = ch, rettype = "fasta", retmode = "text"))
  Sys.sleep(0.5)
  sm <- entrez_summary(db = "nuccore", id = ch)
  if (inherits(sm, "esummary")) sm <- list(sm)
  rows_all <- c(rows_all, lapply(sm, function(x)
    data.frame(accession = g(x$accessionversion), length = g(x$slen),
               title = g(x$title), subtype = g(x$subtype), subname = g(x$subname),
               stringsAsFactors = FALSE)))
  Sys.sleep(0.5)
}
writeLines(fasta_all, "data/raw/partB/Hcontortus_COI.fasta")
write.csv(do.call(rbind, rows_all), "data/raw/partB/Hcontortus_metadata.csv", row.names = FALSE)
cat("DONE\n")
