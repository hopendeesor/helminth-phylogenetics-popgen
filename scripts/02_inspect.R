library(ape)

cat("PART A\n")
for (f in list.files("data/raw/partA", full.names = TRUE)) {
  x <- read.FASTA(f)
  len <- sapply(x, length)
  cat(sprintf("  %-32s n=%.0f  length %.0f-%.0f (median %.0f)\n",
              sub(".fasta", "", basename(f)), length(x), min(len), max(len), median(len)))
}

m <- read.csv("data/raw/partB/Hcontortus_metadata.csv", stringsAsFactors = FALSE)

get_field <- function(subtype, subname, field) {
  k <- strsplit(subtype, "|", fixed = TRUE)[[1]]
  v <- strsplit(subname, "|", fixed = TRUE)[[1]]
  i <- match(field, k)
  if (is.na(i) || i > length(v)) "" else v[i]
}
m$country_raw <- mapply(get_field, m$subtype, m$subname, "country", USE.NAMES = FALSE)
m$host        <- mapply(get_field, m$subtype, m$subname, "host", USE.NAMES = FALSE)
m$country     <- trimws(sub(":.*$", "", m$country_raw))

cat("\nPART B\n")
cat("  records:", nrow(m), "\n")
cat("  length range:", min(m$length), "-", max(m$length), "\n")
cat("  with country:", sum(m$country != ""), "\n")
cat("  with host:", sum(m$host != ""), "\n")
cat("\nRecords per country:\n")
print(sort(table(m$country[m$country != ""]), decreasing = TRUE))
cat("\nTop hosts:\n")
print(head(sort(table(m$host[m$host != ""]), decreasing = TRUE), 10))
cat("\nLength distribution:\n")
print(table(cut(m$length, c(0, 400, 600, 700, 1000, 2000))))

write.csv(m[, c("accession", "length", "country", "host", "title")],
          "data/processed/Hcontortus_metadata_clean.csv", row.names = FALSE)
