library(ape)
m <- read.csv("data/processed/partB_filtered_metadata.csv", stringsAsFactors = FALSE)
pops <- c("Nigeria", "Bangladesh", "Bosnia and Herzegovina")
m <- m[m$country %in% pops, ]

x <- read.FASTA("data/processed/partB_filtered.fasta")
names(x) <- trimws(sub(" .*$", "", names(x)))
cat("Sequences in FASTA:", length(x), "\n")
cat("Rows in metadata:  ", nrow(m), "\n")

idx <- match(m$accession, names(x))
cat("Unmatched accessions:", sum(is.na(idx)), "\n")
m <- m[!is.na(idx), ]
idx <- idx[!is.na(idx)]
x <- x[idx]

chars <- as.character(x)
out <- unlist(lapply(seq_along(chars), function(i)
  c(paste0(">", names(x)[i]), paste(chars[[i]], collapse = ""))))
writeLines(out, "data/processed/partB_3pop.fasta")
write.csv(m, "data/processed/partB_3pop_metadata.csv", row.names = FALSE)
cat("Sequences kept:", nrow(m), "\n")
print(table(m$country))
