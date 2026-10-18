# Phylogeny and population genetics of livestock nematodes from public COI data

A small, fully scripted analysis of mitochondrial cytochrome c oxidase subunit 1 (COI) sequences from public NCBI records, in two parts:

- **Part A:** a maximum-likelihood phylogeny of six livestock strongyle nematodes, rooted on *Ascaris suum*.
- **Part B:** genetic diversity and population differentiation of *Haemonchus contortus* (barber's pole worm, sheep and goats) in Nigeria, Bangladesh and Bosnia and Herzegovina.

**Question:** How genetically diverse is *H. contortus* across countries, and how does it sit among other livestock nematodes?

## Data
All sequences are public GenBank records (COI/cox1), downloaded with the R package rentrez (scripts/01_download.R).

- Part A: 17 sequences from 7 species (*H. contortus*, *Teladorsagia circumcincta*, *Trichostrongylus colubriformis*, *Ostertagia ostertagi*, *Cooperia oncophora*, *Oesophagostomum dentatum*, and the outgroup *Ascaris suum*).
- Part B: 659 *H. contortus* COI records downloaded. Kept sequences of 600 bp or more with a known country, and countries with at least 10 such sequences: 301 sequences from Nigeria (182), Bangladesh (77) and Bosnia and Herzegovina (42), over 587 aligned sites.

## Methods
- Alignment: MAFFT v7 web server (auto strategy, direction adjusted to the first sequence; 200PAM for Part A, 1PAM for Part B). The alignment steps were done on the web server and are not scripted; the aligned files are in data/processed.
- Trimming: Part A kept columns with data in at least 80% of sequences (415 bp). Part B kept columns with data in at least 90% of sequences, then removed sequences with gaps or ambiguous bases in that region (scripts/05_trim.R).
- Tree: IQ-TREE 3.1.4, model chosen by ModelFinder (HKY+F+G4 by BIC), 1000 ultrafast bootstraps, rooted on *A. suum*.
- Diversity (R packages ape and pegas): haplotypes, haplotype diversity, nucleotide diversity (pi), Tajima's D, and pairwise PhiST from an AMOVA on pairwise differences with 999 permutations.
- Structure: PCoA of pairwise distances.

## Results
**Part A.** Each species with more than one sequence formed a well-supported group (bootstrap 90-100). Relationships among genera were weakly supported (bootstrap 24-54), as expected for a short fragment of one mitochondrial gene. Figure: figures/partA_tree.png.

**Part B.** Figure: figures/partB_final.png. Tables: results/final_diversity.csv and results/final_pairwise_phist.csv.

| Population | n | Haplotypes | pi | Tajima's D (p) |
|---|---|---|---|---|
| Bangladesh | 77 | 77 | 0.0292 | -1.02 (0.31) |
| Bosnia and Herzegovina | 42 | 42 | 0.0346 | -1.16 (0.24) |
| Nigeria | 182 | 178 | 0.0321 | -1.40 (0.16) |
| All | 301 | 297 | 0.0399 | -1.11 (0.27) |

Haplotype diversity was about 1 in every population, so it does not distinguish them. Tajima's D was negative but not significant anywhere. Pairwise PhiST (all p < 0.001): Bangladesh vs Bosnia and Herzegovina 0.22, Bangladesh vs Nigeria 0.51, Bosnia and Herzegovina vs Nigeria 0.51.

**Divergent lineage excluded.** In the PCoA, 15 Nigerian sequences labelled *H. contortus* formed a separate cluster, about 65 differences out of 587 sites from the rest (about 23 within the main group). A BLAST search of one representative (OL435511.1) matched *H. placei* mitochondrial genomes at about 98% and *H. contortus* reference mitochondrial genomes at about 90%; the other 14 appeared among its top hits at 96-99%. These sequences are treated as *H. placei*-like and excluded from the main analysis (list in results/excluded_placei_like.csv). With them included (316 sequences), Nigeria's pi was 0.0446 and PhiST between Nigeria and the other countries was about 0.32.

## Limitations
- One short mitochondrial marker (587 bp); no genome-scale data.
- Country is confounded with study: each country's sequences come from a few submitting groups, so PhiST may partly reflect study, host and sampling differences.
- Unequal sample sizes (182, 77, 42).
- The exclusion of the divergent Nigerian lineage rests on BLAST of one representative plus clustering, not a formal species identification; hybridisation between *H. contortus* and *H. placei* can complicate mitochondrial markers.
- Four Bosnian sequences are also divergent (not *H. placei*-like) and were kept; they raise Bosnia's pi (0.0346 vs 0.0268 without them).
- Brazilian sequences (127) were excluded because they are short fragments covering a different region.
- Part A has one to three sequences per species.

## Reproducing
Requires R with ape, pegas, ggplot2, patchwork, ggtree, rentrez and adegenet, and IQ-TREE 3.1.4. Run scripts 01 to 11 in order from the project root, for example:

    Rscript scripts/01_download.R

The MAFFT alignments (between scripts 04 and 05) are made on the web server; the aligned files are included in data/processed.

## Repository layout
- data/raw: sequences as downloaded
- data/processed: filtered, aligned and trimmed data
- scripts: R scripts, numbered in order of use
- results: tree files, statistics tables
- figures: tree and diversity figures
