cat > README.md << 'EOF'
# Mitochondrial COI diversity of *Haemonchus contortus* and phylogeny of livestock strongyle nematodes from public sequence data

**Author:** Hope Ndeesor
**Data retrieved:** 18 October 2026

## Abstract

Gastrointestinal nematodes are a major constraint on small-ruminant production, and *Haemonchus contortus* is among the most economically important. Using only public records from GenBank, this project (i) reconstructed a maximum-likelihood phylogeny of six livestock strongyle nematodes from a 415 bp fragment of the mitochondrial cytochrome c oxidase subunit 1 gene (COI), and (ii) quantified haplotype and nucleotide diversity, neutrality (Tajima's D) and between-country differentiation (PhiST) in 301 *H. contortus* COI sequences from Nigeria, Bangladesh and Bosnia and Herzegovina. Species with multiple sequences were recovered as well-supported clades (ultrafast bootstrap 90-100), whereas relationships among genera were poorly resolved (24-54). Ordination of pairwise distances revealed a cluster of 15 Nigerian sequences labelled *H. contortus* that were about 12% divergent from the remainder and matched *H. placei* mitochondrial genomes at about 98% in a BLAST search; these were excluded from the population analysis. Nucleotide diversity was similar among countries (pi = 0.029-0.035), Tajima's D was negative but not significant in any population, and PhiST indicated substantial differentiation (0.22-0.51; p < 0.001). Because country is confounded with study of origin, these differences should be interpreted cautiously.

## 1. Background and objectives

*Haemonchus contortus* is a blood-feeding abomasal nematode of sheep and goats with a worldwide distribution, and anthelmintic resistance in this species is widespread. Mitochondrial markers such as COI are widely deposited in public databases and are commonly used for species identification and for describing population structure. Public records are, however, heterogeneous in length, coverage and annotation accuracy, so their reuse requires explicit quality control.

The objectives were to:

1. infer the evolutionary relationships among six livestock strongyle nematodes from COI;
2. describe the genetic diversity of *H. contortus* within and among countries using publicly available COI sequences; and
3. assess the quality of the public data, including the detection of mislabelled or divergent records.

## 2. Data

All sequences are GenBank records retrieved with the R package rentrez (`scripts/01_download.R`) using the queries `"<species>"[Organism]` AND (COI[Gene] OR cox1[Gene] OR "cytochrome c oxidase subunit 1"[Title] OR "cytochrome c oxidase subunit I"[Title]) with length limits of 500-1700 bp (Part A) and 300-1700 bp (Part B).

**Part A (phylogeny).** Species: *H. contortus*, *Teladorsagia circumcincta*, *Trichostrongylus colubriformis*, *Ostertagia ostertagi*, *Cooperia oncophora*, *Oesophagostomum dentatum*, with *Ascaris suum* as outgroup. Up to three sequences of at least 600 bp per species, those closest to 650 bp, were used (one for *O. dentatum*, for which one record was available).

**Part B (population genetics).** All 659 *H. contortus* COI records available at the time of retrieval were downloaded. Country was taken from the GenBank `country` qualifier.

| Step | Sequences |
|---|---|
| Records retrieved | 659 |
| With country annotation | 552 |
| Length of at least 600 bp and country known | 353 |
| Countries with at least 10 such sequences (Nigeria, Bangladesh, Bosnia and Herzegovina) | 341 |
| Complete (no gaps or ambiguity) across the trimmed alignment (587 sites) | 316 |
| After exclusion of a divergent *H. placei*-like lineage (Section 4.3) | 301 |

All 127 Brazilian records were shorter than 600 bp (122 of them at most 400 bp) and Turkey fell to three sequences after length filtering, so neither was analysed. Final sample sizes: Nigeria 182, Bangladesh 77, Bosnia and Herzegovina 42. Because GenBank is updated continuously, a new download will not reproduce these counts exactly; the filtered datasets and alignments used here are committed in `data/processed`.

## 3. Methods

**Alignment.** Sequences were aligned with the MAFFT v7 web server (automatic strategy; direction adjusted to the first sequence; nucleotide scoring matrix 200PAM for Part A and 1PAM for Part B). This step was performed manually on the web server and is not scripted; the resulting alignments are provided.

**Alignment trimming.** For Part A, columns with an unambiguous base in at least 80% of sequences were retained (415 bp). Two of the three *A. suum* sequences (BM277881.1, BG733781.1) covered only about 30% of the trimmed region and were removed before trimming. For Part B, columns with a base in at least 90% of sequences were retained (587 sites), and sequences with a gap or ambiguity code within these columns were discarded (`scripts/05_trim.R`, `scripts/06_trim_partA.R`).

**Phylogenetic inference.** A maximum-likelihood tree was inferred with IQ-TREE 3.1.4. The substitution model was selected with ModelFinder (HKY+F+G4 by BIC), branch support was estimated with 1000 ultrafast bootstrap replicates, and the tree was rooted on *A. suum*.

**Population genetic statistics.** Using the R packages ape and pegas, the following were calculated per country and for all sequences pooled: number of haplotypes, haplotype diversity, nucleotide diversity (pi), the number of segregating sites and Tajima's D (p-value from the normal approximation). Pairwise PhiST between countries was estimated by AMOVA on the matrix of pairwise nucleotide differences with 999 permutations. Genetic structure was explored by principal coordinates analysis (PCoA) of pairwise proportions of differing sites.

**Detection of divergent sequences.** Sequences with PC1 greater than 0.01 in the PCoA were flagged. One representative of the divergent cluster (OL435511.1) was compared against the NCBI core nucleotide database with BLASTN (the "somewhat similar sequences" option) on 18 October 2026.

## 4. Results

### 4.1 Phylogeny (Part A)

The alignment comprised 17 sequences and 415 columns, of which 108 were parsimony-informative. Each of the five species represented by more than one sequence formed a clade with high support (*H. contortus* 100, *T. colubriformis* 100, *T. circumcincta* 98, *O. ostertagi* 96, *C. oncophora* 90). Support for relationships among genera was low (24-54). The tree is shown in `figures/partA_tree.png`; the tree file and the IQ-TREE report are `results/partA_ml.treefile` and `results/partA_ml.iqtree`. The outgroup *A. suum* failed the base-composition test (p = 2.8%), which is expected for a distant outgroup.

### 4.2 Population genetics of *H. contortus* (Part B)

Final dataset: 301 sequences, 587 sites, 230 segregating sites, 297 haplotypes. Haplotype diversity was approximately 1.0 in every population and therefore does not discriminate among them.

| Population | n | Segregating sites | Haplotypes | pi | Tajima's D | p |
|---|---|---|---|---|---|---|
| Bangladesh | 77 | 120 | 77 | 0.0292 | -1.016 | 0.31 |
| Bosnia and Herzegovina | 42 | 128 | 42 | 0.0346 | -1.164 | 0.24 |
| Nigeria | 182 | 194 | 178 | 0.0321 | -1.398 | 0.16 |
| All | 301 | 230 | 297 | 0.0399 | -1.113 | 0.27 |

Pairwise PhiST (all p < 0.001, 999 permutations):

| Pair | PhiST |
|---|---|
| Bangladesh vs Bosnia and Herzegovina | 0.224 |
| Bangladesh vs Nigeria | 0.507 |
| Bosnia and Herzegovina vs Nigeria | 0.506 |

The PCoA and the distributions of pairwise differences within each country are shown in `figures/partB_final.png`. Tables: `results/final_diversity.csv`, `results/final_pairwise_phist.csv`.

### 4.3 A divergent lineage labelled *H. contortus*

In the PCoA of the initial 316-sequence dataset, 15 Nigerian sequences (all from goats, from several Nigerian states) formed a cluster clearly separated from all others along PC1. They differed from the remaining sequences at about 70 of 587 sites (about 12%), compared with about 23 sites on average among the remaining sequences. In BLASTN, the representative sequence matched *H. placei* mitochondrial genomes (AP017687.1, NC_029736.1) at 98.3-98.6% identity, but *H. contortus* reference mitochondrial genomes (LS997568.1, NC_010383.2) at only about 90%; the other 14 sequences of the cluster were among its top hits at 96-99%. These sequences are therefore regarded as *H. placei*-like and were excluded from the main analysis (`results/excluded_placei_like.csv`).

Effect of the exclusion (sensitivity analysis, all 316 sequences): pi in Nigeria was 0.0446 (0.0321 after exclusion) and PhiST between Nigeria and the other countries was about 0.32 (0.51 after exclusion). The earlier results are in `results/diversity_by_population.csv` and `results/pairwise_phist.csv`.

Four Bosnian sequences were also displaced along PC1, but not to the extent of the *H. placei*-like cluster; they were retained. Without them pi in Bosnia and Herzegovina would be 0.0268 (`results/diversity_without_divergent.csv`).

## 5. Discussion

COI resolved species boundaries among the strongyles examined but not the order of divergence among genera. This is consistent with limited phylogenetic signal in a 415 bp fragment of a single, fast-evolving mitochondrial gene, and no inference about relationships among genera is drawn.

The divergent Nigerian cluster illustrates that species labels in public databases cannot be assumed correct. Hybridisation between *H. contortus* and *H. placei* has been reported, and a hybrid record appeared among the BLAST hits; mitochondrial introgression could therefore also produce this pattern, so a single marker cannot establish the identity of these sequences. Resolving it would require nuclear markers or genome-scale data.

After exclusion, nucleotide diversity was of a similar magnitude in all three countries, and negative but non-significant values of Tajima's D provide no evidence of recent population expansion or selection. The PhiST values indicate that a considerable proportion of the variance lies among countries, most strongly between Nigeria and the other two. However, the sequences from each country originate from a small number of submissions, so country is confounded with study, host and sampling scheme, and the PhiST values cannot be interpreted as purely geographic structure.

## 6. Limitations

- A single short mitochondrial marker (587 bp); no nuclear or genome-scale data.
- Country is confounded with study of origin; sample sizes are unequal (182, 77, 42).
- Exclusion of the divergent lineage is based on ordination and BLAST of one representative, not on formal species identification.
- Four divergent Bosnian sequences were retained and raise nucleotide diversity there.
- Part A has one to three sequences per species, so the phylogeny tests species monophyly but not intraspecific structure.
- Brazilian sequences (127) were not analysed because they cover a different, shorter region of the gene.
- The alignment step uses a web server and cannot be rerun from the command line as supplied.

## 7. Reproducing the analysis

### 7.1 Requirements

- R (developed with 4.6.1), with the packages rentrez, ape, pegas, ggplot2, patchwork and ggtree
- IQ-TREE 3.1.4
- A Unix-style shell (developed in Git Bash on Windows; macOS and Linux also work) and Git
- Internet access for the download step and for MAFFT

### 7.2 Setup

    git clone https://github.com/hopendeesor/helminth-phylogenetics-popgen
    cd helminth-phylogenetics-popgen

Install the R packages:

    Rscript -e "install.packages(c('rentrez','ape','pegas','ggplot2','patchwork','BiocManager'), repos='https://cloud.r-project.org')"
    Rscript -e "BiocManager::install('ggtree', ask=FALSE, update=FALSE)"

Install IQ-TREE by downloading the release for your system from https://iqtree.github.io/ and unpacking it. Below, `iqtree3` stands for the path to the executable (for example `iqtree/iqtree-3.1.4-Windows/bin/iqtree3.exe`). All commands must be run from the repository root.

### 7.3 Quick start (from the committed alignments)

This reproduces all results without downloading data again or using the MAFFT web server.

    Rscript scripts/05_trim.R
    Rscript scripts/06_trim_partA.R
    iqtree3 -s data/processed/partA_trimmed.fasta -m MFP -B 1000 -bnni -o Ascaris_suum_MK282844.1 --prefix results/partA_ml
    Rscript scripts/07_plot_tree.R
    Rscript scripts/08_diversity.R
    Rscript scripts/09_figures.R
    Rscript scripts/10_divergent.R
    Rscript scripts/11_final_analysis.R

The IQ-TREE run takes about one minute and the AMOVA permutations a few minutes on a four-core laptop.

### 7.4 Full pipeline (from NCBI)

1. Run `Rscript scripts/01_download.R`, then `02_inspect.R`, `03_prepare.R` and `04_subset.R`.
2. Align `data/processed/partA_selected.fasta` and `data/processed/partB_3pop.fasta` on the MAFFT web server (https://mafft.cbrc.jp/alignment/server/) with the settings in Section 3, and save the FASTA output as `data/processed/partA_aligned.fasta` and `data/processed/partB_aligned.fasta`.
3. Continue with the quick start, from `05_trim.R`.

Because GenBank changes over time, a fresh download may differ from the dataset analysed here.

### 7.5 Scripts and outputs

| Script | Purpose | Main outputs |
|---|---|---|
| `01_download.R` | Retrieve COI records from NCBI | `data/raw/` |
| `02_inspect.R` | Summarise lengths, countries and hosts | `data/processed/Hcontortus_metadata_clean.csv` |
| `03_prepare.R` | Apply length and country filters; select Part A sequences | `partA_selected.fasta`, `partB_filtered.*` |
| `04_subset.R` | Restrict to the three countries | `partB_3pop.*` |
| MAFFT (manual) | Multiple sequence alignment | `partA_aligned.fasta`, `partB_aligned.fasta` |
| `05_trim.R` | Trim alignments and remove incomplete sequences | `partB_trimmed.*` |
| `06_trim_partA.R` | Remove low-coverage outgroup sequences and re-trim Part A | `partA_trimmed.fasta` |
| IQ-TREE (command line) | Maximum-likelihood tree with bootstrap | `results/partA_ml.*` |
| `07_plot_tree.R` | Plot the tree | `figures/partA_tree.png` |
| `08_diversity.R` | Diversity statistics and PhiST, all 316 sequences | `results/diversity_by_population.csv`, `pairwise_phist.csv` |
| `09_figures.R` | PCoA and pairwise differences, all 316 sequences | `figures/partB_diversity.png` |
| `10_divergent.R` | Flag divergent sequences; analysis without them | `results/divergent_sequences.csv`, `*_without_divergent.csv` |
| `11_final_analysis.R` | Exclude the *H. placei*-like cluster; final statistics and figure | `data/processed/partB_final.*`, `results/final_*`, `figures/partB_final.png` |

## 8. Repository structure

    data/raw/         sequences as downloaded
    data/processed/   filtered, aligned and trimmed data
    scripts/          R scripts, numbered in order of use
    results/          tree files and statistics tables
    figures/          tree and diversity figures

## 9. References

- Katoh K, Rozewicki J, Yamada KD (2019) MAFFT online service: multiple sequence alignment, interactive sequence choice and visualization. *Briefings in Bioinformatics* 20:1160-1166.
- Minh BQ et al. (2020) IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era. *Molecular Biology and Evolution* 37:1530-1534.
- Wong TKF et al. (2026) IQ-TREE 3. *Molecular Biology and Evolution*. doi:10.1093/molbev/msag117
- Kalyaanamoorthy S et al. (2017) ModelFinder: fast model selection for accurate phylogenetic estimates. *Nature Methods* 14:587-589.
- Minh BQ, Nguyen MAT, von Haeseler A (2013) Ultrafast approximation for phylogenetic bootstrap. *Molecular Biology and Evolution* 30:1188-1195.
- Paradis E (2010) pegas: an R package for population genetics with an integrated-modular approach. *Bioinformatics* 26:419-420.
- Paradis E, Schliep K (2019) ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. *Bioinformatics* 35:526-528.
- Yu G et al. (2018) Two methods for mapping and visualizing associated data on phylogeny using ggtree. *Molecular Biology and Evolution* 35:3041-3043.
- Tajima F (1989) Statistical method for testing the neutral mutation hypothesis by DNA polymorphism. *Genetics* 123:585-595.
- Excoffier L, Smouse PE, Quattro JM (1992) Analysis of molecular variance inferred from metric distances among DNA haplotypes. *Genetics* 131:479-491.
EOF