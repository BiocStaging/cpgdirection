# cpgdirection

Predict the direction of a CpG's effect on gene expression, for CpGs with no
measured eQTM record.

```
+1  higher methylation → HIGHER expression
-1  higher methylation → LOWER  expression
```

## Install

Dependencies first, then the package. `install.packages()` on a local tarball
does **not** resolve dependencies, and the annotation packages are hard imports,
so the order matters.

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("TxDb.Hsapiens.UCSC.hg19.knownGene", "org.Hs.eg.db",
                       "GenomicFeatures", "GenomicRanges", "AnnotationDbi"))
install.packages("data.table")

install.packages("cpgdirection_2.2.3.tar.gz", repos = NULL, type = "source")
```

Budget **~400 MB and several minutes** for `TxDb.Hsapiens.UCSC.hg19.knownGene`
and `org.Hs.eg.db` the first time. They are needed only by
`cpg_direction_universal()`; the catalogue and multi-tissue functions have no
Bioconductor requirement.

Once the repository is public, `BiocManager::install("teindor/cpgdirection")`
does both steps in one line.

All three prediction tables and 930,178 EPIC v2 hg19 probe positions ship
inside the package (~52 MB). No Python, no compilation.

## Use

```r
library(cpgdirection)

cpgd_selftest()        # run once after installing

res <- cpg_expression_direction(my_cpgs)
```

One call runs everything: the measured catalogue, all three tissue models, and
the distance fallback. Read two columns first.

| column | meaning |
|---|---|
| `best_direction` | +1, −1 or `NA` — the answer from the strongest evidence available |
| `best_evidence` | which layer it came from |

`best_evidence` is the one that stops you over-reading the result:

| value | worth | meaning |
|---|---|---|
| `measured` | ~1.00 | a catalogue measured this CpG against this gene |
| `catalogue_consensus` | 0.64–0.87 | tissues that called it agreed |
| `catalogue_single` | 0.55–0.87 | one tissue only, no corroboration |
| `distance_only` | 0.60–0.65 | distance to the TSS, nothing else |
| `tissue_conflict` | — | tissues disagreed; not resolved |
| `no_evidence` | — | no catalogue holds the pair and no distance available |

A `measured` row and a `distance_only` row are not comparable. Sorting on
`best_direction` without reading `best_evidence` will mix them.

Where tissues disagree the answer is `NA`, deliberately. The distance curve is
weaker than either tissue it would be adjudicating, so it is not used to break
the tie.

For one tissue only, pass it:

```r
cpg_expression_direction(my_cpgs, tissue = "blood")
```

Identifiers are stripped to canonical form automatically, so EPIC v2 suffixes
are fine:

```r
cpg_expression_direction(c("cg00000029_TC21", "cg00050692_TC21_DNMT3A"))
```

In the second case the gene is parsed out of the name and used as the target,
which is usually what you want.

**Supply the gene when you know it.** Without one the function falls back to the
nearest transcription start site, which is frequently not the gene you had in
mind:

```r
cpg_expression_direction(cpgs  = c("cg00000029", "cg00050692"),
                         genes = c("RBL2",       "DNMT3A"))
```

The `gene_source` column records which path was taken — `user_supplied`,
`parsed_from_name`, `catalogue` or `nearest_TSS`. Treat `nearest_TSS` rows as
provisional.

A file path also works, and a `data.frame` with `cpg` and optional `gene`
columns:

```r
cpg_expression_direction("my_cpgs.txt")
```

## What you get back

One row per CpG:

| column | meaning |
|---|---|
| `direction` | +1, −1, or `NA` |
| `direction_label` | the same in words |
| `probability_plus1` | calibrated probability of a positive association |
| `confidence` | 0–1 |
| `evidence_tier` | A, B or C |
| `expected_accuracy` | what a direction at this tier is worth |
| `call` | +1, −1 or `ABSTAIN` after thresholds |
| `status` | `DIRECT_eQTM`, `PREDICTED`, `ABSTAIN`, `NO_EXPRESSION_ANCHOR`, `UNMAPPED` |
| `gene_source` | how the target gene was chosen |
| `note` | caveat in plain language |

`summary(res)` gives coverage and the direction split; `cpgd_accuracy_table()`
prints the validation figures behind `expected_accuracy`.

## What a direction is worth

| status | expected accuracy | basis |
|---|---|---|
| `DIRECT_eQTM` | ~1.00 | read from a catalogue, not predicted |
| tier A | 0.77–0.85 | locked test 0.766; independent cohort 0.854 |
| tier B | 0.64–0.85 | locked test 0.642; independent cohort 0.854 |
| baseline | 0.58–0.63 | always guessing −1 |

Tier A means catalogued CpGs lie within 50 kb **and** the target gene is itself
catalogued. Tier B means one of those. Tier C means neither, and the function
always abstains.

## Three things to know before you use this

**It gives direction, not existence.** Every accuracy figure is conditional on
the CpG being an eQTM. The catalogues behind the model contain only significant
associations, so nothing here can tell you whether a CpG affects expression at
all. Feed it a CpG with no real effect and you will get a confident direction
for something that is not there. `ABSTAIN` and `NO_EXPRESSION_ANCHOR` are
correct answers.

**Do not override the abstentions.** Tier C is badly calibrated — that is why it
abstains rather than reporting a number.

**The table is hg19.** EPIC v2 manifests ship GRCh38. Matching by probe
identifier, as this package does, is unaffected. Matching by position is not.

## Tissue

Three tissues ship, each with its own model, calibration and lookup table:

```r
cpg_expression_direction(my_cpgs, tissue = "blood")             # default
cpg_expression_direction(my_cpgs, tissue = "nasal_epithelium")
cpg_expression_direction(my_cpgs, tissue = "solid_tissue")
```

| tissue | cohorts | pairs | CpGs | grouped-CV AUC |
|---|---|---|---|---|
| blood | 4 | 1,663,799 | 641,325 | 0.793 |
| nasal_epithelium | 2 | 1,472,435 | 553,145 | 0.865 |
| solid_tissue | 6 (TCGA) | 849,229 | 454,821 | 0.664 |

**The three are not equally trustworthy, and the package says so in
`expected_accuracy`.**

*Blood* is the only arm with genuine external validation: a pre-registered
locked test set (tier A AUC 0.781, 76.6% accuracy against a 59.2% baseline) and
an independent white-blood-cell cohort (85.4%).

*Nasal epithelium* has the highest cross-validated AUC, but both contributing
cohorts come from one study and share participants, so 0.865 is within-study
consistency, not external validation. Treat it as promising rather than proven.

*Solid tissue* is tumour tissue from TCGA. Two of six tissues (brain, kidney)
did not exceed their own majority baselines under leave-one-cohort-out. It is
provisional, and brain results in particular should not be read as applying to
healthy brain.

The same CpG can return different answers in different tissues. That is the
finding, not a bug: catalogues differ in which CpGs they contain, and base
rates differ from 0.334 to 0.508 across tissues.

## Citation

[Manuscript citation and DOI to be added.]

Analysis code, fitted models and the full derivation are at
[OSF DOI to be added].
