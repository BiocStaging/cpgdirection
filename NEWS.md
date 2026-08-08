# cpgdirection 2.2.4

* **The SMR layer could never fire unless a gene was named.** Its join keyed on
  `requested_gene` and nothing else, so calling
  `cpg_expression_direction(my_cpgs)` — the mode most people use — left the key
  `NA` and the 413,982-pair table contributed nothing. No error, no warning, and
  the printed summary blamed 450K array coverage for the silence.

  The layer is now always reported, with two new columns deciding what may be
  done with it. `smr_gene` names the gene the causal estimate concerns.
  `smr_gene_match` says whether that is the same gene as the rest of the row.

  This matters more than it sounds, because the two layers usually disagree about
  which gene a CpG regulates: of the 85,220 CpGs present in both the SMR table
  and a catalogue, only 58.6% share even one gene, and since the catalogue commits
  to a single target the realised match rate is lower still. `cg00000029` is
  typical — SMR reaches AKTIP and RBL2, the blood catalogue reaches SCP2.

  So the rule is: report always, promote only on match. A non-matching row shows
  its direction, tier and p-value in the `smr_*` columns but cannot supply
  `best_direction`, because that column answers a question about one specific
  pair. `smr_agreement` is `NA` rather than `FALSE` for those rows: two statements
  about different genes cannot contradict each other, and calling it disagreement
  would manufacture a conflict out of compatible facts.

  Where a gene *is* supplied, catalogue-derived fallback keys are not used at all.
  If you ask for CREBBP, a direction for ADCY9 is not an answer to your question.

  No validation figure changes: S1 95.9% / S2 84.9% / S3 70.4% were measured
  against the table directly. What changes is that they now reach the caller.

* **Two-token identifiers lost their gene.** `.cpgd_parse_input()` required more
  than two underscore-separated fields and discarded the first two by position,
  so `"cg00039463_TRAP1"` — the plainest form anyone would type — parsed to a CpG
  with no gene, and the distance layer reported itself unavailable. Gene
  candidates are now taken from whatever follows the token carrying the `cg`
  identifier, with probe-type codes (`TC21`, `BC11`) dropped by name rather than
  by position.

* **The "no SMR evidence" message asserted a cause it had not checked.** It
  attributed every empty result to GoDMC's 450K coverage, including for CpGs
  sitting in the table with four significant rows. It now distinguishes a CpG
  absent from the table from one present but unmatched on gene, and says which.
  New column `smr_in_table` exposes the same distinction.

* Regression tests in `test-smr-layer.R` covering all three. Each would have
  failed on 2.2.3, and none of the existing tests could have caught them: every
  bug returned a plausible answer rather than an error.

* `DESCRIPTION` coverage figures corrected. They still described a single
  gradient-boosted table of 1,603,777 pairs over 604,456 CpGs, which predated the
  nasal and SMR layers. The shipped data are 3,105,697 pairs over 761,445 CpGs
  across four evidence sources. Analysis code and fitted models are archived at
  <https://doi.org/10.17605/OSF.IO/U3VFK>.

# cpgdirection 2.2.3

* `UNMAPPED` was unreachable. An identifier absent from the lookup table came
  back as `NO_EXPRESSION_ANCHOR` instead, because the right join fills
  non-matching rows with `NA` and the code then read a missing `target_gene` as
  "no gene nearby". The two are different facts and a user needs to tell them
  apart: an absent identifier usually means a typo or a probe from a different
  array, while no expression anchor is a real biological answer about a real
  CpG. The distinction is now captured before `status` is overwritten, and
  `UNMAPPED` rows carry a note saying to check the probe ID and array version.

  Found by a test that had been asserting the documented behaviour since 1.0.0
  and had never once been run.

* No change to any prediction. Affects only how absent identifiers are labelled.

# cpgdirection 2.2.2

* Fixes `could not find function "fifelse"`, introduced by me in 2.2.1. Deleting
  the hand-written `NAMESPACE` also deleted its
  `importFrom(data.table, ...)` list, and the replacement declared only `:=`
  and `.SD`. Three functions are called bare inside `[.data.table` expressions,
  where a `::` prefix would be evaluated in the wrong frame: `fifelse`, `fcase`
  and `setorderv`. All three are now declared.

  The scan I used to find them the first time had a hole: it excluded any name
  used both bare and qualified, which is exactly the case here, since
  `data.table::fifelse` appears elsewhere in the package. Corrected and re-run
  against the full data.table export list.

* `cpgd_lookup_info()` gains the `@param tissue` documentation it never had,
  clearing the last `R CMD check` warning.

* No change to any prediction.

# cpgdirection 2.2.1

The first `R CMD check` the package has ever had. Two errors, two warnings and
three notes, all of them real, none of them affecting a prediction.

* **Documentation is now generated by roxygen rather than written by hand.** The
  hand-written `man/*.Rd` and `NAMESPACE` had drifted from the code -- check
  reported codoc mismatches where `cpgd_build_gene_tss()` and `cpgd_selftest()`
  were documented as taking no arguments -- and roxygen refused to touch files it
  had not generated, so every `devtools::document()` was a no-op. Both are
  deleted and regenerated. A new `R/cpgdirection-package.R` carries the imports.

* **Four failing tests, one cause.** They predated the multi-layer default and
  asserted on `gene_source`, `status` and `call`, which exist only when a tissue
  is named; the default returns the `best_*` column set. They now ask for
  `tissue = "blood"` where they mean the per-tissue shape. Two tests added: one
  checking the default column set, one checking that `target_tissue` fills
  blanks without overwriting existing calls.

* **A test of mine was nonsense** -- it called `readLines()` on a directory. Its
  intent, that every label the ladder emits is a known level, is now tested
  directly, along with the ladder ordering.

* **The example errored** because it selected `target_gene` from a result that
  has no such column in the default mode. It now shows both shapes.

* Fixed: non-ASCII characters in `R/report.R` code (permitted in comments, not
  in code); the `cpgd_report()` example writing into the check directory;
  `PUSH_INSTRUCTIONS.md` flagged as a non-standard top-level file; and the
  "declared imports not used" note for the two annotation packages, which are
  reached through `asNamespace()` and so were invisible to check.

* Remaining and expected: installed size 58 Mb, which is the lookup tables.

# cpgdirection 2.2.0

The Mendelian randomisation layer worked but was invisible, which made a working
layer indistinguishable from a broken one. Three changes, no change to any
prediction.

* **New `smr_agreement` column.** SMR frequently loses the evidence ladder to a
  catalogue consensus, and until now it then vanished from the output entirely —
  including where it *disagreed*. In a 65-CpG HPA panel the catalogue reported
  `-1` for cg16302441/POMC in all three tissues while SMR reported `+1` at
  p = 3.8e-8, and nothing in the output said so. A correlational prediction and
  a causal estimate pointing opposite ways is a finding, not a tie to be settled
  silently by precedence. Disagreements are now flagged in the column, in `note`
  and in the printed summary.

* **The printout reports the SMR layer even when it wins nothing**, giving the
  number of CpGs with SMR evidence, how many supplied `best_direction`, and how
  many were outranked. Where a set has no SMR evidence at all it says so and
  says why: GoDMC is a 450K meta-analysis, so EPIC-only probes were never
  assayed and cannot have an mQTL.

* **When `best_direction` contains blanks, the printout now gives the exact call
  that fills them** rather than leaving the reader to find `target_tissue` in the
  help page.

# cpgdirection 2.1.1

* Fixes `factor level [11] is duplicated`, which broke `cpg_expression_direction()`
  in 2.1.0. The print method held the evidence levels as an inline vector that
  had been hand-edited each time a layer was added; the 1.9.0 and 2.1.0 patches
  each appended `distance_targeted` and `catalogue_targeted`, and `factor()`
  rejects a duplicated level.

* The levels now live in a single constant, `CPGD_EVIDENCE`, in ladder order,
  with a `stopifnot(!anyDuplicated(.))` beside it and a regression test. The
  print method also unions in any level it sees but does not recognise, so a
  future layer cannot be silently dropped from the tally. The underlying cause
  was one list maintained in two places; the fix is that there is now one.

* No change to any prediction. 2.1.0 computed results correctly and failed only
  when printing them.

# cpgdirection 2.0.0

**A causal evidence layer, from summary-data Mendelian randomisation.**

The package no longer depends only on what someone happened to measure. A new
layer asks whether a genetic instrument for methylation at a CpG moves
expression of a gene -- so its coverage is not bounded by which pairs a study
chose to test.

* **413,982 CpG-gene pairs, 103,285 CpGs, 14,056 genes.** Exposure: GoDMC blood
  mQTL meta-analysis (~32,000 individuals). Outcome: eQTLGen blood cis-eQTL
  (31,684 individuals). Exposed as `cpgd_smr_directions()` and folded into
  `best_direction`.

* **Validated, not asserted.** Against the directly measured eQTMs already
  shipped, on the pairs where both exist:

  | tier | pairs | validated on | concordance |
  |---|---:|---:|---:|
  | S1 (>1 instrument, all agreeing) | 29,456 | 2,141 | **95.9% ± 0.8** |
  | S2 (single instrument) | 366,838 | 6,008 | 84.9% ± 0.9 |
  | S3 (instruments disagree) | 17,688 | 456 | 70.4% ± 4.2 |

  Majority-class baseline on the same overlap: 59.8%. S1 beats tier A of the
  catalogue models and approaches directly measured evidence. S3 scoring *below*
  S2 is the expected ordering, and a small check on the method: instrument
  disagreement is a warning, not extra evidence.

* **The ladder is reordered by measured accuracy, not by recency**: measured,
  `smr_high`, `catalogue_consensus`, `smr_moderate`, `catalogue_single`,
  `smr_weak`, then the distance layers.

* **It reaches where the distance curves cannot.** Median CpG-gene distance in
  the SMR table is 135 kb, 90th percentile 657 kb -- precisely the range where
  the distance curves flatten and the catalogue models have no support.

* The direction distribution is **49.2% positive**, against 26.2% from the
  catalogue models. SMR is not thresholded on an imbalanced classification
  problem, so it does not shrink towards the majority class.

* **What SMR does not establish.** A significant result is consistent with
  methylation affecting expression and equally with **linkage**: two distinct
  causal variants in LD, one acting on each trait. HEIDI is not run. Multiple
  concordant instruments (S1) argue against linkage without excluding it. Both
  datasets are whole blood. Palindromic SNPs and unresolvable allele mismatches
  were dropped rather than assumed -- 21,698 comparisons refused on those
  grounds, 9 on allele mismatch.

# cpgdirection 1.9.0

* New `target_tissue` argument. Where the tissues disagree, conflicts are
  resolved using a named tissue rather than abstaining, producing
  `distance_targeted` and `catalogue_targeted` rows. This exists because blood
  holds an effective veto over distance calls -- unanimity is required and blood
  never votes `+1` -- which is correct when the question is about blood and
  wrong when it is not. A saliva panel of brain-expressed genes is not asking
  about blood.

  It is tissue **selection, not voting**. Name the tissue on biological grounds
  before looking at the votes. A majority rule tuned after seeing which way each
  tissue fell is selecting on the outcome, and would be indefensible at review;
  a tissue named in advance is a stated assumption that can be argued with.

  The price is carried in the output: targeted rows report the declared tissue's
  accuracy, not the consensus figure, and stay flagged in
  `direction_uncertain` so the sensitivity analysis still covers them. For
  `solid_tissue` that accuracy is 0.55-0.70, tumour tissue, provisional, with
  brain and kidney below their own majority baselines under leave-one-cohort-out
  -- the clause to quote in a limitations section if solid tissue is being used
  as a brain proxy.

# cpgdirection 1.8.1

* States plainly what 1.8.0 made structural: the blood distance curve peaks at
  **0.449 and never reaches 0.5 at any distance**, and a distance call requires
  all three tissue curves to agree. `distance_only` therefore **cannot return
  `+1` for any CpG, ever**. Every direction that layer contributes is `-1` by
  construction. That is a property of the fitted curves, not of any user's
  panel, and it is the whole explanation for a distance-dominated panel looking
  uniformly negative. Now said in the help page and in
  `cpgd_direction_balance()`.

* Corrected the `cpgd_direction_balance()` footer, which still described the
  pre-1.8.0 pooled-curve rule.

# cpgdirection 1.8.0

**The distance layer no longer averages away tissue disagreement.**

Prompted by the observation that a real HPA panel returned 43 directions of
which every single one was `-1`. Investigating that produced a genuine defect,
though not the one expected.

* The three distance curves are not parallel. Blood asymptotes at 0.449 and
  **never crosses 0.5 at any distance**; nasal epithelium reaches 0.594 and
  solid tissue 0.861. Beyond roughly 20 kb they disagree, often sharply.
  Averaging them yielded a pooled probability near 0.50, which the previous
  version read as "no information" and reported as `distance_uninformative`.
  It was not ignorance: blood was saying `-1` and solid tissue `+1`, both
  clearly.

* A distance call now enters `best_direction` only when **all three curves agree
  in sign**. Where they disagree the row is `distance_tissue_conflict`, with new
  `dist_dir_blood`, `dist_dir_nasal` and `dist_dir_solid` columns so the tissue
  relevant to the question can be read directly. This is the rule the catalogue
  layer already used, now applied consistently.

* On the 65-CpG HPA panel: 45 unanimous `-1` calls and 20 tissue conflicts,
  against the previous 43 calls and 22 uninformative. **`best_direction` gains
  no new `+1`** — the change buys transparency about where the tissues part
  company, not extra positives. In 19 of the 20 conflicts solid tissue says
  `+1` while blood says `-1`.

* `min_distance_info` is demoted from the main gate to a floor on the pooled
  probability, default lowered from 0.10 to 0.02. Unanimity across three
  independently derived curves is evidence in its own right and is no longer
  suppressed for having a modest magnitude.

* **Known limitation, documented not silently fixed.** The isotonic fit is
  monotone by construction and in three bins it pulls a raw estimate across 0.5:
  nasal at 41,827 bp (raw 0.5404, fitted 0.4884), nasal at 69,506 bp (raw
  0.5053, fitted 0.4901) and solid tissue at 25,171 bp (raw 0.5216, fitted
  0.4979). The fitted value also falls outside the raw Wilson interval in nine
  bins overall. The fit is the pre-registered method and has not been changed;
  `distance_curves.csv` ships `p_positive_raw` alongside `p_positive` so the
  discrepancy can be inspected.

# cpgdirection 1.7.1

* `cpgd_selftest()` now passes. `fread` returns an empty character field as `""`
  rather than `NA`, while the package returns `NA_character_`, so every
  abstaining row was scored a mismatch on `target_gene` and `evidence_tier` --
  exactly the 12 of 40 that produced the 70% figures. Numeric fields were
  unaffected, which is why `probability_plus1` and `confidence` already read
  100%. Blanks are now normalised before comparison.

* Removed the remaining data.table shallow-copy warning. It came from
  `cpg_direction_universal()`, which built its per-tissue probability columns
  with `[[<-`; the 1.6.1 fix had only covered the equivalent code in the
  multi-tissue assembly.

# cpgdirection 1.7.0

* Two alternative codings for rows whose sign is nearly free to flip.
  `best_direction` still abstains on them, which is honest but leaves nothing to
  analyse. `best_direction_filled` now takes the point estimate for those rows
  and `best_direction_flipped` takes its opposite, with `direction_uncertain`
  flagging which rows differ.

  **They are a sensitivity analysis, not a menu.** The question is whether a
  downstream result survives both codings, not which coding produces one.
  Choosing the column that fits the data better is selecting on the outcome: it
  doubles the effective number of tests and biases towards significance. The
  printout, the HTML report and the help page all say so.

  A row counts as uncertain when the distance prior sits within
  `min_distance_info` of 0.5, or when a call was made below that confidence.

# cpgdirection 1.6.1

* `cpgd_selftest()` no longer reports `NA%` and fails. It compared directions
  with `==` and let `NA` propagate, so a single row missing on one side turned
  the whole score into `NA`. Missing-on-one-side now counts as a mismatch.

* The test now gates on `probability_plus1`, the model's actual output, together
  with gene assignment and evidence tier, rather than on direction. Direction is
  the sign of the probability, so testing it added nothing and was the one field
  whose *reporting* convention had changed. Verified against the shipped tables:
  probability, confidence, tier and TSS distance agree exactly on all 28
  reference CpGs present in the blood table.

* `cpgd_selftest(regenerate = "file.csv")` writes a fresh reference from the
  installed package, for use after a deliberate change of behaviour.

* Removed a spurious data.table shallow-copy warning. The per-tissue assembly
  used `[[<-`, which is base-R assignment and copies the table, so every later
  `:=` warned about a copy the caller had nothing to do with.

# cpgdirection 1.6.0

* `cpg_expression_direction(..., html_out = "results.html")` writes a formatted
  HTML table on completion, into `getwd()` unless given an absolute path, so
  `setwd()` controls where results land. `cpgd_report()` does the same for a
  result you already have. Uses **gt** when installed and falls back to a plain
  styled table when it is not, so a long run never fails at the last step for
  want of an optional package. Direction cells are shaded by sign and
  `best_evidence` by strength, because that is the column readers skip.

* New `cpgd_direction_balance()`. Reports the +1 / -1 balance in the shipped
  tables and, given a result, tests its catalogue-derived calls against the
  relevant base rate with an exact binomial test. Distance-only calls are
  excluded on purpose: the pooled distance curve sits below 0.5 out to ~19 kb
  and above it only past ~300 kb, so a panel of promoter-proximal CpGs cannot
  receive a +1 from that layer at all, and including those calls would make an
  entirely expected result look anomalous.

* Documented a real property of the models rather than a bug: they predict +1
  less often than it occurs. In blood, 25.8% of predictions are +1 against
  39.5% among directly measured pairs, and the same shrinkage towards the
  majority class appears in nasal epithelium (21.3% vs 38.5%) and solid tissue
  (23.8% vs 46.1%). A +1 call is therefore conservative; a -1 call carries less
  weight than its frequency suggests.

# cpgdirection 1.5.1

Three corrections found by running a real 68-CpG HPA panel through 1.5.0.

* The distance layer no longer requires the three tissue curves to agree on a
  side of 0.5. The blood curve asymptotes at 0.449 and never crosses it, while
  nasal and solid rise past it, so beyond roughly 30 kb they disagree by
  construction. Requiring consensus discarded a perfectly well-determined
  distance and reported `no_evidence` for pairs whose TSS distance was known to
  the base pair. The pooled probability is used instead, with a new
  `agree_universal` column recording whether the curves concurred.

* New `best_evidence` level `distance_uninformative`, and a
  `min_distance_info` argument (default 0.10 on the `|p - 0.5| * 2` scale).
  Where the distance prior sits at 0.47-0.51 it is not a weak answer, it is no
  answer, and labelling it `no_evidence` was as misleading as giving it a
  direction would have been.

* `best_confidence` for catalogue rows now averages over the tissues that
  actually **called**, not over every tissue holding a number. Including
  abstainers pulled single confident calls down towards the abstention
  threshold: `catalogue_single` rows were reporting a median confidence of 0.20
  when the calling tissue was well above its tier threshold.

* `best_tier` is now `U` for distance rows and `M` for measured rows rather than
  `NA`, and the printout warns when most directions rest on distance alone,
  since inverse associations are the majority class in every catalogue and those
  calls sit close to the base rate.

# cpgdirection 1.5.0

`cpg_expression_direction()` now does everything by default.

* `tissue = "all"` is the new default. One call runs the measured catalogue, all
  three tissue models, and the distance fallback, and returns them side by side
  with a single **`best_direction`** column holding the answer from the
  strongest evidence available for that CpG. `best_evidence` says which layer it
  came from, and `best_expected_accuracy` what it is worth — those two must be
  read together, since `measured` is worth ~1.00 and `distance_only` ~0.62.

* The evidence ladder, in order: a measured eQTM for the requested gene; a
  cross-tissue consensus; a single calling tissue; distance alone.

* Where tissues return opposite directions, `best_evidence` is
  `tissue_conflict` and `best_direction` is `NA`. The distance curve is
  deliberately not used to break the tie, because it is weaker than either
  tissue it would be adjudicating.

* Passing a specific tissue returns the previous single-tissue result unchanged.

* `cpg_direction_all_tissues()` is now a thin wrapper on the same engine with
  `universal = FALSE`. Its output gains `tier_*` and `best_*` columns; nothing
  was removed.

# cpgdirection 1.4.1

Two bug fixes, both affecting real panels rather than the examples.

* `annotation_mismatch` was TRUE for every CpG with **no** measured record.
  `strsplit(NA, ";")` returns a length-one element holding `NA` rather than an
  empty vector, so the emptiness guard never fired and absence of evidence was
  reported as evidence of conflict. On a 65-CpG HPA panel this flagged 63 CpGs;
  it now flags 5, all real. Absence of a measured record now returns `NA`.

* The mismatch comparison now accepts compound panel labels. A requested gene of
  `LINC02210_CRHR1` or `ENSG00000253356_STAR` is compared on the whole label,
  both separator forms, and its last token, which is where the informative
  symbol sits in Illumina-style annotations.

* `cpg_direction_universal()` no longer requires `genes`. If omitted, the symbol
  is taken from the identifier, as every other function in the package already
  did. Inputs whose symbol matches no gene in the annotation are dropped with a
  count, not an error.

* Duplicate probe identifiers (`cg09527270_TC11` / `_TC12` / `_TC13`) are still
  collapsed, but the printout now says how many were collapsed instead of
  silently reporting a smaller total.

* Removed `inst/extdata/eqtm_direction_lookup_hg19.csv.gz`, an unreferenced
  18 MB holdover from the single-tissue version, superseded by
  `lookup_blood_hg19.csv.gz`. Installed size drops by 18 MB, to about 52 MB.

# cpgdirection 1.4.0

* The four Bioconductor annotation packages moved from `Suggests` to `Imports`,
  with `biocViews:`, so that a single `BiocManager::install()` resolves
  everything. This is what makes the tarball route require the dependencies to
  be present first.

* New `cpgd_gene_tss(refresh = FALSE)`. Universal mode now builds the TSS
  annotation from TxDb on first use and caches it under
  `tools::R_user_dir("cpgdirection", "cache")` — once per machine rather than
  once per session. Previously `cpgd_build_gene_tss()` had to be run by hand and
  its CSV passed in.

No model, lookup table or calibration changed in 1.4.x. Predictions are
identical to 1.3.0.

# cpgdirection 1.3.0

* Universal mode: `cpg_direction_universal()`, `cpgd_build_gene_tss()`,
  `cpgd_cpg_positions()`. Predicts direction from distance alone, for any
  CpG-gene pair, at evidence tier U. Weaker by construction — AUC 0.616 against
  0.793 for the full blood model — and labelled as such in every row.

* Bundled 930,178 EPIC v2 probe positions in hg19.

# cpgdirection 1.2.0

* Multi-tissue: `cpg_direction_all_tissues()` returns one column per tissue plus
  `consensus_direction`, for use where the sampled tissue is a proxy rather than
  the tissue of interest.

* `cpgd_measured_eqtms()` exposes the 25,056 directly measured CpG-gene pairs,
  and results now carry `annotation_mismatch`: a flag for CpGs whose measured
  association is with a gene other than the array annotation. 23.1% of measured
  pairs involve a gene that is not the nearest measured gene.

* Separate lookup tables, models and calibration for blood, nasal epithelium and
  solid tissue.

# cpgdirection 1.1.1

* `cpgd_lookup_info()` used `tissue` in its body without accepting it as an
  argument.

# cpgdirection 1.1.0

* The identifier parser now extracts `cg[0-9]{6,}` from anywhere in a string.
  Prefixed names such as `Zcg00335286_TC21_MC2R`, produced by z-scoring a
  column, were previously dropped without warning.

* New `REQUESTED_GENE_UNAVAILABLE` status and `allow_nearest_gene` argument. When
  a gene is requested and the table holds no pairing for it, the function no
  longer silently returns a direction for the nearest catalogued gene, which is a
  different gene. `cg02079741` annotated to POMC was returning a confident
  direction for DNMT3A.

* `cpgd_parse_check()` shows what the parser extracted from each identifier.

# cpgdirection 1.0.0

* First release. `cpg_expression_direction()`, `cpgd_selftest()`,
  `cpgd_accuracy_table()`, `cpgd_lookup_info()`.
