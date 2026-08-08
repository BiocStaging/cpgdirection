# Pushing this to GitHub

This folder is the package repository, and nothing else. No build artefacts, no
RStudio state, no tarballs, no analysis outputs.

## Before you push

Run these three in RStudio with the working directory set to this folder. I
could not run them — there is no R in my environment — so this is the first time
the package will have been checked by R rather than by inspection.

```r
devtools::document()   # regenerates NAMESPACE and man/ from the roxygen comments
devtools::test()
devtools::check()
```

`document()` may rewrite `man/*.Rd` and `NAMESPACE`. That is expected and
correct: the `.Rd` files here were written by hand to match the roxygen
comments, and roxygen is the authority. Commit whatever it produces.

Expect `check()` to raise notes about the installed size — see below — and
possibly about the Bioconductor annotation packages in `Imports`. Errors and
warnings should be zero; if they are not, send me the output rather than
pushing.

## Creating the repository

```bash
cd "GITHUB_cpgdirection"
git init
git add .
git commit -m "cpgdirection 2.1.1: CpG to gene-expression direction, four evidence layers"
git branch -M main
git remote add origin https://github.com/teindor/cpgdirection.git
git push -u origin main
```

Create the repository on github.com first, **empty** — no README, no .gitignore,
no licence, since this folder already has all three and GitHub's versions would
conflict on the first push.

## The one real problem: size

`inst/extdata/` is 58 MB across five files:

| file | size |
|---|---|
| `lookup_blood_hg19.csv.gz` | 19 MB |
| `lookup_nasal_epithelium_hg19.csv.gz` | 16 MB |
| `lookup_solid_tissue_hg19.csv.gz` | 8.8 MB |
| `cpg_positions_hg19.csv.gz` | 7.7 MB |
| `smr_directions.csv.gz` | 6.8 MB |

Every file is under GitHub's 50 MB warning threshold and well under the 100 MB
hard limit, so **this will push successfully as it stands**. Two consequences to
be aware of:

1. Clones are 58 MB. Tolerable, not pleasant.
2. Git stores every version of every file forever. If you later rebuild the
   lookup tables, the repository grows by another 58 MB and never shrinks. Three
   or four such rebuilds and the repository becomes genuinely awkward.

If the tables are likely to be rebuilt, use Git LFS from the start — retrofitting
it means rewriting history:

```bash
git lfs install
git lfs track "inst/extdata/*.csv.gz"
git add .gitattributes
```

If they are frozen at publication, plain git is fine and simpler.

## What is deliberately not here

- `*.tar.gz` — built artefacts belong in a Release, not in the tree
- `.Rproj.user/`, `.Rhistory`, `.RData` — local state
- the analysis code that *built* the tables (`smr_layer/`, `multi_cohort/`,
  `godmc_scan/`) — that is the paper's deposit, not the package

That last one matters for the paper. Genome Biology requires the analysis source
code as well as the tool. Keep them separate: this repository is the tool, and
the OSF deposit is the analysis.

## After the first push

Two things stop being fictional once the repository exists:

1. `README.md` tells users to run
   `BiocManager::install("teindor/cpgdirection")`. Until now that returned
   HTTP 404.
2. The manuscript can state a real code-availability URL.

Consider tagging the version described in the manuscript:

```bash
git tag -a v2.1.1 -m "version described in the manuscript"
git push origin v2.1.1
```

and connecting the repository to Zenodo, which then mints a DOI for each
release. Genome Biology recommends exactly that: a public repository for
development, plus a static DOI-assigned deposition of the published version.
