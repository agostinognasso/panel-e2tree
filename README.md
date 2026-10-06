# Panel e2tree: ORBIS application

Code and reproducibility materials for the empirical application in

> Massimo Aria, Agostino Gnasso, Carmela Iorio,
> *An explainable decision-support framework for machine-learning models on
> longitudinal data*.

Archived at <https://doi.org/10.5281/zenodo.21301342>.

The paper proposes `panel_e2tree()`, an explanation protocol for tree ensembles fitted to
panel data. The feature representation is split into a between part (unit means) and a
within part (unit-demeaned deviations), and each part gets its own Explainable Ensemble Tree
surrogate. The application is a European firm panel from Bureau van Dijk ORBIS, with the
asset-based solvency ratio (equity over total assets) as the target.

## What is in here, and what is not

The firm-level data are not included. ORBIS is a licensed, proprietary source, so neither
the raw extracts nor the derived `.rds` files can be redistributed. Reproducing the pipeline
end to end requires a valid ORBIS subscription.

What is included is the code, plus the aggregate results it produced: the master key/value
list of every number quoted in the paper, one CSV per table, the diagnostic plots, and the
`sessionInfo()` of the run that produced them. No firm-level record appears anywhere in
this archive. The point is that a reader without an ORBIS licence can
still check that the published numbers come out of this code.

```
code/        the pipeline, from the ORBIS extraction script to the LaTeX macros
data/        DATA_README.md (provenance, universe, ICC table) and an empty landing zone
output/      the aggregate results: paper_numbers.csv, one CSV per table, sessionInfo
reports/     notes on data provenance and on how the target was selected
```

`code/README_reproduce.md` has the run order, the pinned `e2tree` commit, the package list
and the two places where script order matters.

## Quick start

```r
remotes::install_github("agostinognasso/e2tree@644d06f")   # the pinned engine, not CRAN
```

```bash
Rscript code/03_panel_e2tree_orbis.R   # main results and every macro of the main text
Rscript code/numbers_to_tex.R          # manuscript/manuscript_DSS_v2/numbers.tex
```

Without `data/model_data.rds` the first script stops at its `stopifnot()` on the missing
file. Build that file first with `code/01_build_orbis.R`, which needs the ORBIS extracts.

## Reproducibility discipline

The manuscript quotes no number that is typed by hand. Every value in the text is a LaTeX
macro generated from `output/repro/paper_numbers.csv`, which `03_panel_e2tree_orbis.R`
writes. Text and pipeline therefore cannot drift apart. The same holds for the revision runs:
`05`, `06` and `10` write their own CSVs, and `08_numbers_v2.R` and
`10_pooled_target_decomp.R` turn those into `numbers_budget.tex` and
`numbers_pooled_target.tex`.

Fitted objects are cached under `output/repro/cache/`, so a second run is fast. Set
`REPRO_FORCE=1` to recompute from scratch.

## Licence

The code is released under the MIT licence (see `LICENSE`). The ORBIS data it reads are
proprietary to Bureau van Dijk and are covered by their own licence terms, not by this one.

## Citing

Please cite the article. To cite the code itself, use the archived snapshot at
<https://doi.org/10.5281/zenodo.21301342>.
