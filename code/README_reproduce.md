# Reproducing the results

Every quantity in the paper comes out of one scripted pipeline. The scripts fit the models,
write the figures and tables, and emit the numerical macros the text uses. The prose cannot
drift from the computations: `manuscript/manuscript_DSS_v2/numbers.tex` is generated, never
edited by hand.

The firm-level data are not part of this snapshot. ORBIS (Bureau van Dijk) is licensed, so
neither the raw extracts nor the derived `.rds` files can be redistributed. Running the
pipeline from scratch needs a valid ORBIS subscription; `code/00_extract_orbis.sh` documents
exactly which fields are pulled and how. What is shipped is the code plus the aggregate
results it produced, so the published numbers can be checked without the data.

## Software dependency: the `e2tree` package and the bundled panel method

The pipeline uses the e2tree primitives (`e2tree()`, `createDisMatrix()`, `vimp()`,
`get_ensemble_predictions()`, `plot_e2tree()`, proximity/dissimilarity):

- Aria, Gnasso, Iorio, Pandolfo (2024), *Explainable Ensemble Trees*, Computational Statistics 39(1).
- Aria, Gnasso, Iorio, Fokkema (2026), *Extending Explainable Ensemble Trees to regression
  contexts*, Applied Stochastic Models in Business and Industry 42(1).

Install the pinned version used for the paper before running:

```r
remotes::install_github("agostinognasso/e2tree@644d06f")
```

> Use this commit, not CRAN. The numbers in the paper come from the regression e2tree engine
> at `644d06f`. The current CRAN release (`e2tree` 1.2.0) uses an earlier regression engine
> that does not reproduce them: the augmented surrogate degenerates and pooled fidelity drops
> markedly. Pinning to `644d06f` reproduces every number bit for bit (verified).

The panel method, `panel_e2tree()`, is the contribution of this paper. To keep it visible and
the pipeline self-contained, it is bundled here as a standalone source file,
[`code/panel_e2tree.R`](panel_e2tree.R), loaded right after the package (it shadows the
package copy with identical code):

```r
library(e2tree)
source("code/panel_e2tree.R")        # panel_e2tree() and its e2panel methods
source("code/e2tree_split_local.R")  # pure-R candidate-split enumerator for the S3 diagnostic
```

Those two files depend only on the exported e2tree API plus CRAN packages, so the snapshot
stays self-contained: no script reads or writes anything outside the project root.

CRAN packages used across the pipeline: `dplyr`, `tidyr`, `ggplot2`, `ranger`,
`randomForest`, `rpart`, `rpart.plot`, `kernelshap`, `data.table`, `gridExtra`, `patchwork`,
`partitions`. Exact versions are pinned in `../output/repro/sessionInfo.txt`.

## Run

Run everything from the project root. Each script finds the project base as the parent of its
own `code/` folder, so paths resolve whatever the working directory.

```bash
# 1) main results and every macro of the main text (caches fitted objects to
#    output/repro/cache/; set REPRO_FORCE=1 to recompute from scratch)
Rscript code/03_panel_e2tree_orbis.R

# 2) SHAP analysis, between vs within (a few minutes, kernelshap)
Rscript code/17_shap_panel.R

# 3) the two revision runs: common leaf budget, and the firm-age robustness refit
Rscript code/05_budget59.R             # output/budget59/budget59_grid.csv
Rscript code/06_within_no_age.R        # output/budget59/within_no_age.csv

# 4) Proposition 1 on the pooled-target variant (refits S5: its predictions
#    were never cached; POOLED_TARGET_FORCE=1 to refit from scratch)
Rscript code/10_pooled_target_decomp.R

# 5) turn the computed numbers into LaTeX
Rscript code/numbers_to_tex.R          # manuscript/manuscript_DSS_v2/numbers.tex
Rscript code/08_numbers_v2.R           # numbers_budget.tex + supp/tab_budget, tab_noage
Rscript code/supplement_tables.R       # the remaining supp/*.tex

# 6) publication figures (read cached outputs, no refitting)
Rscript code/04_figures_trees.R        # between/within tree diagrams
Rscript code/18_figures_shap.R         # shap_panel_A.png, shap_panel_B.png
Rscript code/07_figures_v2.R           # e2tree_necessity.png, shap_panel_B.png (v2 captions)
```

Order matters in two places. `07_figures_v2.R` must run after `03` and after `05`, because it
overwrites `e2tree_necessity.png` with the version carrying the 59-leaf budget row. It must
also run after `18`, since both write `shap_panel_B.png` and the v2 caption is the published
one. `08_numbers_v2.R` needs `05` and `06` to have run.

Two scripts sit outside this sequence and need the raw ORBIS extracts:
`01_build_orbis.R` rebuilds the panel itself (FAST from the cached pool, or FULL from the raw
extracts), and `02_compare_ideas.R` reproduces the target-selection screening.
`09_leading_contrast.R` is a standalone numerical check of the leading-contrast lemma; it
writes to `output/budget59/` and feeds no macro.

Diagnostic plots from `01`, `02` and `17` go to `output/plots/`. The manuscript figures live
in `manuscript/manuscript_DSS_v2/figures/`.

## Outputs

- `../output/repro/paper_numbers.csv` — master key/value list, the source of `numbers.tex`.
- `../output/repro/*.csv` — one CSV per paper table.
- `../output/repro/cache/` — fitted objects created on the first run. `S1_main.rds` is the
  panel-e2tree object behind the trees and the numbers. Delete the cache, or set
  `REPRO_FORCE=1`, to recompute from scratch.
- `../output/budget59/`, `../output/pooled_target/` — results of the revision runs.
- `../manuscript/manuscript_DSS_v2/figures/*.png` — the figures the manuscript includes.
- `../output/repro/sessionInfo.txt` — the environment of the reported full run.

## Pipeline at a glance

```
data/orbis_eu_pool.rds (cached European pool; or raw .tsv via 00_extract_orbis.sh)
  └─ 01_build_orbis.R ─► data/model_data.rds, harmonized_panel.rds, harmonized_cross.rds
       ├─ 03_panel_e2tree_orbis.R ─► output/repro/*  (paper_numbers.csv, table CSVs,
       │                             cache/S1_main.rds, sessionInfo); sources panel_e2tree.R
       ├─ 05_budget59.R, 06_within_no_age.R ─► output/budget59/*.csv
       ├─ 10_pooled_target_decomp.R ────────► output/pooled_target/*.csv
       └─ 17_shap_panel.R ─────────────────► output/shap_*.csv  (SHAP between vs within)
            ├─ numbers_to_tex.R, 08_numbers_v2.R ─► manuscript_DSS_v2/numbers*.tex
            ├─ supplement_tables.R ──────────────► manuscript_DSS_v2/supp/*.tex
            └─ 04, 18, 07 ──────────────────────► manuscript_DSS_v2/figures/*.png
```
