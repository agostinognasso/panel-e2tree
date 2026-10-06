# DSS_ORBIS — data provenance and files

> **Not in this archive.** The `.rds` files described below (`model_data.rds`,
> `harmonized_panel.rds`, `harmonized_cross.rds`, `orbis_eu_pool.rds`) and the raw
> `.tsv` extracts are derived from Bureau van Dijk ORBIS, a licensed proprietary source,
> and are not redistributed. This file documents how they are built so a reader with an
> ORBIS subscription can rebuild them. The table below describes the files as they exist
> in the authors' working tree.

Application: panel e2tree on a European firm panel (Bureau van Dijk ORBIS,
Dec-2025 snapshot). Target: the asset-based solvency ratio, shareholders' funds
over total assets, a high-ICC structural buffer used in early-warning and distress
analysis.

## Reproducibility chain

```
ORBIS raw .txt (external, licensed)         code/00_extract_orbis.sh
  Industry_Global_financials_and_ratios-EUR   ───────────────►  data/raw/output/orbis_eu_*.tsv
  Legal_info, Industry_classifications                            (NOT shipped; see licensing)
        │
        │  code/01_build_orbis.R  (FULL build, reads the .tsv)
        ▼
  data/orbis_eu_pool.rds   full balanced complete-case European pool
        │
        │  code/01_build_orbis.R  (FAST build, reloads the pool; identical output)
        ▼
  data/model_data.rds  +  harmonized_panel.rds  +  harmonized_cross.rds
```

`01_build_orbis.R` runs either way: if the raw `.tsv` extracts are present it does
the FULL build and (re)writes `orbis_eu_pool.rds`; if they are absent it does the
FAST build from `orbis_eu_pool.rds`. Both produce a **bit-identical**
`model_data.rds` (seed fixed at 1).

## Files in this folder

| File | Rows | What it is |
|---|---|---|
| `model_data.rds` | 8,000 | **Main modelling input.** Balanced panel: 800 firms × 10 years (2013-2022), complete-case. Columns: `bvdid`, `country`, `sector` (NACE section), `year`, target `solvency_asset`, + 11 predictors. |
| `harmonized_panel.rds` | 8,000 | Same panel + Mundlak between/within columns (`*_bw`, `*_wn`) for every model variable. |
| `harmonized_cross.rds` | 800 | Per-firm cross-section: long-run means + structural labels (country, NACE, legal form, age, n_years). |
| `orbis_eu_pool.rds` | 2,278,440 | **Full balanced complete-case European pool** (227,844 firms × 10 years) the sample is drawn from. Kept so the analysis reproduces without the 22 GB raw files, and for robustness / re-sampling. |
| `raw/output/orbis_eu_*.tsv` | — | Pre-filtered raw ORBIS extracts. **Deleted to save space and because ORBIS is licensed** (see below). Regenerate with `code/00_extract_orbis.sh`. |

## Universe and construction

- **Countries**: EU-27 + EFTA/EEA (CH, NO, IS, LI) + UK = 32 ISO-2 prefixes.
  (After the balanced complete-case filter, Southern Europe — IT, ES, PT —
  dominates the sample: this reflects ORBIS coverage of unconsolidated SME
  accounts with full ratio reporting, not a sampling choice.)
- **Accounts**: unconsolidated (`U*`), 12-month fiscal years, EUR-denominated
  (cross-country comparable), closing years 2013-2022.
- **Balanced panel**: only firms with all 10 years (lifts target ICC above the
  0.8 necessity threshold and gives a clean within decomposition).
- **Winsorisation**: 1%/99% on all heavy-tailed ratios.

## Target and necessity (ICC)

Variance is decomposed one-way (between firms vs within firm over time). ICC is
the between share of total variance, σ²B/(σ²B+σ²W) — the single convention used
across the pipeline and the manuscript. The method (panel e2tree) is needed
precisely when the target ICC is high (> 0.8):

| | ICC (full pool) | ICC (sample) |
|---|---|---|
| **solvency_asset (target)** | **0.826** | **0.839** |
| predictors (roa, margins, liquidity, efficiency, …) | 0.41 – 0.72 | 0.39 – 0.75 |
| ln_assets, firm_age (structural controls) | ~0.97 | ~0.97 |

The target is strongly between-dominated while the predictors carry real within
variation — the regime in which a pooled SHAP / single surrogate conflates the
two. Interest cover, the other distress candidate, was rejected (ICC ~ 0.35).
See `output/icc_by_var.csv`.

## Licensing

ORBIS / Bureau van Dijk data are proprietary and licensed. The raw extracts
(`data/raw/output/orbis_eu_*.tsv`) and the derived `.rds` files must not be
redistributed in a public replication package. Reproduction requires a valid ORBIS
subscription; `code/00_extract_orbis.sh` documents exactly how the extracts are
produced from it. Only the code is publicly shareable.
