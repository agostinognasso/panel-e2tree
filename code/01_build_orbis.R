# 01 — Building the ORBIS firm-year panel (Idea 4: early-warning of distress)
#
# Replaces the country/life-expectancy application with a corporate one built on
# Bureau van Dijk ORBIS (snapshot Dec-2025). Universe: European unconsolidated
# accounts (EU-27 + EFTA/EEA + UK; consolidation code U*, 12-month fiscal years,
# 2013-2022), EUR-denominated so countries are directly comparable. The firm
# country (BvD ID prefix) enters as a time-invariant structural (between-only)
# predictor, exactly the kind of level effect the between e2tree should isolate.
#
# Target (regression, high ICC, verified empirically):
#   solvency_asset = Solvency ratio, asset-based (%) = shareholders' funds / total
#   assets. The firm's capitalisation buffer: a structural, highly persistent
#   level (who is well- vs under-capitalised) with genuine year-on-year drift.
#   Measured ICC ~ 0.83 (> 0.8 threshold) on this panel, the regime in which a
#   pooled SHAP/surrogate conflates between & within. (Interest cover, the other
#   distress candidate, was rejected: its ICC is only ~0.35, within-dominated.)
#
# Provenance: the heavy ORBIS .txt files (tens of GB, held on a licensed local
#   export pointed to by ORBIS_ROOT) are not shipped. They are pre-filtered to
#   the European universe by code/00_extract_orbis.sh into
#   data/raw/output/orbis_eu_*.tsv. This script reduces those extracts to a
#   tractable, complete-case modelling panel and draws a reproducible firm sample
#   (seed fixed) small enough for the e2tree O(N^2) proximity matrix downstream.
#   Re-running 00 + 01 regenerates data/model_data.rds bit-for-bit.
#
# Output (under data/):
#   data/model_data.rds       complete-case firm-year modelling panel (main input)
#   data/harmonized_panel.rds sampled panel + Mundlak between/within (_bw/_wn)
#   data/harmonized_cross.rds per-firm cross-section (long-run means + labels)
#   data/orbis_eu_pool.rds    FULL complete-case pool before sampling (provenance)
# Output (under output/):
#   output/icc_by_var.csv     ICC (between share of variance) for target+features
#   output/coverage_by_var.csv, descriptives.csv, firms_by_country.csv
#   output/plots/data_*.png

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(ggplot2)
})

# ROOT = parent of code/: every path stays inside the package.
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
ROOT <- if (length(script_path)) dirname(dirname(normalizePath(script_path))) else normalizePath("..")
SRC  <- file.path(ROOT, "data", "raw", "output")
DATA <- file.path(ROOT, "data")
OUT  <- file.path(ROOT, "output")
PLT  <- file.path(ROOT, "output", "plots")
for (d in c(DATA, OUT, PLT)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
set.seed(1)

# ---- parameters -----------------------------------------------------------
YEAR_MIN  <- 2013L    # recent decade: best, most comparable ORBIS coverage
YEAR_MAX  <- 2022L
MIN_YEARS <- 10L      # require the full balanced decade 2013-2022 (10 obs/firm):
                      #   established firms -> persistent capital structure, and a
                      #   clean balanced within decomposition. Lifts target ICC to
                      #   ~0.81 (> 0.8 threshold); >=7 would give ~0.795.
N_FIRMS   <- 800L     # sampled firms -> e2tree dissimilarity matrix stays O(N^2)-feasible
WINS_P    <- 0.01     # winsorisation tail (1% / 99%)

winsor <- function(x, p = WINS_P) {
  q <- quantile(x, c(p, 1 - p), na.rm = TRUE)
  pmin(pmax(x, q[1]), q[2])
}

# ---- model variables -------------------------------------------------------
# Target = solvency ratio (equity/assets, %): high-ICC structural buffer.
# gearing is dropped from predictors (debt/equity is the near-mechanical
# complement of equity/assets); icover enters as a debt-service predictor.
target <- "solvency_asset"
feats  <- c("roa_pbt", "ebitda_margin", "cf_margin", "icover",
            "current_ratio", "liquidity_ratio", "wc_gap", "stock_turnover",
            "ln_assets", "firm_age", "rev_per_emp")
model_vars <- c(target, feats)

# ---- 1-6. build the complete-case European pool ----------------------------
# Two entry points:
#  (a) FULL  : if the raw ORBIS extracts (data/raw/output/orbis_eu_*.tsv) are
#              present, stream them and rebuild the pool from scratch.
#  (b) FAST  : otherwise reload the shipped data/orbis_eu_pool.rds, the same
#              balanced complete-case pool, so the analysis reproduces without
#              re-reading 22 GB (and without the licensed raw ORBIS files, which
#              are not redistributable). Regenerate the .tsv with 00_extract_orbis.sh.
fin_tsv  <- file.path(SRC, "orbis_eu_financials.tsv")
pool_rds <- file.path(DATA, "orbis_eu_pool.rds")

if (file.exists(fin_tsv)) {
  # FULL build: raw extracts present, rebuild the pool from them.
  fin_cols <- c("bvdid", "consol", "closdate", "total_assets", "turnover",
                "cashflow", "collection_days", "credit_days", "oprev_per_emp",
                "interest_cover", "gearing", "ebitda_margin", "current_ratio",
                "liquidity_ratio", "stock_turnover", "roa_pbt", "solvency_asset")
  fin <- fread(fin_tsv, sep = "\t", quote = "", select = fin_cols,
               na.strings = c("", "n.a.", "n.s."), showProgress = FALSE)
  fin[, year    := as.integer(substr(closdate, 1, 4))]
  fin[, country := substr(bvdid, 1, 2)]            # structural (between) attribute

  # one row per firm-year; basic validity filters
  fin <- fin[year >= YEAR_MIN & year <= YEAR_MAX &
               is.finite(total_assets) & total_assets > 0 &
               is.finite(turnover) & turnover > 0]
  setorder(fin, bvdid, year, -closdate)
  fin <- unique(fin, by = c("bvdid", "year"))      # keep latest closing per firm-year

  # derived predictors (economic meaning)
  fin[, `:=`(
    ln_assets   = log(total_assets),
    cf_margin   = 100 * cashflow / turnover,       # cash flow / operating revenue (%)
    wc_gap      = collection_days - credit_days,   # working-capital squeeze (days)
    rev_per_emp = oprev_per_emp,                   # operating revenue per employee (th)
    icover      = interest_cover                   # interest coverage = EBIT / interest paid
  )]

  # attach firm age, sector, status (small lookup tables)
  legal <- fread(file.path(SRC, "orbis_eu_legal.tsv"), sep = "\t", quote = "",
                 na.strings = c("", "n.a."), showProgress = FALSE)
  legal[, incorp_year := as.integer(substr(incorp_date, 1, 4))]
  legal <- unique(legal, by = "bvdid")[, .(bvdid, status, status_date, legal_form,
                                           incorp_year, listed)]
  nace <- fread(file.path(SRC, "orbis_eu_nace.tsv"), sep = "\t", quote = "",
                na.strings = c("", "n.a."), showProgress = FALSE)
  nace <- unique(nace, by = "bvdid")[, .(bvdid, nace_section, nace_core, bvd_sector)]
  fin <- merge(fin, legal, by = "bvdid", all.x = TRUE)
  fin <- merge(fin, nace,  by = "bvdid", all.x = TRUE)
  fin[, firm_age := pmax(year - incorp_year, 0L)]

  # winsorise heavy tails, then complete-case + balanced (>= MIN_YEARS) panel
  wins_vars <- c("icover", "gearing", "ebitda_margin", "current_ratio",
                 "liquidity_ratio", "cf_margin", "wc_gap", "stock_turnover",
                 "roa_pbt", "rev_per_emp", "solvency_asset")
  for (v in intersect(wins_vars, names(fin))) fin[[v]] <- winsor(fin[[v]])

  ok <- rowSums(is.finite(as.matrix(fin[, ..model_vars]))) == length(model_vars)
  pool <- fin[ok]
  pool[, nyr := .N, by = bvdid]
  pool <- pool[nyr >= MIN_YEARS]
  pool <- pool[, if (sd(solvency_asset) > 1e-8) .SD, by = bvdid]  # within variation in target
  saveRDS(as.data.frame(pool), pool_rds)           # provenance / robustness intermediate
} else {
  # FAST build: no raw extracts, restart from the saved pool.
  if (!file.exists(pool_rds))
    stop("Neither raw extracts nor orbis_eu_pool.rds found; run code/00_extract_orbis.sh first.")
  pool <- as.data.table(readRDS(pool_rds))
}

# ---- 7. reproducible firm sample (keeps e2tree O(N^2) tractable) -----------
all_firms <- unique(pool$bvdid)
sel_firms <- if (length(all_firms) > N_FIRMS) sample(all_firms, N_FIRMS) else all_firms
panel <- as.data.frame(pool[bvdid %in% sel_firms])
panel$year_c <- panel$year - round(mean(c(YEAR_MIN, YEAR_MAX)))

# ---- 8. Mundlak between/within decomposition ------------------------------
hpanel <- panel %>%
  group_by(bvdid) %>%
  mutate(across(all_of(model_vars),
                list(bw = ~mean(.x, na.rm = TRUE),
                     wn = ~.x - mean(.x, na.rm = TRUE)),
                .names = "{.col}_{.fn}")) %>%
  ungroup()
saveRDS(hpanel, file.path(DATA, "harmonized_panel.rds"))

# ---- 9. per-firm cross-section (long-run means + structural labels) --------
cross <- panel %>%
  group_by(bvdid, country, nace_section, bvd_sector, legal_form) %>%
  summarise(across(all_of(model_vars), ~mean(.x, na.rm = TRUE)),
            n_years = n(), firm_age = max(firm_age),
            .groups = "drop")
saveRDS(cross, file.path(DATA, "harmonized_cross.rds"))

# ---- 10. the complete-case modelling panel (main input) -------------------
md <- panel %>%
  transmute(bvdid, country, sector = nace_section, year,
            solvency_asset, roa_pbt, ebitda_margin, cf_margin, icover,
            current_ratio, liquidity_ratio, wc_gap, stock_turnover,
            ln_assets, firm_age, rev_per_emp)
saveRDS(md, file.path(DATA, "model_data.rds"))

# ---- 11. ICC (between share of variance) for target + features ------------
# ICC = between share of the observed variance, sigma_B^2 / (sigma_B^2 + sigma_W^2)
# with sigma_B^2 + sigma_W^2 = total variance. Single convention used throughout the
# paper: the same estimator as manuscript icc_of() in 03_panel_e2tree_orbis.R and
# supplement_tables.R, and the additive decomposition the Proposition-1 bound relies
# on. Computed directly (no aov over hundreds of thousands of factor levels). The full
# pool value is the structural population fact; the sampled panel is reported alongside.
fast_icc <- function(x, id) {
  dt <- data.table(x = as.numeric(x), id = id)
  dt <- dt[is.finite(x)]
  g  <- dt[, .(ni = .N, mi = mean(x), ssw = sum((x - mean(x))^2)), by = id]
  N  <- sum(g$ni); k <- nrow(g); gm <- sum(g$ni * g$mi) / N
  if (k < 2 || N <= k) return(NA_real_)
  vb <- sum(g$ni * (g$mi - gm)^2) / N     # between component (share of total variance)
  vw <- sum(g$ssw) / N                    # within component
  vb / (vb + vw)
}
icc_tbl <- data.frame(
  variable = model_vars,
  ICC_pool   = round(vapply(model_vars,
                            function(v) fast_icc(pool[[v]], pool$bvdid), numeric(1)), 3),
  ICC_sample = round(vapply(model_vars,
                            function(v) fast_icc(md[[v]], md$bvdid), numeric(1)), 3),
  row.names = NULL
)
write.csv(icc_tbl, file.path(OUT, "icc_by_var.csv"), row.names = FALSE)

# ---- 12. coverage + descriptives ------------------------------------------
cov_tbl <- data.frame(
  variable = model_vars,
  n_obs    = vapply(model_vars, function(v) sum(!is.na(md[[v]])), integer(1)),
  pct_obs  = round(100 * vapply(model_vars, function(v) mean(!is.na(md[[v]])), numeric(1)), 1),
  n_firms  = vapply(model_vars, function(v) length(unique(md$bvdid[!is.na(md[[v]])])), integer(1)),
  row.names = NULL
)
write.csv(cov_tbl, file.path(OUT, "coverage_by_var.csv"), row.names = FALSE)

desc <- md %>% select(all_of(model_vars)) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "x") %>%
  group_by(variable) %>%
  summarise(mean = round(mean(x, na.rm = TRUE), 2), sd = round(sd(x, na.rm = TRUE), 2),
            p05 = round(quantile(x, .05, na.rm = TRUE), 2),
            median = round(median(x, na.rm = TRUE), 2),
            p95 = round(quantile(x, .95, na.rm = TRUE), 2), .groups = "drop")
write.csv(desc, file.path(OUT, "descriptives.csv"), row.names = FALSE)

# ---- 13. diagnostic plots -------------------------------------------------
th <- theme_minimal(base_size = 12)
p_icc <- icc_tbl %>% mutate(variable = reorder(variable, ICC_pool)) %>%
  ggplot(aes(variable, ICC_pool)) +
  geom_col(fill = "#1b9e77", width = .7) +
  geom_hline(yintercept = 0.8, linetype = 2, color = "#d95f02") +
  coord_flip() + ylim(0, 1) + th +
  labs(x = NULL, y = "ICC (between share of variance, full pool)")
ggsave(file.path(PLT, "data_icc.png"), p_icc, width = 7, height = 5, dpi = 130)

p_n <- md %>% count(year) %>%
  ggplot(aes(year, n)) + geom_col(fill = "#2c7fb8", width = .7) + th +
  labs(x = NULL, y = "n. firms")
ggsave(file.path(PLT, "data_coverage.png"), p_n, width = 7, height = 4, dpi = 130)

# ---- 14. firms per country ------------------------------------------------
cty_tbl <- md %>% distinct(bvdid, country) %>% count(country, name = "n_firms") %>%
  arrange(desc(n_firms))
write.csv(cty_tbl, file.path(OUT, "firms_by_country.csv"), row.names = FALSE)
