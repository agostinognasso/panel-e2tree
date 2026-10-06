# 06_within_no_age.R — robustness run prompted by the v2 verification, not by CI:
# firm_age (demeaned) is numerically identical to centred calendar time, so the
# within tree has a pure time channel among its split variables. This refits the
# whole protocol with firm_age removed, to see whether the within fidelity of
# 0.477 survives without it.
#
# S10 in 03_panel_e2tree_orbis.R already drops three variables at once
# (ln_assets, firm_age, rev_per_emp) and gives 0.480. This run isolates firm_age.
#
# Nothing under output/repro/ is touched. Results go to output/budget59/.
# Usage: Rscript code/06_within_no_age.R

suppressPackageStartupMessages({ library(dplyr); library(e2tree) })

BASE <- normalizePath(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))), ".."))
source(file.path(BASE, "code", "panel_e2tree.R"))
source(file.path(BASE, "code", "e2tree_split_local.R"))

OUT <- file.path(BASE, "output", "budget59")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

r2_honest <- function(obs, pred) 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)
n_leaves  <- function(tr) sum(tr$tree$terminal)
r2_within_of <- function(y, pred, g) {
  yd <- y - stats::ave(y, g); pd <- pred - stats::ave(pred, g)
  1 - sum((yd - pd)^2) / sum(yd^2)
}
r2_between_of <- function(y, pred, g) {
  yb <- stats::ave(y, g); pb <- stats::ave(pred, g)
  1 - sum((yb - pb)^2) / sum((yb - mean(yb))^2)
}

OUTCOME <- "solvency_asset"; UNIT <- "bvdid"; TIME <- "year"
md    <- readRDS(file.path(BASE, "data", "model_data.rds"))
feats <- c("roa_pbt","ebitda_margin","cf_margin","icover","current_ratio",
           "liquidity_ratio","wc_gap","stock_turnover","ln_assets","firm_age","rev_per_emp")
md <- md[complete.cases(md[, c(UNIT, TIME, OUTCOME, feats)]), ]

SET_B <- list(impTotal = 0.10, maxDec = 1e-6, n = 2, level = 5)
SET_W <- list(impTotal = 0.05, maxDec = 1e-7, n = 5, level = 6)

# --- the demeaned-firm_age == centred-time identity, restated for the record ---
dev_age <- md$firm_age - stats::ave(md$firm_age, md[[UNIT]])
dev_yr  <- md[[TIME]]  - stats::ave(md[[TIME]],  md[[UNIT]])
log_msg("max |demeaned firm_age - centred year| = ", max(abs(dev_age - dev_yr)),
        "  (distinct values: ", length(unique(round(dev_age, 10))), ")")

runs <- list(
  full     = feats,
  no_age   = setdiff(feats, "firm_age"),
  no_age_size = setdiff(feats, c("firm_age", "ln_assets"))   # the SHAP 9-predictor set
)

rows <- list()
for (nm in names(runs)) {
  f <- runs[[nm]]
  log_msg("fitting '", nm, "' with ", length(f), " features ...")
  t0 <- Sys.time()
  m <- suppressWarnings(
    panel_e2tree(reformulate(f, OUTCOME), data = md, unit = UNIT, time = TIME,
                 engine = "ranger", ntree = 500,
                 setting_between = SET_B, setting_within = SET_W, seed = 123))
  p <- m$predictions$.panel
  rows[[nm]] <- data.frame(
    run = nm, n_features = length(f),
    between_fidelity = m$between$fidelity, within_fidelity = m$within$fidelity,
    panel_r2 = m$r2_panel, within_signal = m$within_signal,
    between_leaves = n_leaves(m$between$tree),
    within_leaves  = if (is.null(m$within$tree)) NA_integer_ else n_leaves(m$within$tree),
    r2_within  = r2_within_of(md[[OUTCOME]], p, md[[UNIT]]),
    r2_between = r2_between_of(md[[OUTCOME]], p, md[[UNIT]]),
    within_vars = paste(m$within$variables, collapse = ", "),
    minutes = as.numeric(round(difftime(Sys.time(), t0, units = "mins"), 2)),
    stringsAsFactors = FALSE)
  log_msg(sprintf("  %-12s within fidelity %.3f (%d leaves), panel R2 %.3f, r2_within %+.3f",
                  nm, m$within$fidelity, rows[[nm]]$within_leaves, m$r2_panel, rows[[nm]]$r2_within))
  write.csv(do.call(rbind, rows), file.path(OUT, "within_no_age.csv"), row.names = FALSE)
}

log_msg("=== done ===")
print(do.call(rbind, rows)[, setdiff(names(rows[[1]]), "within_vars")], digits = 4)
for (nm in names(rows)) cat("\n", nm, " within vars: ", rows[[nm]]$within_vars, "\n", sep = "")
