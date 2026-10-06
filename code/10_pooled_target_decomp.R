# 10_pooled_target_decomp.R — Proposition 1 applied to the pooled-target variant.
#
# The paper reports one aggregate score for the pooled-target variant
# (R2 = 0.582, fidelity of the additive reconstruction to the predictions of the
# fitted pooled random forest) and then observes that an aggregate score cannot
# establish fidelity to the pooled model's within-firm variation. This script
# closes that gap: it recovers, on the same evaluation sample,
#
#   z_it = fhat(x_it)          the original pooled ensemble's predictions
#   g_it = b(xbar_i) + w(xtilde_it)   the additive surrogate reconstruction
#
# and applies the Proposition 1 definitions to the pair (z, g):
#
#   zbar_i, gbar_i          unit means of both series, on the evaluation sample
#   ztilde, gtilde          the respective within-unit deviations
#   eta_B = S_B/(S_B+S_W)   between share of the explanandum's variability
#   R2_total   = 1 - SSE/(S_B+S_W)
#   R2_between = 1 - sum_i T_i (zbar_i - gbar_i)^2 / S_B
#   R2_within  = 1 - sum (ztilde - gtilde)^2 / S_W
#
# and verifies the identity R2_total = eta_B R2_b + (1 - eta_B) R2_w.
#
# r2_within_of()/r2_between_of() are the helpers of 03_panel_e2tree_orbis.R:
# they are the Proposition 1 definitions (the observation-level form carries the
# T_i weights implicitly), so the numbers here are produced by the same code
# that produced decomp_r2_within / decomp_r2_between against the outcome.
#
# Three (z, g) pairs are reported:
#   (a) z = pooled predictions, g = pooled-target reconstruction   <- the paragraph
#   (b) z = pooled predictions, g = outcome-target reconstruction  <- the bridge, 0.563
#   (c) z = observed outcome,   g = pooled-target reconstruction   <- for completeness
#
# Nothing under output/repro/ is touched. Results go to output/pooled_target/.
# Usage: Rscript code/10_pooled_target_decomp.R          (POOLED_TARGET_FORCE=1 to refit)

suppressPackageStartupMessages({ library(dplyr); library(e2tree) })

BASE <- normalizePath(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))), ".."))
stopifnot(file.exists(file.path(BASE, "data", "model_data.rds")))

OUT <- file.path(BASE, "output", "pooled_target")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
FIT <- file.path(OUT, "pooled_target_fit.rds")
FORCE <- identical(Sys.getenv("POOLED_TARGET_FORCE"), "1")

source(file.path(BASE, "code", "panel_e2tree.R"))

log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

# --- helpers copied verbatim from 03_panel_e2tree_orbis.R ------------------
icc_of <- function(x, g) {
  m <- tapply(x, g, mean); ng <- tapply(x, g, length)
  vb <- sum(ng * (m - mean(x))^2) / length(x)
  vw <- sum((x - m[match(g, names(m))])^2) / length(x)
  vb / (vb + vw)
}
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

# --- Proposition 1, spelled out and self-checked ---------------------------
prop1 <- function(z, g, unit, label) {
  stopifnot(length(z) == length(g), length(z) == length(unit))
  zbar <- stats::ave(z, unit); gbar <- stats::ave(g, unit)
  ztil <- z - zbar;            gtil <- g - gbar
  S_B  <- sum((zbar - mean(z))^2)        # = sum_i T_i (zbar_i - zbar)^2
  S_W  <- sum(ztil^2)
  eta  <- S_B / (S_B + S_W)
  r2_t <- 1 - sum((z - g)^2) / (S_B + S_W)
  r2_b <- 1 - sum((zbar - gbar)^2) / S_B
  r2_w <- 1 - sum((ztil - gtil)^2) / S_W
  data.frame(pair = label,
             eta_B = eta, icc_z = icc_of(z, unit),
             r2_total = r2_t, r2_between = r2_b, r2_within = r2_w,
             r2_total_rebuilt = eta * r2_b + (1 - eta) * r2_w,
             r2_between_helper = r2_between_of(z, g, unit),
             r2_within_helper  = r2_within_of(z, g, unit),
             sd_z = stats::sd(z), sd_g = stats::sd(g),
             sd_z_within = stats::sd(ztil), sd_g_within = stats::sd(gtil),
             stringsAsFactors = FALSE)
}

# --- data, identical to S0 of 03_panel_e2tree_orbis.R ----------------------
OUTCOME <- "solvency_asset"; UNIT <- "bvdid"; TIME <- "year"
md    <- readRDS(file.path(BASE, "data", "model_data.rds"))
feats <- c("roa_pbt","ebitda_margin","cf_margin","icover","current_ratio",
           "liquidity_ratio","wc_gap","stock_turnover","ln_assets","firm_age","rev_per_emp")
stopifnot(all(c(UNIT, TIME, OUTCOME, feats) %in% names(md)))
md <- md[complete.cases(md[, c(UNIT, TIME, OUTCOME, feats)]), ]
form_solv <- reformulate(feats, OUTCOME)
SET_B <- list(impTotal = 0.10, maxDec = 1e-6, n = 2, level = 5)
SET_W <- list(impTotal = 0.05, maxDec = 1e-7, n = 5, level = 6)
log_msg("data: ", nrow(md), " obs, ", length(unique(md[[UNIT]])), " firms")

# --- S5 of the published pipeline, refitted (predictions were not cached) ---
if (!FORCE && file.exists(FIT)) {
  log_msg("loading cached pooled-target fit from ", basename(FIT))
  keep <- readRDS(FIT)
} else {
  log_msg("refitting the pooled-target variant (slow step: two dissimilarity matrices)")
  t0 <- Sys.time()
  set.seed(123)
  pooled_ens <- ranger::ranger(form_solv, data = md, num.trees = 500)
  mp <- panel_e2tree(form_solv, data = md, unit = UNIT, time = TIME,
                     engine = "ranger", ntree = 500,
                     setting_between = SET_B, setting_within = SET_W,
                     target = "pooled", pooled_ensemble = pooled_ens, seed = 123)
  log_msg("fit done in ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
  keep <- list(                              # everything needed downstream, no D
    predictions = mp$predictions,
    r2_panel = mp$r2_panel, cor2 = mp$outcome_var_recovered,
    vs_outcome = mp$vs_outcome,
    between_fid = mp$between$fidelity, within_fid = mp$within$fidelity,
    between_leaves = n_leaves(mp$between$tree),
    within_leaves  = n_leaves(mp$within$tree),
    between_vars = mp$between$variables, within_vars = mp$within$variables,
    between_oob = mp$between$ensemble$r.squared,
    within_oob  = mp$within$ensemble$r.squared,
    unit = mp$predictions$unit)
  saveRDS(keep, FIT)
  log_msg("fit cached to ", FIT)
}

pr <- keep$predictions
stopifnot(nrow(pr) == nrow(md), !anyDuplicated(pr$.row),
          identical(as.character(pr$unit), as.character(md[[UNIT]])[pr$.row]))

# --- gate: the refit must reproduce the published S5 scalars ----------------
pub <- c(pooledtarget_r2 = 0.582, pooledtarget_within_fid = 0.666,
         between_fid = 0.722, icc_pooled_preds = 0.843)
got <- c(pooledtarget_r2 = round(keep$r2_panel, 3),
         pooledtarget_within_fid = round(keep$within_fid, 3),
         between_fid = round(keep$between_fid, 3),
         icc_pooled_preds = round(icc_of(pr$.pooled, pr$unit), 3))
print(rbind(published = pub, refit = got))
stopifnot(isTRUE(all.equal(pub, got, tolerance = 0)))
log_msg("gate passed: the refit reproduces the published pooled-target scalars")

# --- the three Proposition 1 decompositions --------------------------------
m_main <- readRDS(file.path(BASE, "output", "repro", "cache", "S1_main.rds"))
stopifnot(nrow(m_main$predictions) == nrow(pr))
g_out  <- m_main$predictions$.panel[match(pr$.row, m_main$predictions$.row)]
stopifnot(!anyNA(g_out))

res <- rbind(
  prop1(pr$.pooled, pr$.panel, pr$unit, "fhat vs pooled-target reconstruction"),
  prop1(pr$.pooled, g_out,     pr$unit, "fhat vs outcome-target reconstruction"),
  prop1(pr$outcome, pr$.panel, pr$unit, "outcome vs pooled-target reconstruction"))

# identity check and consistency with the pipeline helpers
stopifnot(max(abs(res$r2_total - res$r2_total_rebuilt)) < 1e-10,
          max(abs(res$r2_between - res$r2_between_helper)) < 1e-10,
          max(abs(res$r2_within  - res$r2_within_helper))  < 1e-10)
# the aggregate score of the first row is the published 0.582
stopifnot(abs(res$r2_total[1] - keep$r2_panel) < 1e-10,
          abs(res$r2_total[2] - r2_honest(pr$.pooled, g_out)) < 1e-10)
# eta_B of the explanandum is its ICC: the paper's \iccPooledPreds is the weight
stopifnot(abs(res$eta_B[1] - res$icc_z[1]) < 1e-12,
          round(res$eta_B[1], 3) == pub[["icc_pooled_preds"]])
log_msg("Proposition 1 identity holds to 1e-10 on all three pairs; eta_B = ICC confirmed")

# how much of the pooled model's within variation the reconstruction even moves
within_scale <- data.frame(
  sd_fhat_within   = stats::sd(pr$.pooled - stats::ave(pr$.pooled, pr$unit)),
  sd_recon_within  = stats::sd(pr$.panel  - stats::ave(pr$.panel,  pr$unit)),
  cor_within       = stats::cor(pr$.pooled - stats::ave(pr$.pooled, pr$unit),
                                pr$.panel  - stats::ave(pr$.panel,  pr$unit)),
  cor2_within      = stats::cor(pr$.pooled - stats::ave(pr$.pooled, pr$unit),
                                pr$.panel  - stats::ave(pr$.panel,  pr$unit))^2,
  cor_between      = stats::cor(stats::ave(pr$.pooled, pr$unit),
                                stats::ave(pr$.panel,  pr$unit)),
  cor2_between     = stats::cor(stats::ave(pr$.pooled, pr$unit),
                                stats::ave(pr$.panel,  pr$unit))^2)

print(res[, c("pair", "eta_B", "r2_total", "r2_between", "r2_within")], digits = 4)
print(within_scale, digits = 4)

write.csv(res, file.path(OUT, "pooled_target_prop1.csv"), row.names = FALSE)
write.csv(within_scale, file.path(OUT, "pooled_target_within_scale.csv"), row.names = FALSE)

# --- macros for the manuscript ---------------------------------------------
# Only the macros this run adds. \pooledTargetR, \bridgeR, \iccPooledPreds and
# \pooledTargetWithinFid are already defined in numbers.tex by the v1 pipeline;
# re-emitting them here would be a duplicate \newcommand (a LaTeX error) and a
# second place for the same number to drift. The gate above is what ties the
# two files together: it fails if this refit stops reproducing them.
fmt3 <- function(x) sprintf("%.3f", x)
fmt2 <- function(x) sprintf("%.2f", x)
nm <- c(
  pooledTargetEtaB       = fmt3(res$eta_B[1]),
  pooledTargetWithinW    = fmt3(1 - res$eta_B[1]),
  pooledTargetBetween    = fmt3(res$r2_between[1]),
  pooledTargetWithin     = fmt3(res$r2_within[1]),
  bridgeBetween          = fmt3(res$r2_between[2]),
  bridgeWithin           = fmt3(res$r2_within[2]),
  pooledTargetVsYR       = fmt3(res$r2_total[3]),
  pooledTargetVsYBetween = fmt3(res$r2_between[3]),
  pooledTargetVsYWithin  = fmt3(res$r2_within[3]),
  pooledModelWithinSd    = fmt2(within_scale$sd_fhat_within),
  pooledTargetWithinSd   = fmt2(within_scale$sd_recon_within),
  pooledTargetWithinCor  = fmt2(within_scale$cor_within),
  pooledTargetBetweenFid = fmt3(keep$between_fid),
  pooledTargetBetweenLeaves = as.character(keep$between_leaves),
  pooledTargetWithinLeaves  = as.character(keep$within_leaves))
TEX <- file.path(BASE, "manuscript", "manuscript_DSS_v2", "numbers_pooled_target.tex")
writeLines(c(
  "% AUTO-GENERATED by code/10_pooled_target_decomp.R from output/pooled_target/*.csv",
  "% Do not edit by hand: re-run code/10_pooled_target_decomp.R.",
  "% Holds ONLY the Proposition 1 decomposition of the pooled-target variant;",
  "% numbers.tex and numbers_budget.tex hold the rest. All three are \\input",
  "% by manuscript.tex: none of them is a newer version of another.",
  sprintf("\\newcommand{\\%s}{%s}", names(nm), nm)), TEX)
write.csv(data.frame(key = names(nm), value = unname(nm)),
          file.path(OUT, "pooled_target_numbers.csv"), row.names = FALSE)
cat("written:", TEX, "\n"); writeLines(sprintf("  \\%s = %s", names(nm), nm))
log_msg("written to ", OUT)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
