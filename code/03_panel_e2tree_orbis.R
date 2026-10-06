# 03_panel_e2tree_orbis.R — single source for every number of the ORBIS
# application of the "Panel e2tree" paper, on the European firm panel:
#   outcome = solvency_asset (equity/total assets, %)   unit = bvdid   time = year
#
# Usage:   Rscript code/03_panel_e2tree_orbis.R     (from the project root)
# Output:  output/repro/*.csv  (one file per table/figure) + paper_numbers.csv
#          manuscript/manuscript_DSS_v2/figures/{e2tree_necessity,proximity_decomp}.png
#          output/repro/sessionInfo.txt
# Each step caches in output/repro/cache/ and is skipped on re-run
# (delete the cache or set REPRO_FORCE=1 to recompute).

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
BASE <- if (length(script_path)) dirname(dirname(normalizePath(script_path))) else normalizePath("..")
stopifnot(file.exists(file.path(BASE, "data", "model_data.rds")))

OUT   <- file.path(BASE, "output", "repro")
CACHE <- file.path(OUT, "cache")
FEN   <- file.path(BASE, "manuscript", "manuscript_DSS_v2", "figures")
# e2tree_necessity.png is written here and then overwritten by 07_figures_v2.R,
# which adds the 59-leaf budget row. Keep the order 03 -> 07.
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)
dir.create(FEN,   recursive = TRUE, showWarnings = FALSE)

FORCE <- identical(Sys.getenv("REPRO_FORCE"), "1")

suppressPackageStartupMessages(library(e2tree))
source(file.path(BASE, "code", "panel_e2tree.R"))
source(file.path(BASE, "code", "e2tree_split_local.R"))

step <- function(name, expr) {
  f <- file.path(CACHE, paste0(name, ".rds"))
  if (!FORCE && file.exists(f)) return(readRDS(f))
  val <- expr; saveRDS(val, f); val
}
numbers <- list()
# put(key, value, digits): `digits` rounds the value and records the number of
# decimals the LaTeX macro must print, so 0.58 is written "0.580" and the text
# never shows a result with fewer decimals than the others. The value stays
# numeric: the tables and figures at the bottom still compute with it.
put <- function(key, value, digits = NULL) {
  if (!is.null(digits)) { value <- round(value, digits); attr(value, "digits") <- digits }
  numbers[[key]] <<- value
}
fmt_number <- function(v) {
  d <- attr(v, "digits")
  if (is.null(d)) paste(v, collapse = " ") else paste(sprintf("%.*f", d, v), collapse = " ")
}

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

# --- Proposition 1: regularity constant c_B and threshold ICC*, operationalised.
# The proposition is stated for a variance-reduction (CART) surrogate; c_B is
# calibrated by the leading between contrast -- the best between split at the root,
# where all N units are still mixed -- on the unit-mean data {(xbar_i, ybar_i)}:
#   c_B = K * Delta_B* / sigma_B^2 ,   ICC* = sigma_B^2 / (sigma_B^2 + Delta_B*).
# var_parts() returns the ANOVA between/within variance components (sigma_B^2, sigma_W^2);
# leading_between_gain() returns Delta_B*/sigma_B^2, the share of between SS the root split cuts.
var_parts <- function(x, g) {
  m <- tapply(x, g, mean); ng <- tapply(x, g, length)
  c(vb = sum(ng * (m - mean(x))^2) / length(x),
    vw = sum((x - m[match(g, names(m))])^2) / length(x))
}
leading_between_gain <- function(um, ycol, xcols) {
  f  <- reformulate(intersect(xcols, names(um)), ycol)
  ct <- rpart::rpart.control(minsplit = 2, minbucket = 1, cp = 0,
                             maxcompete = 0, maxsurrogate = 0, xval = 0)
  fr <- rpart::rpart(f, data = um, method = "anova", control = ct)$frame
  dev <- stats::setNames(fr$dev, rownames(fr))
  as.numeric((dev[["1"]] - dev[["2"]] - dev[["3"]]) / dev[["1"]])   # root split, node 1 -> {2,3}
}

SET_B <- list(impTotal = 0.10, maxDec = 1e-6, n = 2, level = 5)
SET_W <- list(impTotal = 0.05, maxDec = 1e-7, n = 5, level = 6)

fit_single_e2 <- function(d, ylab, setting, seed = 7) {
  set.seed(seed)
  form <- reformulate(setdiff(names(d), ylab), ylab)
  ens  <- ranger::ranger(form, data = d, num.trees = 500)
  D    <- createDisMatrix(ens, data = d, label = ylab,
                          parallel = list(active = FALSE, no_cores = 1))
  tr   <- e2tree(form, d, D, ens, setting)
  fit  <- as.numeric(predict(tr, newdata = d)$fit)
  ep   <- get_ensemble_predictions(ens, d, type = "regression")
  list(tree = tr, ensemble = ens,
       fid_vs_ens = cor(fit, ep)^2, cor2_vs_y = cor(fit, d[[ylab]])^2,
       r2_vs_y = r2_honest(d[[ylab]], fit), leaves = n_leaves(tr),
       vars = unique(na.omit(tr$tree$variable)), fit = fit)
}

# ---------------------------------------------------------------------------
# S0 — data and descriptives
# ---------------------------------------------------------------------------
OUTCOME <- "solvency_asset"; UNIT <- "bvdid"; TIME <- "year"
md    <- readRDS(file.path(BASE, "data", "model_data.rds"))
feats <- c("roa_pbt","ebitda_margin","cf_margin","icover","current_ratio",
           "liquidity_ratio","wc_gap","stock_turnover","ln_assets","firm_age","rev_per_emp")
stopifnot(all(c(UNIT, TIME, OUTCOME, feats) %in% names(md)))
md <- md[complete.cases(md[, c(UNIT, TIME, OUTCOME, feats)]), ]

put("n_obs",       nrow(md))
put("n_firms",     length(unique(md[[UNIT]])))
put("n_countries", length(unique(md$country)))
put("year_min",    min(md[[TIME]])); put("year_max", max(md[[TIME]]))
put("icc_solv",    icc_of(md[[OUTCOME]], md[[UNIT]]), 3)
form_solv <- reformulate(feats, OUTCOME)

# population (full-pool) ICC of the target, the structural necessity fact
poolf <- file.path(BASE, "data", "orbis_eu_pool.rds")
if (file.exists(poolf)) {
  pl <- readRDS(poolf)
  put("icc_solv_poolwide", icc_of(pl[[OUTCOME]], pl[[UNIT]]), 3)
  put("n_firms_pool", format(length(unique(pl[[UNIT]])), big.mark = "{,}"))
  rm(pl)
}

# ---------------------------------------------------------------------------
# S1 — main result: panel_e2tree (target = outcome, within = unit)
# ---------------------------------------------------------------------------
m_main <- step("S1_main", {
  panel_e2tree(form_solv, data = md, unit = UNIT, time = TIME,
               engine = "ranger", ntree = 500,
               setting_between = SET_B, setting_within = SET_W, seed = 123)
})
put("between_fidelity", m_main$between$fidelity, 3)
put("within_fidelity",  m_main$within$fidelity, 3)
put("panel_r2",         m_main$r2_panel, 3)
put("panel_cor2",       m_main$outcome_var_recovered, 3)
put("within_signal",    m_main$within_signal, 3)
put("between_leaves",   n_leaves(m_main$between$tree))
put("within_leaves",    n_leaves(m_main$within$tree))
put("between_vars",     paste(m_main$between$variables, collapse = ", "))
put("within_vars",      paste(m_main$within$variables, collapse = ", "))
put("within_oob_r2",    m_main$within$ensemble$r.squared, 3)
put("gap_bound",        (1 - numbers[["icc_solv"]]) * numbers[["within_oob_r2"]], 3)
put("decomp_r2_within",  r2_within_of(md[[OUTCOME]], m_main$predictions$.panel, md[[UNIT]]), 3)
put("decomp_r2_between", r2_between_of(md[[OUTCOME]], m_main$predictions$.panel, md[[UNIT]]), 3)

# tree-reading anchors quoted in the walk-through of Section "Application":
# between root threshold, lowest between leaf (%), deepest within leaf (|points|)
tb_fr <- m_main$between$tree$tree
tw_fr <- m_main$within$tree$tree
put("between_root_thr",   as.numeric(sub("^.*<=\\s*", "", tb_fr$splitLabel[1])), 1)
put("between_low_leaf",   round(min(as.numeric(tb_fr$pred[tb_fr$terminal]))))
put("within_low_leaf_abs", round(abs(min(as.numeric(tw_fr$pred[tw_fr$terminal])))))

# Proposition 1 made numerical: leading-contrast constant c_B and threshold ICC*
vp_solv     <- var_parts(md[[OUTCOME]], md[[UNIT]])
R2root_solv <- leading_between_gain(m_main$unit_means, OUTCOME, feats)
Kbud_solv   <- n_leaves(m_main$between$tree)
put("sigma_w2_solv", round(vp_solv[["vw"]]))
put("delta_b_solv",  round(R2root_solv * vp_solv[["vb"]]))
put("c_b_solv",      Kbud_solv * R2root_solv, 1)
put("icc_star_solv", 1 / (1 + R2root_solv), 3)

# ---------------------------------------------------------------------------
# S2 — shortcut (a): pooled single e2tree on raw features (+ capacity-matched)
# ---------------------------------------------------------------------------
slim <- function(f, y, g) {
  out <- f[c("fid_vs_ens", "cor2_vs_y", "r2_vs_y", "leaves", "vars")]
  out$r2_within  <- r2_within_of(y, f$fit, g)
  out$r2_between <- r2_between_of(y, f$fit, g)
  out
}
s2 <- step("S2_pooled", {
  pooled  <- fit_single_e2(md[, c(OUTCOME, feats)], OUTCOME, SET_W, seed = 7)
  matched <- fit_single_e2(md[, c(OUTCOME, feats)], OUTCOME,
                           list(impTotal = 0.02, maxDec = 1e-8, n = 2, level = 7), seed = 7)
  list(pooled = slim(pooled, md[[OUTCOME]], md[[UNIT]]),
       matched = slim(matched, md[[OUTCOME]], md[[UNIT]]))
})
put("pooled_r2_vs_y",   s2$pooled$r2_vs_y, 3)
put("pooled_cor2_vs_y", s2$pooled$cor2_vs_y, 3)
put("pooled_fid_vs_ens",s2$pooled$fid_vs_ens, 3)
put("pooled_leaves",    s2$pooled$leaves)
put("pooled_r2_within", s2$pooled$r2_within, 3)
put("pooled_r2_between",s2$pooled$r2_between, 3)
put("matched_r2_vs_y",  s2$matched$r2_vs_y, 3)
put("matched_leaves",   s2$matched$leaves)
put("matched_r2_within",s2$matched$r2_within, 3)
put("matched_r2_between",s2$matched$r2_between, 3)
put("decomp_leaves_total", numbers[["between_leaves"]] + numbers[["within_leaves"]])

# ---------------------------------------------------------------------------
# S3 — shortcut (c): single e2tree on Mundlak-augmented features + bridge evidence
# ---------------------------------------------------------------------------
s3 <- step("S3_augmented", {
  aug <- md %>%
    group_by(.data[[UNIT]]) %>%
    mutate(across(all_of(feats), list(bw = ~mean(.x), wn = ~.x - mean(.x)))) %>%
    ungroup() %>%
    select(all_of(OUTCOME), ends_with("_bw"), ends_with("_wn")) %>%
    as.data.frame()
  keep <- vapply(aug, function(z) sd(z) > 1e-8, logical(1)) | names(aug) == OUTCOME
  aug  <- aug[, keep, drop = FALSE]
  set.seed(7)
  form <- reformulate(setdiff(names(aug), OUTCOME), OUTCOME)
  ens  <- ranger::ranger(form, data = aug, num.trees = 500)
  D    <- createDisMatrix(ens, data = aug, label = OUTCOME,
                          parallel = list(active = FALSE, no_cores = 1))
  tr   <- e2tree(form, aug, D, ens, SET_W)
  fit  <- as.numeric(predict(tr, newdata = aug)$fit)
  ep   <- get_ensemble_predictions(ens, aug, type = "regression")
  X    <- aug[, setdiff(names(aug), OUTCOME), drop = FALSE]
  S    <- e2_candidate_splits(X, max_cat = 10, max_thresholds = 256)$S
  n    <- nrow(aug)
  imp0 <- e2tree:::eImpurity(D, seq_len(n), S)
  dec0 <- sum(D) / (n * (n - 1)) - imp0
  top10 <- sub(" .*$", "", names(sort(dec0, decreasing = TRUE))[1:10])
  sv_root <- vapply(seq_len(200), function(b) {
    ti <- ranger::treeInfo(ens, b); ti$splitvarName[ti$nodeID == 0][1]
  }, character(1))
  sv_all <- unlist(lapply(seq_len(100), function(b) {
    ti <- ranger::treeInfo(ens, b)$splitvarName; ti[!is.na(ti)]
  }))
  list(fid_vs_ens = if (sd(fit) < 1e-12) NA_real_ else cor(fit, ep)^2,
       cor2_vs_y = if (sd(fit) < 1e-12) NA_real_ else cor(fit, aug[[OUTCOME]])^2,
       r2_vs_y = r2_honest(aug[[OUTCOME]], fit),
       r2_within = r2_within_of(md[[OUTCOME]], fit, md[[UNIT]]),
       r2_between = r2_between_of(md[[OUTCOME]], fit, md[[UNIT]]),
       vars = unique(na.omit(tr$tree$variable)), leaves = n_leaves(tr),
       top10_bw_share = mean(grepl("_bw$", top10)),
       ens_root_bw_share = mean(grepl("_bw$", sv_root), na.rm = TRUE),
       ens_share_bw = mean(grepl("_bw$", sv_all)))
})
put("aug_r2_vs_y",      s3$r2_vs_y, 3)
put("aug_r2_within",    s3$r2_within, 3)
put("aug_r2_between",   s3$r2_between, 3)
put("aug_splits_bw",    sum(grepl("_bw$", s3$vars)))
put("aug_splits_wn",    sum(grepl("_wn$", s3$vars)))
put("aug_leaves",       s3$leaves)
put("aug_top10_bw",     round(100 * s3$top10_bw_share))   # percentage of the top-10 root candidates
put("ens_root_bw_share",s3$ens_root_bw_share, 3)
put("ens_share_bw",     s3$ens_share_bw, 3)

# ---------------------------------------------------------------------------
# S4 — shortcut (b): PERMANOVA decomposition of a pooled proximity
# ---------------------------------------------------------------------------
permanova_share <- function(D, groups) {
  N <- nrow(D); D2 <- D^2
  SST <- sum(D2[upper.tri(D2)]) / N; SSW <- 0
  for (g in unique(groups)) {
    idx <- which(groups == g); ng <- length(idx)
    if (ng > 1) { sub <- D2[idx, idx]; SSW <- SSW + sum(sub[upper.tri(sub)]) / ng }
  }
  c(between = (SST - SSW) / SST, within = SSW / SST)
}
s4 <- step("S4_permanova", {
  set.seed(123)
  ens_pool <- randomForest::randomForest(form_solv, md, ntree = 400, proximity = TRUE)
  sh_real  <- permanova_share(1 - ens_pool$proximity, md[[UNIT]])
  sim_share <- function(icc_target) {              # negative control: proximity between-share vs true ICC
    nc <- 60; Tt <- 12
    A <- runif(nc, 0, 10); Bbar <- runif(nc, 0, 5)
    d <- expand.grid(i = seq_len(nc), t = seq_len(Tt))
    d$A <- A[d$i]; d$B <- Bbar[d$i] + rnorm(nrow(d), 0, 2)
    d <- d %>% group_by(i) %>% mutate(Bw = B - mean(B)) %>% ungroup()
    wB <- sqrt(icc_target); wW <- sqrt(1 - icc_target)
    d$y <- 50 + 6*wB*scale(d$A)[,1] + 6*wW*scale(d$Bw)[,1] + rnorm(nrow(d), 0, 0.8)
    e <- randomForest::randomForest(y ~ A + B + Bw, d, ntree = 400, proximity = TRUE)
    s <- permanova_share(1 - e$proximity, d$i)
    c(icc_target = icc_target, icc_emp = icc_of(d$y, d$i), prox_between = unname(s["between"]))
  }
  sims <- as.data.frame(do.call(rbind, lapply(c(.3, .5, .7, .85, .95), sim_share)))
  list(real_between = unname(sh_real["between"]), sims = sims)
})
write.csv(s4$sims, file.path(OUT, "proximity_decomp.csv"), row.names = FALSE)
put("prox_between_real",  s4$real_between, 3)
put("prox_between_range", paste(sprintf("%.2f", range(s4$sims$prox_between)), collapse = "-"))
put("prox_icc_range",     paste(sprintf("%.2f", range(s4$sims$icc_emp)), collapse = "-"))

# ---------------------------------------------------------------------------
# S5 — target = "pooled": explaining a given pooled ensemble
# ---------------------------------------------------------------------------
s5 <- step("S5_pooled_target", {
  set.seed(123)
  pooled_ens <- ranger::ranger(form_solv, data = md, num.trees = 500)
  mp <- panel_e2tree(form_solv, data = md, unit = UNIT, time = TIME,
                     engine = "ranger", ntree = 500,
                     setting_between = SET_B, setting_within = SET_W,
                     target = "pooled", pooled_ensemble = pooled_ens, seed = 123)
  fhat <- mp$predictions$.pooled
  list(r2_vs_pooled = mp$r2_panel, cor2_vs_pooled = mp$outcome_var_recovered,
       r2_vs_outcome = mp$vs_outcome$r2,
       between_fid = mp$between$fidelity, within_fid = mp$within$fidelity,
       bridge_r2 = r2_honest(fhat, m_main$predictions$.panel),
       icc_fhat = icc_of(fhat, md[[UNIT]]))
})
put("pooledtarget_r2",        s5$r2_vs_pooled, 3)
put("pooledtarget_within_fid",s5$within_fid, 3)
put("bridge_r2",              s5$bridge_r2, 3)
put("icc_pooled_preds",       s5$icc_fhat, 3)

# ---------------------------------------------------------------------------
# S6 — within = "twoway": common period effects separated (captures 2020 COVID)
# ---------------------------------------------------------------------------
s6 <- step("S6_twoway", {
  mt <- panel_e2tree(form_solv, data = md, unit = UNIT, time = TIME,
                     engine = "ranger", ntree = 500,
                     setting_between = SET_B, setting_within = SET_W,
                     within = "twoway", seed = 123)
  tm <- mt$time_means[order(mt$time_means$.time), c(".time", OUTCOME)]
  list(r2_panel = mt$r2_panel, within_fid = mt$within$fidelity,
       within_vars = mt$within$variables, time_effects = tm,
       trend_cor = cor(tm$.time, tm[[OUTCOME]]))
})
write.csv(s6$time_effects, file.path(OUT, "twoway_time_effects.csv"), row.names = FALSE)
put("twoway_r2",         s6$r2_panel, 3)
put("twoway_within_fid", s6$within_fid, 3)
put("twoway_trend_cor",  s6$trend_cor, 3)

# ---------------------------------------------------------------------------
# S7 — out-of-time validation: train <= 2020, explain 2021-2022
# ---------------------------------------------------------------------------
s7 <- step("S7_holdout", {
  tr_idx <- md[[TIME]] <= 2020
  m_tr <- panel_e2tree(form_solv, data = md[tr_idx, ], unit = UNIT, time = TIME,
                       engine = "ranger", ntree = 500,
                       setting_between = SET_B, setting_within = SET_W, seed = 123)
  te <- md[!tr_idx, ]
  pr <- suppressWarnings(predict(m_tr, newdata = te))
  list(n_test = nrow(te), oot_r2 = r2_honest(te[[OUTCOME]], pr$.panel),
       oot_cor2 = cor(pr$.panel, te[[OUTCOME]])^2,
       oot_r2_within  = r2_within_of (te[[OUTCOME]], pr$.panel, te[[UNIT]]),
       oot_r2_between = r2_between_of(te[[OUTCOME]], pr$.panel, te[[UNIT]]))
})
put("oot_n", s7$n_test); put("oot_r2", s7$oot_r2, 3); put("oot_cor2", s7$oot_cor2, 3)
put("oot_train_max", 2020); put("oot_test_min", 2021)  # the split years of tr_idx above
if (!is.null(s7$oot_r2_within)) {   # absent only in pre-extension caches
  put("oot_r2_within",  s7$oot_r2_within, 3)
  put("oot_r2_between", s7$oot_r2_between, 3)
}

# ---------------------------------------------------------------------------
# S8 — simulation with known roles (method validation; synthetic, generic)
# ---------------------------------------------------------------------------
s8 <- step("S8_simulation", {
  set.seed(123)
  nc <- 60; Tt <- 12
  sim <- expand.grid(i = seq_len(nc), t = seq_len(Tt))
  Ai <- runif(nc, 0, 10); Bbar <- runif(nc, 0, 5)
  sim$A <- Ai[sim$i]; sim$B <- Bbar[sim$i] + rnorm(nrow(sim), 0, 2)
  sim$Z1 <- rnorm(nrow(sim)); sim$Z2 <- runif(nrow(sim), 0, 10)
  sim <- sim %>% group_by(i) %>% mutate(Bw = B - mean(B)) %>% ungroup()
  sim$y <- 50 + 2.5 * sim$A + 3 * sim$Bw + rnorm(nrow(sim), 0, 1.5)
  icc_sim <- icc_of(sim$y, sim$i)
  S  <- as.data.frame(sim[, c("y", "A", "B", "Z1", "Z2")])
  ps <- fit_single_e2(S, "y", SET_W, seed = 7)
  ms <- panel_e2tree(y ~ A + B + Z1 + Z2, data = as.data.frame(sim),
                     unit = "i", time = "t", engine = "ranger", ntree = 400,
                     setting_between = SET_B, setting_within = SET_W, seed = 7)
  # between x within interaction stress test (additivity violation)
  inter_one <- function(gamma) {
    s2 <- sim
    s2$y <- 50 + 2.5*s2$A + 3*s2$Bw*(1 + gamma*scale(s2$A)[,1]) + rnorm(nrow(s2), 0, 1.5)
    mi <- suppressWarnings(panel_e2tree(y ~ A + B + Z1 + Z2, data = as.data.frame(s2),
                     unit = "i", time = "t", engine = "ranger", ntree = 400,
                     setting_between = SET_B, setting_within = SET_W, seed = 7))
    data.frame(gamma = gamma, r2_panel = mi$r2_panel)
  }
  inter <- do.call(rbind, lapply(c(0, 0.5, 1), inter_one))
  list(icc_sim = icc_sim,
       pooled = list(vars = ps$vars, fid = ps$fid_vs_ens),
       between = list(vars = ms$between$variables, fid = ms$between$fidelity),
       within = list(vars = ms$within$variables, fid = ms$within$fidelity),
       r2_panel_sim = ms$r2_panel, interaction = inter)
})
write.csv(s8$interaction, file.path(OUT, "sim_interaction.csv"), row.names = FALSE)
sim_tab <- data.frame(
  explanation = c("pooled (single tree)", "between", "within"),
  variables = c(paste(sort(s8$pooled$vars), collapse = ", "),
                paste(sort(s8$between$vars), collapse = ", "),
                paste(sort(s8$within$vars), collapse = ", ")),
  fidelity = round(c(s8$pooled$fid, s8$between$fid, s8$within$fid), 3))
write.csv(sim_tab, file.path(OUT, "sim_known_roles.csv"), row.names = FALSE)
put("icc_sim", s8$icc_sim, 2)
put("sim_pooled_vars",  paste(sort(s8$pooled$vars), collapse = ", "))
put("sim_between_vars", paste(sort(s8$between$vars), collapse = ", "))
put("sim_within_vars",  paste(sort(s8$within$vars), collapse = ", "))
put("sim_pooled_fid",   s8$pooled$fid, 3)
put("sim_between_fid",  s8$between$fid, 3)
put("sim_within_fid",   s8$within$fid, 3)
put("sim_inter_r2_g0",  s8$interaction$r2_panel[1], 3)
put("sim_inter_r2_g05", s8$interaction$r2_panel[2], 3)
put("sim_inter_r2_g1",  s8$interaction$r2_panel[3], 3)

# ---------------------------------------------------------------------------
# S9 — within importance across seeds (stability of the ranking)
# ---------------------------------------------------------------------------
s9 <- step("S9_vimp_seeds", {
  one <- function(seed) {
    m <- suppressWarnings(
      panel_e2tree(form_solv, data = md, unit = UNIT, time = TIME,
                   engine = "ranger", ntree = 500,
                   setting_between = SET_B, setting_within = SET_W, seed = seed))
    vi <- m$within$varimp$vimp
    if (is.null(vi)) return(NULL)
    vi$rel <- 100 * vi$MeanImpurityDecrease / sum(vi$MeanImpurityDecrease)
    data.frame(seed = seed, Variable = vi$Variable, rel = vi$rel)
  }
  do.call(rbind, lapply(c(11, 22, 33, 44, 55), one))
})
vimp_seeds <- s9 %>% group_by(Variable) %>%
  summarise(mean_rel = mean(rel), sd_rel = sd(rel), n_seeds = n(), .groups = "drop") %>%
  arrange(desc(mean_rel))
write.csv(vimp_seeds, file.path(OUT, "within_vimp_seeds.csv"), row.names = FALSE)
for (i in 1:min(3, nrow(vimp_seeds)))
  put(paste0("vimp_top", i), sprintf("%s (%.1f+-%.1f%%)",
      vimp_seeds$Variable[i], vimp_seeds$mean_rel[i], vimp_seeds$sd_rel[i]))

# ---------------------------------------------------------------------------
# S10 — robustness: reduced eight-feature set (within-relevant features only)
# ---------------------------------------------------------------------------
s10 <- step("S10_reduced", {
  feats8 <- c("roa_pbt","ebitda_margin","cf_margin","icover",
              "current_ratio","liquidity_ratio","wc_gap","stock_turnover")
  m8 <- suppressWarnings(
    panel_e2tree(reformulate(feats8, OUTCOME), data = md, unit = UNIT, time = TIME,
                 engine = "ranger", ntree = 500,
                 setting_between = SET_B, setting_within = SET_W, seed = 123))
  list(within_fid = m8$within$fidelity, within_signal = m8$within_signal)
})
put("reduced_within_fid", s10$within_fid, 3)

# proximity negative-control figure (Fig. shortcuts panel b)
th0 <- theme_minimal(base_size = 12) + theme(panel.grid.minor = element_blank())
p_prox <- ggplot(s4$sims, aes(icc_emp, prox_between)) +
  geom_line(colour = "black", linewidth = .8) + geom_point(size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  ylim(0, 1) + xlim(0, 1) + th0 +
  labs(x = "true ICC of the outcome", y = "between share of pooled proximity")
ggsave(file.path(FEN, "proximity_decomp.png"), p_prox, width = 6.4, height = 4.2, dpi = 140)

# ---------------------------------------------------------------------------
# Tables and figures
# ---------------------------------------------------------------------------
nec <- data.frame(
  approach = c("(a) pooled e2tree, raw features",
               "(a') pooled, capacity-matched",
               "(c) single e2tree, Mundlak-augmented",
               "(d) two-representation (between+within)"),
  R2_vs_y = c(numbers[["pooled_r2_vs_y"]], numbers[["matched_r2_vs_y"]],
              numbers[["aug_r2_vs_y"]], numbers[["panel_r2"]]),
  R2_within = c(numbers[["pooled_r2_within"]], numbers[["matched_r2_within"]],
                numbers[["aug_r2_within"]], numbers[["decomp_r2_within"]]),
  R2_between = c(numbers[["pooled_r2_between"]], numbers[["matched_r2_between"]],
                 numbers[["aug_r2_between"]], numbers[["decomp_r2_between"]]),
  leaves = c(numbers[["pooled_leaves"]], numbers[["matched_leaves"]],
             s3$leaves, numbers[["decomp_leaves_total"]]))
write.csv(nec, file.path(OUT, "necessity.csv"), row.names = FALSE)

fid_tab <- data.frame(
  component = c("pooled (single tree)", "between", "within", "panel (between+within)"),
  metric = c("outcome variance recovered (R2)", "fidelity vs full ensemble",
             "fidelity vs full ensemble", "outcome variance recovered (R2)"),
  value = c(numbers[["pooled_r2_vs_y"]], numbers[["between_fidelity"]],
            numbers[["within_fidelity"]], numbers[["panel_r2"]]))
write.csv(fid_tab, file.path(OUT, "fidelity_table.csv"), row.names = FALSE)

th <- theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 11), panel.grid.minor = element_blank())
nec_long <- nec %>% dplyr::filter(!is.na(R2_vs_y)) %>%
  mutate(approach = factor(approach, levels = rev(approach))) %>%
  tidyr::pivot_longer(c(R2_vs_y, R2_within), names_to = "metric", values_to = "r2") %>%
  mutate(metric = factor(ifelse(metric == "R2_vs_y", "total R² (observed solvency)",
                                "within-component R²"),
                         levels = c("total R² (observed solvency)", "within-component R²")))
p_nec <- ggplot(nec_long, aes(approach, r2, fill = r2 > 0)) +
  geom_col(width = .65) + geom_hline(yintercept = 0, colour = "grey40", linewidth = .3) +
  geom_text(aes(label = sprintf("%.2f", r2), hjust = ifelse(r2 >= 0, -0.15, 1.1)), size = 3.2) +
  coord_flip() + facet_wrap(~metric) +
  scale_fill_manual(values = c(`TRUE` = "#1b9e77", `FALSE` = "#d95f02"), guide = "none") +
  ylim(-0.6, 1) + th +
  labs(x = NULL, y = "honest R² (1 − SSE/SST) vs observed solvency")
ggsave(file.path(FEN, "e2tree_necessity.png"), p_nec, width = 10, height = 3.8, dpi = 140)

# SHAP quantities quoted in Section "Between dominance degrades model-agnostic
# explanations": read from the CSVs written by code/17_shap_panel.R (if present),
# so those numbers are macro-generated like every other quantity in the text.
shap_sum <- file.path(BASE, "output", "shap_panel_summary.csv")
shap_pan <- file.path(BASE, "output", "shap_panel.csv")
shap_sim <- file.path(BASE, "output", "shap_sim.csv")
if (all(file.exists(shap_sum, shap_pan, shap_sim))) {
  ss <- read.csv(shap_sum); sp <- read.csv(shap_pan); sm <- read.csv(shap_sim)
  raw <- sp[sp$stage == "raw pooled", ]
  raw <- raw[order(-raw$shap_importance_pct), ]
  wn  <- sp[sp$stage == "within representation", ]
  wn  <- wn[order(-wn$shap_importance_pct), ]
  put("shap_icc_mean", ss$mean_shap_icc_raw, 2)
  put("shap_icc_min",  floor(min(raw$shap_icc, na.rm = TRUE) * 100) / 100, 2)  # floor keeps "at x or above" true
  put("shap_raw_top1", round(raw$shap_importance_pct[1]))
  put("shap_raw_top2", round(raw$shap_importance_pct[2]))
  put("shap_raw_top3", round(raw$shap_importance_pct[3]))
  put("shap_wn_top1",  round(wn$shap_importance_pct[1]))
  put("shap_wn_top2",  round(wn$shap_importance_pct[2]))
  put("shap_wn_top3",  round(wn$shap_importance_pct[3]))
  put("shap_sim_raw",  sm$within_share_pct[1], 1)
  put("shap_sim_aug",  sm$within_share_pct[2], 1)
  put("shap_sim_sep",  sm$within_share_pct[3], 1)
}

# Proposition 1 on the simulation: regenerate the identical synthetic panel (seed 123,
# matching S8) and read off its leading-contrast c_B and threshold ICC*.
sim_cb <- local({
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    old <- get(".Random.seed", .GlobalEnv); on.exit(assign(".Random.seed", old, .GlobalEnv))
  }
  set.seed(123); nc <- 60; Tt <- 12
  s  <- expand.grid(i = seq_len(nc), t = seq_len(Tt))
  Ai <- runif(nc, 0, 10); Bbar <- runif(nc, 0, 5)
  s$A <- Ai[s$i]; s$B <- Bbar[s$i] + rnorm(nrow(s), 0, 2)
  s$Z1 <- rnorm(nrow(s)); s$Z2 <- runif(nrow(s), 0, 10)
  s$Bw <- ave(s$B, s$i, FUN = function(z) z - mean(z))
  s$y  <- 50 + 2.5 * s$A + 3 * s$Bw + rnorm(nrow(s), 0, 1.5)
  um   <- aggregate(s[c("y", "A", "B", "Z1", "Z2")], by = list(.unit = s$i), FUN = mean)
  R2   <- leading_between_gain(um, "y", c("A", "B", "Z1", "Z2"))
  list(cB = 12 * R2, iccstar = 1 / (1 + R2))          # K = 12 (ICC* is K-free)
})
put("c_b_sim",      sim_cb$cB, 1)
put("icc_star_sim", sim_cb$iccstar, 3)

mast <- data.frame(key = names(numbers),
                   value = vapply(numbers, fmt_number, character(1)))
write.csv(mast, file.path(OUT, "paper_numbers.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
