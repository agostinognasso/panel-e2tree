# 05_budget59.R — Step 4A of the Carmela Iorio revision: pooled single e2tree
# grown under a leaf budget matched to the decomposed procedure (12 + 47 = 59).
#
# e2tree exposes no leaf-budget argument: tree size is controlled only through the
# stopping rules (impTotal, maxDec, n, level). This script therefore
#   (i)  rebuilds exactly the pooled ensemble of fit_single_e2() (seed 7, ranger,
#        500 trees, md[, c(OUTCOME, feats)]),
#   (ii) computes the dissimilarity matrix once and caches it,
#   (iii) sweeps the stopping rules over that single D, restoring the RNG state
#         to the post-createDisMatrix state before every e2tree() call, so each
#         configuration is bit-identical to what 03_panel_e2tree_orbis.R would
#         have produced by refitting from scratch,
#   (iv) reports the effective leaf count and maximum depth of every configuration,
#        outcome-recovery decomposition (r2_within / r2_between).
#
# The first two rows of the grid reproduce the published `pooled` (14 leaves) and
# `matched` (18 leaves) configurations: they are the validation that reusing one
# cached D is equivalent to refitting.
#
# Nothing under output/repro/ is touched. Results go to output/budget59/.
# Usage: Rscript code/05_budget59.R

suppressPackageStartupMessages({ library(dplyr); library(e2tree) })

BASE <- normalizePath(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))), ".."))
stopifnot(file.exists(file.path(BASE, "data", "model_data.rds")))

OUT <- file.path(BASE, "output", "budget59")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# D is ~512 MB: cache it outside the project tree.
DCACHE <- Sys.getenv("BUDGET59_DCACHE", unset = file.path(tempdir(), "D_pooled.rds"))

log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

# --- helpers copied verbatim from 03_panel_e2tree_orbis.R ------------------
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
max_depth <- function(tr) max(floor(log2(tr$tree$node)))

# --- data, identical to S0/S2 ---------------------------------------------
OUTCOME <- "solvency_asset"; UNIT <- "bvdid"; TIME <- "year"
md    <- readRDS(file.path(BASE, "data", "model_data.rds"))
feats <- c("roa_pbt","ebitda_margin","cf_margin","icover","current_ratio",
           "liquidity_ratio","wc_gap","stock_turnover","ln_assets","firm_age","rev_per_emp")
md <- md[complete.cases(md[, c(UNIT, TIME, OUTCOME, feats)]), ]
d   <- md[, c(OUTCOME, feats)]
y   <- md[[OUTCOME]]; g <- md[[UNIT]]
log_msg("data: ", nrow(d), " obs, ", length(unique(g)), " firms, ", length(feats), " features")

# --- ensemble + dissimilarity matrix, computed once ------------------------
form <- reformulate(setdiff(names(d), OUTCOME), OUTCOME)
set.seed(7)                                    # same seed as fit_single_e2()
ens <- ranger::ranger(form, data = d, num.trees = 500)
log_msg("ranger fitted (OOB R2 = ", round(ens$r.squared, 4), ")")

if (file.exists(DCACHE)) {
  log_msg("loading cached D from ", DCACHE)
  cached <- readRDS(DCACHE)
  D <- cached$D; rng_state <- cached$rng_state
  rm(cached)
} else {
  log_msg("computing dissimilarity matrix (this is the slow step)...")
  t0 <- Sys.time()
  D <- createDisMatrix(ens, data = d, label = OUTCOME,
                       parallel = list(active = FALSE, no_cores = 1))
  log_msg("D done in ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min; dim ",
          paste(dim(D), collapse = " x "))
  rng_state <- .Random.seed                    # state fit_single_e2 would hand to e2tree()
  saveRDS(list(D = D, rng_state = rng_state), DCACHE, compress = FALSE)
  log_msg("D cached to ", DCACHE)
}

ep <- get_ensemble_predictions(ens, d, type = "regression")

# --- grid of stopping rules ------------------------------------------------
# Two published configurations first: they are the validation that reusing one
# cached D is equivalent to refitting from scratch.
# Then the sweep. A first pass showed that below impTotal = 0.01 the impurity
# threshold stops binding altogether and the tree size is set by the depth cap
# (level) and the minimum node size (n). The sweep is therefore over (level, n)
# with impTotal held at a non-binding 1e-4 and maxDec at the `matched` value.
grid <- list(
  list(tag = "pooled (published)",  impTotal = 0.05, maxDec = 1e-7, n = 5, level = 6),
  list(tag = "matched (published)", impTotal = 0.02, maxDec = 1e-8, n = 2, level = 7)
)
sweep <- expand.grid(level = c(8, 10, 11, 12, 13, 14, 15, 16),
                     n = c(2, 5, 10, 15, 20, 25, 30, 40, 50),
                     KEEP.OUT.ATTRS = FALSE)
sweep <- sweep[order(sweep$level, sweep$n), ]
for (k in seq_len(nrow(sweep)))
  grid[[length(grid) + 1]] <- list(tag = "budget sweep", impTotal = 1e-4,
                                   maxDec = 1e-8, n = sweep$n[k], level = sweep$level[k])

# An impurity ladder at fixed (level, n), so the claim that the total-impurity
# threshold stops binding below a certain value is backed by this file rather
# than asserted in the text.
for (it in c(0.05, 0.02, 0.01, 5e-3, 1e-3, 1e-4, 1e-5, 1e-6))
  grid[[length(grid) + 1]] <- list(tag = "impurity ladder", impTotal = it,
                                   maxDec = 1e-8, n = 2, level = 10)

rows <- list(); trees <- list()
for (i in seq_along(grid)) {
  cfg <- grid[[i]]
  setting <- cfg[c("impTotal", "maxDec", "n", "level")]
  assign(".Random.seed", rng_state, envir = .GlobalEnv)   # exact refit equivalence
  t0 <- Sys.time()
  tr <- e2tree(form, d, D, ens, setting)
  el <- as.numeric(round(difftime(Sys.time(), t0, units = "mins"), 2))
  fit <- as.numeric(predict(tr, newdata = d)$fit)
  rows[[i]] <- data.frame(
    tag = cfg$tag, impTotal = cfg$impTotal, maxDec = cfg$maxDec,
    n = cfg$n, level = cfg$level,
    leaves = n_leaves(tr), depth = max_depth(tr),
    r2_vs_y = r2_honest(y, fit), cor2_vs_y = cor(fit, y)^2,
    fid_vs_ens = cor(fit, ep)^2,
    r2_within = r2_within_of(y, fit, g), r2_between = r2_between_of(y, fit, g),
    n_vars = length(unique(na.omit(tr$tree$variable))),
    minutes = el, stringsAsFactors = FALSE)
  trees[[i]] <- list(setting = setting, leaves = n_leaves(tr), depth = max_depth(tr),
                     vars = unique(na.omit(tr$tree$variable)), fit = fit)
  log_msg(sprintf("%-20s impTotal=%-8g -> %3d leaves, depth %2d, r2=%.3f, within=%+.3f, between=%.3f (%.1f min)",
                  cfg$tag, cfg$impTotal, n_leaves(tr), max_depth(tr),
                  rows[[i]]$r2_vs_y, rows[[i]]$r2_within, rows[[i]]$r2_between, el))
  res <- do.call(rbind, rows)
  write.csv(res, file.path(OUT, "budget59_grid.csv"), row.names = FALSE)
  saveRDS(trees, file.path(OUT, "budget59_trees.rds"))
}

log_msg("=== grid complete ===")
print(do.call(rbind, rows), digits = 4)
log_msg("written to ", file.path(OUT, "budget59_grid.csv"))
