# supplement_tables.R — generate LaTeX table fragments for the Supplementary
# Material from the pipeline outputs, so no number in the supplement is typed by
# hand (same discipline as numbers_to_tex.R for the main text). Fragments are
# written to manuscript/manuscript_DSS_v2/supp/ and \input by supplementary.tex.
# Run after code/03_panel_e2tree_orbis.R and code/17_shap_panel.R.
suppressWarnings(suppressPackageStartupMessages(library(rpart)))

args <- commandArgs(trailingOnly = FALSE)
sp <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
BASE <- if (length(sp)) dirname(dirname(normalizePath(sp))) else normalizePath("..")
OUT <- file.path(BASE, "output"); REP <- file.path(OUT, "repro")
SUPP <- file.path(BASE, "manuscript", "manuscript_DSS_v2", "supp")
dir.create(SUPP, recursive = TRUE, showWarnings = FALSE)

esc  <- function(s) gsub("_", "\\\\_", s)                 # escape underscores for LaTeX text
tt   <- function(s) sprintf("\\texttt{%s}", esc(s))       # monospace variable names
row  <- function(cells) paste0(paste(cells, collapse = " & "), " \\\\")
emit <- function(file, colspec, header, body) {
  writeLines(c(sprintf("\\begin{tabular}{%s}", colspec), "\\toprule",
               row(header), "\\midrule", body, "\\bottomrule", "\\end{tabular}"),
             file.path(SUPP, file))
}
fmtnum <- function(x, d = 2) formatC(round(as.numeric(x), d), format = "f", digits = d)

## ---- S1: descriptive statistics -----------------------------------------
d <- read.csv(file.path(OUT, "descriptives.csv"), stringsAsFactors = FALSE)
emit("tab_descriptives.tex", "lrrrrr",
     c("Variable", "Mean", "SD", "P05", "Median", "P95"),
     apply(d, 1, function(r) row(c(tt(r[["variable"]]),
       fmtnum(r[["mean"]]), fmtnum(r[["sd"]]), fmtnum(r[["p05"]]),
       fmtnum(r[["median"]]), fmtnum(r[["p95"]])))))

## ---- S1: firms by country ------------------------------------------------
fc <- read.csv(file.path(OUT, "firms_by_country.csv"), stringsAsFactors = FALSE)
fc <- fc[order(-fc$n_firms), ]
emit("tab_country.tex", "lr", c("Country (ISO-2)", "Firms"),
     apply(fc, 1, function(r) row(c(r[["country"]], r[["n_firms"]]))))

## ---- S2: per-variable ICC, recomputed consistently on both samples -------
icc_of <- function(x, g) {
  m <- tapply(x, g, mean); ng <- tapply(x, g, length)
  vb <- sum(ng * (m - mean(x))^2) / length(x)
  vw <- sum((x - m[match(g, names(m))])^2) / length(x); vb / (vb + vw)
}
vars <- c("solvency_asset", "roa_pbt", "ebitda_margin", "cf_margin", "icover",
          "current_ratio", "liquidity_ratio", "wc_gap", "stock_turnover",
          "rev_per_emp", "ln_assets", "firm_age")
pl <- readRDS(file.path(BASE, "data", "orbis_eu_pool.rds"))
md <- readRDS(file.path(BASE, "data", "model_data.rds"))
feats11 <- setdiff(vars, "solvency_asset")
md <- md[complete.cases(md[, c("bvdid", "year", vars)]), ]      # the modelling sample
icc_row <- function(v) {
  ip <- icc_of(pl[[v]], pl$bvdid); is <- icc_of(md[[v]], md$bvdid)
  role <- if (v == "solvency_asset") "target"
          else if (v %in% c("ln_assets", "firm_age")) "structural control" else "predictor"
  row(c(tt(v), role, fmtnum(ip, 3), fmtnum(is, 3)))
}
emit("tab_icc.tex", "llrr",
     c("Variable", "Role", "ICC (pool)", "ICC (sample)"), vapply(vars, icc_row, ""))

## ---- S3: between-split gain profile behind c_B / ICC* --------------------
# variance-reduction CART on the unit means {(xbar_i, ybar_i)}: the object of
# Proposition 1. Each split's variance reduction as a share of total between SS.
um <- aggregate(md[c("solvency_asset", feats11)], by = list(.unit = md$bvdid), FUN = mean)
ct <- rpart.control(minsplit = 2, minbucket = 1, cp = 0, maxcompete = 0,
                    maxsurrogate = 0, xval = 0)
fit <- rpart(reformulate(feats11, "solvency_asset"), data = um, method = "anova", control = ct)
cpt <- fit$cptable; nsp <- min(11, max(cpt[, "nsplit"]))
pr  <- prune(fit, cp = cpt[max(which(cpt[, "nsplit"] <= nsp)), "CP"])
fr  <- pr$frame; dev <- setNames(fr$dev, rownames(fr)); ssb <- dev[["1"]]
internal <- as.integer(rownames(fr))[fr$var != "<leaf>"]
gains <- vapply(internal, function(m) (dev[[as.character(m)]] -
           dev[[as.character(2*m)]] - dev[[as.character(2*m+1)]]) / ssb, numeric(1))
R2 <- sort(gains, decreasing = TRUE); cum <- cumsum(R2)
emit("tab_cb_profile.tex", "rrr",
     c("Split (by gain)", "Share of between SS", "Cumulative"),
     mapply(function(i, r, c) row(c(i, fmtnum(r, 3), fmtnum(c, 3))),
            seq_along(R2), R2, cum))

## ---- S5: SHAP attribution table (raw pooled vs within) -------------------
sh <- read.csv(file.path(OUT, "shap_panel.csv"), stringsAsFactors = FALSE)
raw <- sh[sh$stage == "raw pooled", ]
wn  <- sh[sh$stage == "within representation", ]
wn$base <- sub("_wn$", "", wn$feature)
raw <- raw[order(-raw$shap_importance_pct), ]
shap_row <- function(f) {
  rr <- raw[raw$feature == f, ]; ww <- wn[wn$base == f, ]
  row(c(tt(f), fmtnum(rr$shap_importance_pct, 1), fmtnum(rr$shap_icc, 3),
        if (nrow(ww)) fmtnum(ww$shap_importance_pct, 1) else "--"))
}
emit("tab_shap.tex", "lrrr",
     c("Feature", "Raw pooled (\\%)", "Attribution ICC", "Within repr.\\ (\\%)"),
     vapply(raw$feature, shap_row, ""))

## ---- S5: within variable importance across seeds -------------------------
vi <- read.csv(file.path(REP, "within_vimp_seeds.csv"), stringsAsFactors = FALSE)
vi <- vi[order(-vi$mean_rel), ]
emit("tab_vimp.tex", "lrr",
     c("Variable", "Mean importance (\\%)", "SD (\\%)"),
     apply(vi, 1, function(r) row(c(tt(r[["Variable"]]),
       fmtnum(r[["mean_rel"]], 1), fmtnum(r[["sd_rel"]], 1)))))

## ---- S5: two-way period effects ------------------------------------------
te <- read.csv(file.path(REP, "twoway_time_effects.csv"), stringsAsFactors = FALSE)
emit("tab_timeeffects.tex", "rr", c("Year", "Period effect $\\hat\\tau_t$"),
     apply(te, 1, function(r) row(c(r[[".time"]], fmtnum(r[["solvency_asset"]], 2)))))

## ---- S5: interaction stress test -----------------------------------------
it <- read.csv(file.path(REP, "sim_interaction.csv"), stringsAsFactors = FALSE)
emit("tab_interaction.tex", "rr",
     c("Interaction strength $\\gamma$", "Panel $R^2$"),
     apply(it, 1, function(r) row(c(fmtnum(r[["gamma"]], 1), fmtnum(r[["r2_panel"]], 3)))))

## ---- S5: fidelity summary ------------------------------------------------
fd <- read.csv(file.path(REP, "fidelity_table.csv"), stringsAsFactors = FALSE)
emit("tab_fidelity.tex", "llr", c("Component", "Metric", "Value"),
     apply(fd, 1, function(r) row(c(esc(r[["component"]]), r[["metric"]],
       fmtnum(r[["value"]], 3)))))

## ---- S7: session info (verbatim copy) ------------------------------------
file.copy(file.path(REP, "sessionInfo.txt"), file.path(SUPP, "sessioninfo.txt"),
          overwrite = TRUE)
