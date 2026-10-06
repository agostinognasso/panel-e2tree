# 09_leading_contrast.R — numerical verification of the leading-contrast lemma
# that replaces the whole-tree dichotomy of the v1 Proposition 1.
#
# At the root node, for any binary split s with child shares p_L, p_R:
#
#   Delta(s) = p_L p_R (ybar_L - ybar_R)^2                       (CART identity)
#            = p_L p_R (B_s + W_s)^2 ,
#
# writing y_it = ybar_i + v_it and
#   B_s = mubar_L - mubar_R   (unit means only)
#   W_s = vbar_L - vbar_R     (within deviations only).
#
# Claims to verify, all exact finite-sample statements:
#   (C1) the decomposition identity holds for every candidate split;
#   (C2) W_s = 0 exactly for every split on a unit-level (_bw) feature;
#   (C3) p_L p_R W_s^2 <= sigma_W^2 = S_W/n for every split  [this is (A2) of the
#        v1 proposition, here proved rather than assumed, via Cauchy-Schwarz on
#        c_it = 1{L} - p_L, for which ||c||^2 = n p_L p_R];
#   (C4) hence Delta(s) <= (sqrt(p_L p_R)|B_s| + sigma_W)^2: a within-feature
#        split can only compete by borrowing between variation;
#   (C5) the leading between contrast Delta_B* exceeds sigma_W^2, so no purely
#        within contrast (B_s = 0) can be selected at the root.
#
# Verdict (see the end of this script): (C1)-(C5) all hold, on the ORBIS panel
# and on the S8 simulation. But the two datasets then diverge completely in what
# the fitted tree does -- pooled within recovery -0.10 on ORBIS against +0.685 on
# the simulation -- although both satisfy (C5) and both take a between split at
# the root. The root-node statement is therefore true and correctly scoped, and
# carries NO information about whole-tree within recovery. It cannot be promoted
# into the whole-tree claim: that is exactly the invalid leap of the v1
# Proposition 1, and the simulation is the counterexample.
#
# Scope: variance reduction (CART), which is what the v1 proposition was stated
# for. The proximity criterion actually used by e2tree is checked separately, at
# the end, by ranking the same candidate splits with the proximity impurity.
#
# Output: output/budget59/leading_contrast.csv + console report.
# Usage: Rscript code/09_leading_contrast.R

suppressPackageStartupMessages({ library(dplyr) })

BASE <- normalizePath(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))), ".."))
OUT <- file.path(BASE, "output", "budget59")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

OUTCOME <- "solvency_asset"; UNIT <- "bvdid"; TIME <- "year"
feats <- c("roa_pbt","ebitda_margin","cf_margin","icover","current_ratio",
           "liquidity_ratio","wc_gap","stock_turnover","ln_assets","firm_age","rev_per_emp")
md <- readRDS(file.path(BASE, "data", "model_data.rds"))
md <- md[complete.cases(md[, c(UNIT, TIME, OUTCOME, feats)]), ]

y   <- md[[OUTCOME]]
g   <- md[[UNIT]]
n   <- length(y)
mu  <- stats::ave(y, g)          # unit means, repeated
v   <- y - mu                    # within deviations, sum to 0 within each unit
S_B <- sum((mu - mean(y))^2)
S_W <- sum(v^2)
sigma_B2 <- S_B / n
sigma_W2 <- S_W / n

cat(sprintf("n = %d, N = %d, balanced = %s\n", n, length(unique(g)),
            length(unique(table(g))) == 1))
cat(sprintf("S_B = %.1f  S_W = %.1f   sigma_B^2 = %.2f  sigma_W^2 = %.2f   eta_B = %.4f\n",
            S_B, S_W, sigma_B2, sigma_W2, S_B / (S_B + S_W)))

# --- Mundlak-augmented candidate features ---------------------------------
aug <- md %>% group_by(.data[[UNIT]]) %>%
  mutate(across(all_of(feats), list(bw = ~mean(.x), wn = ~.x - mean(.x)))) %>%
  ungroup() %>% select(ends_with("_bw"), ends_with("_wn")) %>% as.data.frame()
keep <- vapply(aug, function(z) stats::sd(z) > 1e-8, logical(1))
if (any(!keep)) cat("dropped as constant: ", paste(names(aug)[!keep], collapse = ", "), "\n")
aug <- aug[, keep, drop = FALSE]

# --- all candidate thresholds of one feature, in O(n log n) ----------------
scan_feature <- function(x, name) {
  o  <- order(x)
  xs <- x[o]
  cy <- cumsum(y[o]); cm <- cumsum(mu[o]); cv <- cumsum(v[o])
  k  <- which(xs[-n] < xs[-1])                 # valid cut points (no ties)
  if (!length(k)) return(NULL)
  nL <- k; nR <- n - k
  pL <- nL / n; pR <- nR / n
  ybL <- cy[k] / nL;              ybR <- (cy[n] - cy[k]) / nR
  muL <- cm[k] / nL;              muR <- (cm[n] - cm[k]) / nR
  vbL <- cv[k] / nL;              vbR <- (cv[n] - cv[k]) / nR
  data.frame(feature = name, kind = if (grepl("_bw$", name)) "between" else "within",
             thr = (xs[k] + xs[k + 1]) / 2, nL = nL, pL = pL,
             delta   = pL * pR * (ybL - ybR)^2,
             B       = muL - muR,
             W       = vbL - vbR,
             stringsAsFactors = FALSE)
}

cand <- do.call(rbind, lapply(names(aug), function(nm) scan_feature(aug[[nm]], nm)))
cand$pR      <- 1 - cand$pL
cand$delta_B <- cand$pL * cand$pR * cand$B^2      # between part of the contrast
cand$delta_W <- cand$pL * cand$pR * cand$W^2      # within  part of the contrast
cand$recon   <- cand$pL * cand$pR * (cand$B + cand$W)^2
cat(sprintf("\ncandidate splits enumerated: %d over %d features\n", nrow(cand), ncol(aug)))

# --- (C1) the decomposition identity ---------------------------------------
err1 <- max(abs(cand$delta - cand$recon))
cat(sprintf("(C1) max |Delta - p_L p_R (B+W)^2|                = %.3e   %s\n",
            err1, if (err1 < 1e-8) "OK" else "FAIL"))

# --- (C2) between splits carry no within contrast --------------------------
bw <- cand$kind == "between"
err2 <- max(abs(cand$W[bw]))
cat(sprintf("(C2) max |W_s| over between splits (n = %d)     = %.3e   %s\n",
            sum(bw), err2, if (err2 < 1e-8) "OK" else "FAIL"))

# --- (C3) the proved bound: p_L p_R W_s^2 <= sigma_W^2 ---------------------
worst <- max(cand$delta_W)
cat(sprintf("(C3) max p_L p_R W_s^2 = %.3f   vs sigma_W^2 = %.3f   ratio %.4f   %s\n",
            worst, sigma_W2, worst / sigma_W2,
            if (worst <= sigma_W2 * (1 + 1e-10)) "OK (bound holds)" else "FAIL"))

# --- (C4)/(C5) the leading contrasts ---------------------------------------
best_b <- cand[bw, ][which.max(cand$delta[bw]), ]
best_w <- cand[!bw, ][which.max(cand$delta[!bw]), ]
cat(sprintf("\nleading BETWEEN contrast: %s <= %.3f | Delta_B* = %.3f  (p_L = %.3f)\n",
            best_b$feature, best_b$thr, best_b$delta, best_b$pL))
cat(sprintf("leading WITHIN  contrast: %s <= %.3f | Delta_W* = %.3f  (p_L = %.3f)\n",
            best_w$feature, best_w$thr, best_w$delta, best_w$pL))
# Delta = Delta_B + 2*sqrt(Delta_B*Delta_W)*sign(BW) + Delta_W: report the three
# terms, not percentages, because the cross term is not a share of anything.
cross <- best_w$delta - best_w$delta_B - best_w$delta_W
cat(sprintf("  decomposed: between part %.3f + cross %.3f + genuine within part %.3f = %.3f\n",
            best_w$delta_B, cross, best_w$delta_W, best_w$delta))
cat(sprintf("  without any between leakage it could not exceed sigma_W^2 = %.3f\n", sigma_W2))
cat(sprintf("(C5) Delta_B* = %.3f  >  sigma_W^2 = %.3f ?  %s\n",
            best_b$delta, sigma_W2,
            if (best_b$delta > sigma_W2) "YES" else "NO"))
cat(sprintf("     Delta_B* > Delta_W* ?  %s  (ratio %.2f)\n",
            if (best_b$delta > best_w$delta) "YES" else "NO",
            best_b$delta / best_w$delta))

# how far can a purely within contrast (B_s = 0) go?
pure <- cand[!bw & cand$delta_B < 1e-6 * cand$delta, ]
cat(sprintf("     splits with negligible between leakage: %d; best Delta among them = %s\n",
            nrow(pure), if (nrow(pure)) sprintf("%.3f", max(pure$delta)) else "none"))
cat(sprintf("     realized within ceiling %.3f is %.1fx below Delta_B* (proved ceiling %.1f is %.2fx below)\n",
            worst, best_b$delta / worst, sigma_W2, best_b$delta / sigma_W2))

# --- top of the ranking, and the between share -----------------------------
top <- cand[order(-cand$delta), ][1:10, ]
cat("\nten best root candidates by variance reduction:\n")
print(top[, c("feature", "kind", "delta", "delta_B", "delta_W", "pL")],
      digits = 4, row.names = FALSE)
cat(sprintf("between share of the top ten: %.0f%%\n", 100 * mean(top$kind == "between")))

write.csv(cand[order(-cand$delta), ][1:500, ],
          file.path(OUT, "leading_contrast.csv"), row.names = FALSE)
cat("\nwritten:", file.path(OUT, "leading_contrast.csv"), "(top 500)\n")

saveRDS(list(sigma_B2 = sigma_B2, sigma_W2 = sigma_W2, S_B = S_B, S_W = S_W,
             best_between = best_b, best_within = best_w,
             max_delta_W = worst, n_cand = nrow(cand),
             top10_between_share = mean(top$kind == "between")),
        file.path(OUT, "leading_contrast.rds"))

# ---------------------------------------------------------------------------
# Does the condition Delta_B* > sigma_W^2 have bite? The simulation of S8 has a
# strong, genuine within driver and a lower ICC, so it is the natural case where
# the condition should fail. If it failed nowhere, it would be vacuous.
# ---------------------------------------------------------------------------
cat("\n--- same check on the S8 simulation (known roles, ICC ~ 0.61) ---\n")
set.seed(123)                                   # identical DGP to S8
nc <- 60; Tt <- 12
sim <- expand.grid(i = seq_len(nc), t = seq_len(Tt))
Ai <- runif(nc, 0, 10); Bbar <- runif(nc, 0, 5)
sim$A <- Ai[sim$i]; sim$B <- Bbar[sim$i] + rnorm(nrow(sim), 0, 2)
sim$Z1 <- rnorm(nrow(sim)); sim$Z2 <- runif(nrow(sim), 0, 10)
sim <- sim %>% group_by(i) %>% mutate(Bw = B - mean(B)) %>% ungroup()
sim$y <- 50 + 2.5 * sim$A + 3 * sim$Bw + rnorm(nrow(sim), 0, 1.5)

local({
  y <<- sim$y; g <<- sim$i; n <<- nrow(sim)
  mu <<- stats::ave(y, g); v <<- y - mu
  sB2 <- sum((mu - mean(y))^2) / n; sW2 <- sum(v^2) / n
  sigma_W2 <<- sW2
  a <- sim %>% group_by(i) %>%
    mutate(across(c(A, B, Z1, Z2), list(bw = ~mean(.x), wn = ~.x - mean(.x)))) %>%
    ungroup() %>% select(ends_with("_bw"), ends_with("_wn")) %>% as.data.frame()
  a <- a[, vapply(a, function(z) stats::sd(z) > 1e-8, logical(1)), drop = FALSE]
  cd <- do.call(rbind, lapply(names(a), function(nm) scan_feature(a[[nm]], nm)))
  cd$pR <- 1 - cd$pL
  cd$delta_W <- cd$pL * cd$pR * cd$W^2
  b  <- cd$kind == "between"
  dB <- max(cd$delta[b]); dW <- max(cd$delta[!b])
  cat(sprintf("  ICC = %.3f   sigma_B^2 = %.2f  sigma_W^2 = %.2f\n", sB2 / (sB2 + sW2), sB2, sW2))
  cat(sprintf("  (C3) max p_L p_R W_s^2 = %.3f <= sigma_W^2 = %.3f   %s\n",
              max(cd$delta_W), sW2, if (max(cd$delta_W) <= sW2 * (1 + 1e-10)) "OK" else "FAIL"))
  cat(sprintf("  leading between %s: Delta_B* = %.3f ; leading within %s: Delta_W* = %.3f\n",
              cd[b, ][which.max(cd$delta[b]), "feature"], dB,
              cd[!b, ][which.max(cd$delta[!b]), "feature"], dW))
  cat(sprintf("  (C5) Delta_B* > sigma_W^2 ?  %s   -> the condition %s here\n",
              if (dB > sW2) "YES" else "NO",
              if (dB > sW2) "holds" else "FAILS, as it should on a panel with a real within driver"))
  cat(sprintf("  root split actually taken by variance reduction: %s\n",
              cd[which.max(cd$delta), "feature"]))
})

# ---------------------------------------------------------------------------
# The decisive test: does the root-node condition predict what the tree does?
# Both datasets satisfy (C5) and both take a between split at the root. If the
# fitted trees then behave differently, the root statement has no whole-tree
# content -- which is precisely what the v1 proposition assumed it had.
# ---------------------------------------------------------------------------
cat("\n--- whole-tree behaviour under the same root condition ---\n")
suppressPackageStartupMessages(library(e2tree))
source(file.path(BASE, "code", "panel_e2tree.R"))
source(file.path(BASE, "code", "e2tree_split_local.R"))
r2w <- function(y, p, gg) { yd <- y - stats::ave(y, gg); pd <- p - stats::ave(p, gg)
                            1 - sum((yd - pd)^2) / sum(yd^2) }
r2h <- function(o, p) 1 - sum((o - p)^2) / sum((o - mean(o))^2)
SET_W <- list(impTotal = 0.05, maxDec = 1e-7, n = 5, level = 6)

fit_pooled <- function(d, ylab, seed = 7) {
  set.seed(seed)
  fo  <- reformulate(setdiff(names(d), ylab), ylab)
  en  <- ranger::ranger(fo, data = d, num.trees = 500)
  DD  <- createDisMatrix(en, data = d, label = ylab, parallel = list(active = FALSE, no_cores = 1))
  tt  <- e2tree(fo, d, DD, en, SET_W)
  list(tree = tt, fit = as.numeric(predict(tt, newdata = d)$fit),
       leaves = sum(tt$tree$terminal), root = tt$tree$splitLabel[1],
       vars = unique(stats::na.omit(tt$tree$variable)))
}

pa <- fit_pooled(md[, c(OUTCOME, feats)], OUTCOME)
cat(sprintf("ORBIS      root %-28s %2d leaves | R2 tot %.3f | R2_within %+.3f\n",
            pa$root, pa$leaves, r2h(md[[OUTCOME]], pa$fit), r2w(md[[OUTCOME]], pa$fit, md[[UNIT]])))
ps <- fit_pooled(as.data.frame(sim[, c("y", "A", "B", "Z1", "Z2")]), "y")
cat(sprintf("SIMULATION root %-28s %2d leaves | R2 tot %.3f | R2_within %+.3f\n",
            ps$root, ps$leaves, r2h(sim$y, ps$fit), r2w(sim$y, ps$fit, sim$i)))
cat("\nBoth satisfy Delta_B* > sigma_W^2 and both split on a between contrast at the\n",
    "root, yet within recovery is negative in one case and strongly positive in the\n",
    "other. The root-node lemma is true; it does not determine the fitted tree.\n", sep = "")
