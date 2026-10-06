# 02 — Comparative screening of the 4 candidate ORBIS applications
#
# For each idea we ask two questions the paper cares about:
#   (1) Feasibility / necessity: is the target's ICC > 0.8 (between-dominated)?
#       and how much complete-case balanced panel data survives?
#   (2) Interestingness: once decomposed a la Mundlak, do the between drivers
#       (what explains the firm's level) and the within drivers (what explains
#       its year-on-year moves) actually differ, and is there real within signal?
#       If between and within tell the same story, the decomposition adds little;
#       if they diverge, the panel e2tree is compelling.
#
# Each idea is evaluated on its own natural complete-case, balanced (10-year)
# universe, not one pre-filtered for another target.
#
# Inputs : data/raw/output/orbis_eu_*.tsv  (rich 43-col extract from 00_extract_orbis.sh)
# Outputs: output/idea_comparison.csv            headline table
#          output/idea_target_icc.csv            ICC of every candidate target
#          output/idea_between_within_drivers.csv top drivers per idea/dimension
#          output/plots/idea_*.png
#          reports/IDEA_COMPARISON.md is written by hand from these.

suppressPackageStartupMessages({
  library(data.table); library(randomForest); library(ggplot2)
})
set.seed(1)

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
ROOT <- if (length(script_path)) dirname(dirname(normalizePath(script_path))) else normalizePath("..")
SRC  <- file.path(ROOT, "data", "raw", "output")
OUT  <- file.path(ROOT, "output")
PLT  <- file.path(OUT, "plots")
dir.create(PLT, recursive = TRUE, showWarnings = FALSE)

YEAR_MIN <- 2013L; YEAR_MAX <- 2022L; MIN_YEARS <- 10L
N_SAMPLE <- 1000L            # firms sampled per idea for the RF interestingness step
WINS_P   <- 0.01

winsor <- function(x, p = WINS_P) {
  q <- quantile(x, c(p, 1 - p), na.rm = TRUE); pmin(pmax(x, q[1]), q[2])
}
# ICC = between share of the observed variance, sigma_B^2 / (sigma_B^2 + sigma_W^2)
# with sigma_B^2 + sigma_W^2 = total. Single convention across the pipeline, matching
# manuscript icc_of() and 01_build_orbis.R.
fast_icc <- function(x, id) {
  dt <- data.table(x = as.numeric(x), id = id)[is.finite(x)]
  g  <- dt[, .(ni = .N, mi = mean(x), ssw = sum((x - mean(x))^2)), by = id]
  N <- sum(g$ni); k <- nrow(g); if (k < 2 || N <= k) return(NA_real_)
  gm  <- sum(g$ni * g$mi) / N
  vb <- sum(g$ni * (g$mi - gm)^2) / N     # between component (share of total variance)
  vw <- sum(g$ssw) / N                    # within component
  vb / (vb + vw)
}
rf_rsq <- function(form, dat) {                 # OOB R^2
  rf <- randomForest(form, data = dat, ntree = 300)
  list(rsq = tail(rf$rsq, 1), imp = rf$importance[, 1])
}

# ---- 1. load + derive a rich firm-year base table -------------------------
fin_cols <- c("bvdid","consol","closdate","tang_fixed","intang_fixed","total_assets",
              "employees","turnover","gross_profit","pbt","cost_employees","added_value",
              "cashflow","roa_pbt","profit_margin","gross_margin","ebitda_margin",
              "interest_cover","stock_turnover","collection_days","credit_days",
              "current_ratio","liquidity_ratio","solvency_asset","gearing","oprev_per_emp")
fin <- fread(file.path(SRC, "orbis_eu_financials.tsv"), sep = "\t", quote = "",
             select = fin_cols, na.strings = c("", "n.a.", "n.s."), showProgress = FALSE)
fin[, year := as.integer(substr(closdate, 1, 4))]
fin[, country := substr(bvdid, 1, 2)]
fin <- fin[year >= YEAR_MIN & year <= YEAR_MAX &
             is.finite(total_assets) & total_assets > 0 &
             is.finite(turnover) & turnover > 0]
setorder(fin, bvdid, year, -closdate)
fin <- unique(fin, by = c("bvdid", "year"))

legal <- fread(file.path(SRC, "orbis_eu_legal.tsv"), sep = "\t", quote = "",
               na.strings = c("", "n.a."), showProgress = FALSE)
legal[, incorp_year := as.integer(substr(incorp_date, 1, 4))]
legal <- unique(legal, by = "bvdid")[, .(bvdid, incorp_year)]
fin <- merge(fin, legal, by = "bvdid", all.x = TRUE)

fin[, `:=`(
  ln_assets    = log(total_assets),
  tangibility  = 100 * tang_fixed   / total_assets,
  intang_int   = 100 * intang_fixed / total_assets,
  labor_share  = 100 * cost_employees / turnover,
  cf_margin    = 100 * cashflow / turnover,
  va_per_emp   = ifelse(is.finite(employees) & employees > 0, added_value / employees, NA_real_),
  rev_per_emp  = oprev_per_emp,
  wc_gap       = collection_days - credit_days,
  icover       = interest_cover,
  firm_age     = pmax(year - incorp_year, 0L)
)]

# ---- 2. define the four ideas ---------------------------------------------
ideas <- list(
  `Idea1_ROA` = list(
    target = "roa_pbt",
    feats  = c("ln_assets","tangibility","intang_int","gross_margin","labor_share",
               "va_per_emp","stock_turnover","collection_days","gearing",
               "current_ratio","firm_age")),
  `Idea2_Leverage` = list(
    target = "gearing",
    feats  = c("ln_assets","tangibility","intang_int","roa_pbt","ebitda_margin",
               "current_ratio","liquidity_ratio","icover","stock_turnover",
               "firm_age","va_per_emp")),
  `Idea3_Productivity` = list(
    target = "va_per_emp",
    feats  = c("ln_assets","tangibility","intang_int","labor_share","ebitda_margin",
               "gearing","stock_turnover","collection_days","current_ratio",
               "firm_age","cf_margin")),
  `Idea4_Solvency` = list(
    target = "solvency_asset",
    feats  = c("roa_pbt","ebitda_margin","cf_margin","icover","current_ratio",
               "liquidity_ratio","wc_gap","stock_turnover","ln_assets",
               "firm_age","va_per_emp"))
)

# ---- 3. evaluate each idea on its own balanced complete-case universe ------
res <- list(); drv <- list()
for (nm in names(ideas)) {
  tg <- ideas[[nm]]$target; ft <- ideas[[nm]]$feats
  vars <- c(tg, ft)
  d <- fin[, c("bvdid","country","year", vars), with = FALSE]
  d <- d[complete.cases(d)]
  d[, nyr := .N, by = bvdid]
  d <- d[nyr >= MIN_YEARS]                          # balanced decade
  d <- d[, if (sd(get(tg)) > 1e-8) .SD, by = bvdid] # need within variation in target
  if (length(unique(d$bvdid)) < 30) next   # too few firms to be informative

  # feasibility on the FULL idea-specific pool
  icc_t <- fast_icc(d[[tg]], d$bvdid)
  icc_f <- sapply(ft, function(v) fast_icc(d[[v]], d$bvdid))

  # winsorise, then sample firms for the RF interestingness step
  for (v in vars) d[[v]] <- winsor(d[[v]])
  fs <- unique(d$bvdid); if (length(fs) > N_SAMPLE) fs <- sample(fs, N_SAMPLE)
  s  <- d[bvdid %in% fs]

  # Mundlak transforms
  s[, paste0(vars, "_bw") := lapply(.SD, mean), .SDcols = vars, by = bvdid]
  s[, paste0(vars, "_wn") := lapply(.SD, function(z) z - mean(z)), .SDcols = vars, by = bvdid]

  bw <- unique(s[, c(paste0(vars, "_bw")), with = FALSE])     # one row per firm
  wn <- s[, c(paste0(vars, "_wn")), with = FALSE]
  names(bw) <- vars; names(wn) <- vars
  # drop within-degenerate predictors
  keepw <- c(tg, ft[sapply(ft, function(v) sd(wn[[v]]) > 1e-8)])
  wn <- wn[, ..keepw]

  fb <- as.formula(paste(tg, "~ ."))
  pooled <- rf_rsq(fb, as.data.frame(s[, c(tg, ft), with = FALSE]))
  betw   <- rf_rsq(fb, as.data.frame(bw))
  with_  <- rf_rsq(fb, as.data.frame(wn))

  # driver divergence between the between and within importances (shared feats)
  shared <- intersect(names(betw$imp), names(with_$imp)); shared <- setdiff(shared, tg)
  rcorr <- suppressWarnings(cor(rank(betw$imp[shared]), rank(with_$imp[shared]),
                                method = "spearman"))
  top3b <- names(sort(betw$imp[shared], decreasing = TRUE))[1:min(3, length(shared))]
  top3w <- names(sort(with_$imp[shared], decreasing = TRUE))[1:min(3, length(shared))]
  jac   <- length(intersect(top3b, top3w)) / length(union(top3b, top3w))

  res[[nm]] <- data.table(
    idea = nm, target = tg,
    n_firms = length(unique(d$bvdid)), n_obs = nrow(d),
    ICC_target = round(icc_t, 3),
    pred_ICC_med = round(median(icc_f, na.rm = TRUE), 3),
    R2_pooled = round(pooled$rsq, 3),
    R2_between = round(betw$rsq, 3),
    R2_within  = round(with_$rsq, 3),
    drivers_rankcorr = round(rcorr, 2),    # low = between & within driven differently
    top3_overlap = round(jac, 2)           # low = different top drivers
  )
  drv[[nm]] <- data.table(
    idea = nm,
    between_top3 = paste(top3b, collapse = ", "),
    within_top3  = paste(top3w, collapse = ", ")
  )
}

comp <- rbindlist(res)
fwrite(comp, file.path(OUT, "idea_comparison.csv"))
fwrite(rbindlist(drv), file.path(OUT, "idea_between_within_drivers.csv"))

# ---- 4. plots --------------------------------------------------------------
th <- theme_minimal(base_size = 12)
p1 <- ggplot(comp, aes(reorder(idea, ICC_target), ICC_target)) +
  geom_col(aes(fill = ICC_target >= 0.8), width = .7) +
  geom_hline(yintercept = 0.8, linetype = 2) +
  scale_fill_manual(values = c(`TRUE` = "#1b9e77", `FALSE` = "#d95f02"), guide = "none") +
  coord_flip() + ylim(0, 1) + th +
  labs(x = NULL, y = "ICC (between share of variance)")
ggsave(file.path(PLT, "idea_icc.png"), p1, width = 7.5, height = 4.5, dpi = 130)

cm <- melt(comp[, .(idea, R2_between, R2_within)], id.vars = "idea")
p2 <- ggplot(cm, aes(idea, value, fill = variable)) +
  geom_col(position = "dodge", width = .7) + th +
  scale_fill_manual(values = c("#7570b3", "#e7298a"),
                    labels = c("between (levels)", "within (changes)"), name = NULL) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
  labs(x = NULL, y = expression(R^2))
ggsave(file.path(PLT, "idea_between_within_R2.png"), p2, width = 7.5, height = 4.5, dpi = 130)
