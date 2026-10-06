# 17 — Between dominance also degrades model-agnostic explanations (SHAP)
#  The between/within problem is not specific to e2tree. On high-ICC panel data,
#  SHAP attributions are themselves between-dominated: SHAP explains "what
#  distinguishes the firms", nearly mute on "what moves them" over time.
#  Decomposing the representation a la Mundlak rescues SHAP too.
#
#  (1) Simulation with known roles (A=between, B=within, Z=noise).
#  (2) Real ORBIS solvency panel: ICC of the SHAP attributions is high (they vary
#      across firms, not over time) -> SHAP is mute on the within. The separated
#      within representation surfaces the within drivers (profitability/working-
#      capital deviations).
#
# Output: output/shap_sim.csv, output/shap_panel.csv, output/plots/shap_panel.png
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(kernelshap); library(randomForest)
})
set.seed(123)
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
ROOT <- if (length(script_path)) dirname(dirname(normalizePath(script_path))) else normalizePath("..")
DATA <- file.path(ROOT,"data"); OUT <- file.path(ROOT,"output")
PLT  <- file.path(ROOT,"output","plots"); FEN <- PLT
dir.create(PLT, recursive = TRUE, showWarnings = FALSE)

UNIT <- "bvdid"; OUTCOME <- "solvency_asset"
between_share <- function(x, grp){
  m <- tapply(x, grp, mean); ng <- tapply(x, grp, length)
  vb <- sum(ng*(m - mean(x))^2) / length(x)
  vw <- sum((x - m[match(grp, names(m))])^2) / length(x)
  if (vb + vw < 1e-12) return(NA_real_); vb/(vb+vw)
}
shap_mat <- function(rf, Xexpl, Xbg){
  ks <- kernelshap(rf, X = Xexpl, bg_X = Xbg,
                   pred_fun = function(object, X) predict(object, X), verbose = FALSE)
  as.data.frame(ks$S)
}
mean_abs <- function(S) sort(colMeans(abs(S)), decreasing = TRUE)
pct <- function(v) round(100*v/sum(v), 1)

# (1) Simulation with known roles (generic method validation) -----------------
sim_shap <- function(icc=0.9, nc=60, Tt=12){
  A <- runif(nc,0,10); Bbar <- runif(nc,0,5); Z1b <- runif(nc,0,5); Z2b <- runif(nc,0,5)
  d <- expand.grid(i=1:nc, t=1:Tt)
  d$A  <- A[d$i]; d$B <- Bbar[d$i] + rnorm(nrow(d),0,2)
  d$Z1 <- Z1b[d$i] + rnorm(nrow(d),0,2); d$Z2 <- Z2b[d$i] + rnorm(nrow(d),0,2)
  d <- d %>% group_by(i) %>% mutate(Bw = B - mean(B)) %>% ungroup()
  d$y <- 50 + 6*sqrt(icc)*scale(d$A)[,1] + 6*sqrt(1-icc)*scale(d$Bw)[,1] + rnorm(nrow(d),0,0.6)
  d <- d %>% group_by(i) %>%
    mutate(across(c(A,B,Z1,Z2), list(bw=~mean(.x), wn=~.x-mean(.x)))) %>% ungroup() %>% as.data.frame()
  dropc <- function(x) x[, sapply(x, function(z) sd(z) > 1e-8), drop=FALSE]
  bgN <- 40; exN <- min(400, nrow(d))
  bg <- d[sample(nrow(d), bgN), ]; ex <- d[sample(nrow(d), exN), ]
  Xr <- c("A","B","Z1","Z2")
  rf_r <- randomForest(reformulate(Xr,"y"), d, ntree=400); mr <- mean_abs(shap_mat(rf_r, ex[,Xr], bg[,Xr]))
  Xa <- dropc(d[, grep("_(bw|wn)$", names(d), value=TRUE)]); Xa_n <- names(Xa)
  rf_a <- randomForest(reformulate(Xa_n,"y"), cbind(y=d$y, Xa), ntree=400); ma <- mean_abs(shap_mat(rf_a, ex[,Xa_n], bg[,Xa_n]))
  dW <- d %>% group_by(i) %>% mutate(yw = y - mean(y)) %>% ungroup() %>% as.data.frame()
  Xw <- dropc(dW[, grep("_wn$", names(dW), value=TRUE)]); Xw_n <- names(Xw)
  rf_w <- randomForest(reformulate(Xw_n,"yw"), cbind(yw=dW$yw, Xw), ntree=400)
  exW <- dW[match(rownames(ex), rownames(dW)), ]; bgW <- dW[match(rownames(bg), rownames(dW)), ]
  mw <- mean_abs(shap_mat(rf_w, exW[,Xw_n], bgW[,Xw_n]))
  topname <- function(v) names(v)[1]
  share_of <- function(v, nm) if (nm %in% names(v)) unname(pct(v)[match(nm, names(v))]) else 0
  rank_of  <- function(v, nm) if (nm %in% names(v)) match(nm, names(v)) else NA_integer_
  data.frame(
    representation = c("raw features {A,B,Z1,Z2}","Mundlak-augmented, single model","within representation (separated)"),
    top_feature = c(topname(mr), topname(ma), topname(mw)),
    within_driver = c("B","B_wn","B_wn"),
    within_share_pct = c(share_of(mr,"B"), share_of(ma,"B_wn"), share_of(mw,"B_wn")),
    within_rank = c(rank_of(mr,"B"), rank_of(ma,"B_wn"), rank_of(mw,"B_wn")),
    recovered = c(ifelse(rank_of(mr,"B")==1,"yes","no"),
                  ifelse(rank_of(ma,"B_wn")==1,"yes","no"),
                  ifelse(rank_of(mw,"B_wn")==1,"yes","no")),
    stringsAsFactors = FALSE)
}
sim_tab <- sim_shap(icc=0.9)
write.csv(sim_tab, file.path(OUT,"shap_sim.csv"), row.names=FALSE)

# (2) Real ORBIS solvency panel ----------------------------------------------
panel <- readRDS(file.path(DATA,"harmonized_panel.rds"))
# within-relevant drivers: drop near-time-invariant size/age (their _wn ~ 0)
drv <- c("roa_pbt","ebitda_margin","cf_margin","icover","current_ratio",
         "liquidity_ratio","wc_gap","stock_turnover","rev_per_emp")

# (2a) Raw features: pooled RF + SHAP + ICC of the attributions
raw <- panel %>% select(all_of(UNIT), all_of(OUTCOME), all_of(drv)) %>%
  filter(if_all(everything(), ~!is.na(.x))) %>% as.data.frame()
icc_y <- between_share(raw[[OUTCOME]], raw[[UNIT]])
rf_raw <- randomForest(reformulate(drv, OUTCOME), raw[,c(OUTCOME,drv)], ntree=400)
set.seed(7)
# balanced explain subsample: up to 4 years/firm capped at ~400 rows, bg=40
ex_idx <- raw %>% mutate(.r=row_number()) %>% group_by(.data[[UNIT]]) %>%
  slice_sample(n=4) %>% ungroup() %>% pull(.r)
if (length(ex_idx) > 400) ex_idx <- sort(sample(ex_idx, 400)) else ex_idx <- sort(ex_idx)
bg_idx <- sample(nrow(raw), 40)
S_raw <- shap_mat(rf_raw, raw[ex_idx, drv], raw[bg_idx, drv])
grp_ex <- raw[[UNIT]][ex_idx]
shap_icc <- sapply(drv, function(f) between_share(S_raw[[f]], grp_ex))
mean_shap_icc <- mean(shap_icc, na.rm=TRUE)
imp_raw <- pct(colMeans(abs(S_raw)))
panel_tab <- data.frame(feature = drv, shap_importance_pct = unname(imp_raw[drv]),
                        shap_icc = round(unname(shap_icc[drv]), 3)) %>%
  arrange(desc(shap_importance_pct))

# (2b) separated within representation: SHAP recovers the within drivers
wn <- paste0(drv,"_wn")
dw <- panel %>% select(all_of(UNIT), all_of(OUTCOME), all_of(wn)) %>%
  filter(if_all(everything(), ~!is.na(.x))) %>%
  group_by(.data[[UNIT]]) %>% mutate(y_wn = .data[[OUTCOME]] - mean(.data[[OUTCOME]])) %>%
  ungroup() %>% as.data.frame()
keep_wn <- wn[sapply(dw[,wn], function(z) sd(z) > 1e-8)]
rf_wn <- randomForest(reformulate(keep_wn,"y_wn"), dw[,c("y_wn",keep_wn)], ntree=400)
set.seed(7)
exw_idx <- sample(nrow(dw), min(500, nrow(dw))); bgw_idx <- sample(nrow(dw), 40)
S_wn <- shap_mat(rf_wn, dw[exw_idx, keep_wn], dw[bgw_idx, keep_wn])
imp_wn <- sort(pct(colMeans(abs(S_wn))), decreasing=TRUE)
within_tab <- data.frame(feature = names(imp_wn), shap_importance_pct = unname(imp_wn))

out_panel <- bind_rows(
  panel_tab %>% transmute(stage="raw pooled", feature, shap_importance_pct, shap_icc),
  within_tab %>% transmute(stage="within representation", feature, shap_importance_pct, shap_icc=NA_real_))
write.csv(out_panel, file.path(OUT,"shap_panel.csv"), row.names=FALSE)
write.csv(data.frame(icc_y=round(icc_y,3), mean_shap_icc_raw=round(mean_shap_icc,3)),
          file.path(OUT,"shap_panel_summary.csv"), row.names=FALSE)

# Figure: 2 panels
have_patch <- requireNamespace("patchwork", quietly=TRUE)
th <- theme_minimal(base_size=12) +
  theme(plot.title=element_text(face="bold", size=11), panel.grid.minor=element_blank())
pA_df <- panel_tab %>% mutate(feature=factor(feature, levels=feature[order(shap_icc)]))
pA <- ggplot(pA_df, aes(feature, shap_icc, fill=shap_icc)) +
  geom_col(width=.66) +
  geom_hline(yintercept=icc_y, linetype="dashed", colour="grey30") +
  annotate("text", x=nrow(pA_df)-0.3, y=icc_y-0.02,
           label=sprintf("outcome ICC = %.2f", icc_y), hjust=1, vjust=0.5, size=3.0, colour="grey30") +
  coord_flip() + ylim(0,1) +
  scale_fill_gradient(low="#7fb3d5", high="#1b4f72", guide="none") + th +
  labs(x=NULL, y="ICC of SHAP attribution")
pB_df <- sim_tab %>% mutate(representation=factor(representation, levels=representation),
                            lab=sprintf("%.1f%%", within_share_pct))
pB <- ggplot(pB_df, aes(representation, within_share_pct, fill=recovered)) +
  geom_col(width=.6) + geom_text(aes(label=lab), vjust=-0.4, size=3.4) +
  scale_fill_manual(values=c(no="#d95f02", yes="#1b9e77"), name="within driver\nrecovered") +
  scale_x_discrete(labels=function(x) gsub(", ", ",\n", gsub(" \\(", "\n(", x))) +
  th + theme(axis.text.x=element_text(size=8)) + ylim(0, max(pB_df$within_share_pct)*1.25) +
  labs(x=NULL, y="SHAP share of true within driver (%)")
if (have_patch){ library(patchwork); p <- pA + pB + patchwork::plot_layout(widths=c(1,1))
} else { p <- gridExtra::arrangeGrob(pA, pB, ncol=2) }
ggsave(file.path(PLT,"shap_panel.png"), p, width=11, height=4.4, dpi=140)
ggsave(file.path(FEN,"shap_panel.png"), p, width=11, height=4.4, dpi=140)
