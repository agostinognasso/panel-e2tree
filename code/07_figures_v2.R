# 07_figures_v2.R — regenerates the two figures the v2 revision re-captioned.
# No refitting: reads output/repro/necessity.csv (Fig. 2a), output/budget59/
# budget59_grid.csv (the common-budget row) and output/shap_sim.csv (Fig. 3b).
#
# Fig. 2a  manuscript/manuscript_DSS_v2/figures/e2tree_necessity.png
#          - approach labels carry the actual leaf count
#          - y axis "Outcome recovery R^2 (1 - SSE/SST)"
#          - the common-budget pooled configuration (a'') is included
# Fig. 3b  manuscript/manuscript_DSS_v2/figures/shap_panel_B.png
#          - no recovered no/yes legend, single colour
#          - y axis "SHAP importance share of within driver B (%)"
#
# The v1 figures under manuscript/figures/ are left untouched.
# Usage: Rscript code/07_figures_v2.R

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })

BASE <- normalizePath(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))), ".."))
EN <- file.path(BASE, "manuscript", "manuscript_DSS_v2", "figures")
source(file.path(BASE, "code", "palette.R"))
stopifnot(dir.exists(EN))

th <- theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 11),
        panel.grid.minor = element_blank())

# --------------------------------------------------------------------------
# Fig. 2a — outcome recovery by explanation strategy, at actual leaf counts
# --------------------------------------------------------------------------
nec <- read.csv(file.path(BASE, "output", "repro", "necessity.csv"))
grd <- read.csv(file.path(BASE, "output", "budget59", "budget59_grid.csv"))

BUDGET <- nec$leaves[grep("^\\(d\\)", nec$approach)]          # 59
# Selection rule, stated in the paper: among the stopping-rule configurations
# yielding at most `BUDGET` terminal nodes, take the one with the most terminal
# nodes; ties are broken towards the higher total recovery, i.e. in favour of the
# pooled competitor.
cand <- grd[grd$tag == "budget sweep" & grd$leaves <= BUDGET, ]
sel  <- cand[order(-cand$leaves, -cand$r2_vs_y), ][1, ]
cat(sprintf("common-budget row: %d leaves, depth %d, level=%g n=%g -> r2 %.3f, within %+.3f, between %.3f\n",
            sel$leaves, sel$depth, sel$level, sel$n, sel$r2_vs_y, sel$r2_within, sel$r2_between))

nec2 <- rbind(
  nec[1:2, ],
  data.frame(approach = sprintf("(a\u2033) pooled, %d-leaf budget", BUDGET),
             R2_vs_y = round(sel$r2_vs_y, 3), R2_within = round(sel$r2_within, 3),
             R2_between = round(sel$r2_between, 3), leaves = sel$leaves),
  nec[3:4, ])
nec2$approach[1] <- "(a) pooled e2tree, raw features"
nec2$approach[2] <- "(a\u2032) pooled, more permissive configuration"
nec2$label <- sprintf("%s (%d leaves)", nec2$approach, nec2$leaves)

nec_long <- nec2 %>%
  mutate(label = factor(label, levels = rev(label))) %>%
  pivot_longer(c(R2_vs_y, R2_within), names_to = "metric", values_to = "r2") %>%
  mutate(metric = factor(ifelse(metric == "R2_vs_y", "total", "within-unit component"),
                         levels = c("total", "within-unit component")))

p_nec <- ggplot(nec_long, aes(label, r2, fill = r2 > 0)) +
  geom_col(width = .65) +
  geom_hline(yintercept = 0, colour = PAL$grid, linewidth = .3) +
  geom_text(aes(label = sprintf("%.2f", r2), hjust = ifelse(r2 >= 0, -0.15, 1.1)), size = 3.2) +
  coord_flip() + facet_wrap(~metric) +
  # sign is already carried by the bar direction and the printed value, so the
  # fill is reinforcement, not the only cue: no legend.
  scale_fill_manual(values = c(`TRUE` = PAL$teal, `FALSE` = PAL$terracotta), guide = "none") +
  ylim(-0.6, 1) + th +
  labs(x = NULL, y = expression("Outcome recovery " * R^2 * " (1 " - " SSE/SST)"))
ggsave(file.path(EN, "e2tree_necessity.png"), p_nec, width = 10.4, height = 4.2, dpi = 140)
cat("written: ", file.path(EN, "e2tree_necessity.png"), "\n")

# --------------------------------------------------------------------------
# Fig. 3b — SHAP importance share of the true within driver
# --------------------------------------------------------------------------
sim_tab <- read.csv(file.path(BASE, "output", "shap_sim.csv"))
pB_df <- sim_tab %>%
  mutate(representation = factor(representation, levels = representation),
         lab = sprintf("%.1f%%", within_share_pct))
pB <- ggplot(pB_df, aes(representation, within_share_pct)) +
  geom_col(width = .6, fill = PAL$teal) +
  geom_text(aes(label = lab), vjust = -0.4, size = 4) +
  scale_x_discrete(labels = function(x) gsub(", ", ",\n", gsub(" \\(", "\n(", x))) +
  ylim(0, max(pB_df$within_share_pct) * 1.25) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(size = 9)) +
  labs(x = NULL, y = "SHAP importance share of within driver B (%)")
ggsave(file.path(EN, "shap_panel_B.png"), pB, width = 6.4, height = 4.4, dpi = 150)
cat("written: ", file.path(EN, "shap_panel_B.png"), "\n")
