# Publication SHAP figure for the ORBIS solvency application: the two panels the
# manuscript uses as separate files (17_shap_panel.R draws the same content only
# as a single combined diagnostic in output/plots/). Reads the saved CSV outputs
# of 17_shap_panel.R and re-renders with the manuscript labels; no refitting.
# Output: manuscript/figures/{shap_panel_A,shap_panel_B}.png
suppressPackageStartupMessages({library(dplyr); library(ggplot2)})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
ROOT <- if (length(script_path)) dirname(dirname(normalizePath(script_path))) else normalizePath("..")
DATA <- file.path(ROOT, "output")
EN   <- file.path(ROOT, "manuscript", "manuscript_DSS_v2", "figures")
source(file.path(ROOT, "code", "palette.R"))
dir.create(EN, recursive = TRUE, showWarnings = FALSE)

panel_tab <- read.csv(file.path(DATA,"shap_panel.csv")) %>%
  filter(stage=="raw pooled") %>% filter(!is.na(shap_icc))
icc_y  <- read.csv(file.path(DATA,"shap_panel_summary.csv"))$icc_y
sim_tab <- read.csv(file.path(DATA,"shap_sim.csv"))

th <- theme_minimal(base_size=13) +
  theme(plot.title=element_text(face="bold", size=12),
        plot.subtitle=element_text(size=11), panel.grid.minor=element_blank())

# Panel A
pA_df <- panel_tab %>% mutate(feature=factor(feature, levels=feature[order(shap_icc)]))
# One flat fill, not a ramp: every ICC here sits between 0.86 and 0.96, so a
# gradient over that range would fade the lowest bar to near-nothing and make a
# high ICC read as a low one. The bar length already carries the value.
pA <- ggplot(pA_df, aes(feature, shap_icc)) +
  geom_col(width=.66, fill=PAL$teal) +
  geom_hline(yintercept=icc_y, linetype="dashed", colour=PAL$ink) +
  annotate("label", x=1.6, y=0.99, label=sprintf("outcome ICC = %.2f", icc_y),
           hjust=1, vjust=0.5, size=3.4, colour="grey20",
           fill="white", label.size=0, label.r=unit(0.1,"lines")) +
  coord_flip() + ylim(0,1) +
  th +
  labs(x=NULL, y="ICC of SHAP attribution")
ggsave(file.path(EN,"shap_panel_A.png"), pA, width=7.2, height=4.4, dpi=150)

# Panel B
pB_df <- sim_tab %>%
  mutate(representation=factor(representation, levels=representation),
         lab=sprintf("%.1f%%", within_share_pct))
pB <- ggplot(pB_df, aes(representation, within_share_pct, fill=recovered)) +
  geom_col(width=.6) + geom_text(aes(label=lab), vjust=-0.4, size=4) +
  scale_fill_manual(values=c(no=PAL$terracotta, yes=PAL$teal), name="within driver\nrecovered") +
  scale_x_discrete(labels=function(x) gsub(", ", ",\n", gsub(" \\(", "\n(", x))) +
  th + theme(axis.text.x=element_text(size=9)) +
  ylim(0, max(pB_df$within_share_pct)*1.25) +
  labs(x=NULL, y="SHAP share of true within driver (%)")
ggsave(file.path(EN,"shap_panel_B.png"), pB, width=6.4, height=4.4, dpi=150)
