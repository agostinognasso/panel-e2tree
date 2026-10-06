# Draw the e2tree surrogate diagrams (between / within) for the ORBIS solvency
# application. These figures are produced here only: 03_panel_e2tree_orbis.R fits
# and caches the panel model but does not draw the trees. We read the cached
# object S1_main.rds (the same e2panel object behind the paper's numbers, so the
# trees match the reported between/within splits exactly) and render the plots
# without refitting -- keeping figure styling decoupled from the costly fit.
# Output: manuscript/figures/{e2tree_diagram,panele2_within_diagram,panele2_within_diagram_top}.png
suppressPackageStartupMessages({library(e2tree); library(rpart.plot)})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
ROOT <- if (length(script_path)) dirname(dirname(normalizePath(script_path))) else normalizePath("..")
EN <- file.path(ROOT, "manuscript", "manuscript_DSS_v2", "figures")
source(file.path(ROOT, "code", "palette.R"))
dir.create(EN, recursive = TRUE, showWarnings = FALSE)
S1 <- readRDS(file.path(ROOT, "output", "repro", "cache", "S1_main.rds"))

# Display-only truncation to the top `d` levels of splits: internal nodes at
# depth d become leaves showing their node mean, deeper branches are dropped.
# Node numbers follow the binary-heap convention, so depth = floor(log2(node)),
# and (with compete/surrogate splits = 0) split rows map 1:1 to internal nodes
# in frame order. Metrics in the paper are computed on the FULL tree, never on
# this display copy. (rpart::snip.rpart hangs on the e2tree object, hence the
# manual surgery.)
truncate_rpart <- function(rp, d){
  fr  <- rp$frame
  dep <- floor(log2(as.integer(rownames(fr))))
  internal      <- fr$var != "<leaf>"
  becomes_leaf  <- dep == d & internal
  keep          <- dep <= d
  keep_internal <- internal & !becomes_leaf & keep
  rp$splits <- rp$splits[keep_internal[internal], , drop=FALSE]
  fr$var[becomes_leaf] <- "<leaf>"                 # var is character; keep it so
  fr$ncompete[becomes_leaf] <- 0L; fr$nsurrogate[becomes_leaf] <- 0L
  rp$frame <- fr[keep, , drop=FALSE]
  rp
}

# box.palette: one teal ramp for the between tree, whose node values are solvency
# LEVELS and run one way; a two-sided ramp for the within tree, whose node values
# are demeaned deviations, so zero has to read as neutral and the two signs as
# opposite. Both ramps stay light: rpart.plot prints the node text in black.
draw_tree <- function(comp, main, file, w, h, tweak,
                      box.palette = PAL$box_teal, diverging = FALSE,
                      extra=101, split.cex=0.95, space=0.4, mtop=2.2, max_depth=Inf){
  rp <- plot_e2tree(comp$tree, comp$ensemble)
  rp$frame$yval[abs(rp$frame$yval) < 1e-9] <- 0   # clean 0 for demeaned within
  if (is.finite(max_depth)) rp <- truncate_rpart(rp, max_depth)
  # Pass either box.palette or box.col, never both, since rpart.plot derives the
  # one from the other. Hence do.call with only the argument that applies.
  # (The full within tree warns "labs do not fit even at cex 0.15": it is dense
  # enough that the split labels overlap at this size, and it did so with the
  # earlier grey palette too. The main text uses the four-level version.)
  fill <- if (diverging) list(box.col = diverging_box_col(rp$frame$yval))
          else list(box.palette = box.palette)
  png(file.path(EN, file), width=w, height=h, res=150)
  par(mar=c(0,0,0.4,0), xpd=NA)          # no title: only a thin top margin
  do.call(rpart.plot::rpart.plot, c(list(
    rp, type=2, extra=extra, branch=.3,
    fallen.leaves=TRUE, roundint=FALSE, tweak=tweak,
    split.cex=split.cex, faclen=0, varlen=0, gap=0, space=space,
    shadow.col=NULL, main=NULL),      # figure titles removed (kept in captions)
    fill))
  dev.off()
}

# Between tree: firm means, what distinguishes structurally well- vs under-capitalised firms.
draw_tree(S1$between,
          "Between tree: what distinguishes firms' structural solvency (firm means)",
          "e2tree_diagram.png", w=2600, h=1500, tweak=1.05)

# Within tree: demeaned deviations, what moves a firm's solvency over time.
# Full tree: supplementary material (and the object all reported metrics use).
draw_tree(S1$within,
          "Within tree: what moves solvency over time (within-firm deviations)",
          "panele2_within_diagram.png", w=3800, h=1800, tweak=0.92,
          diverging=TRUE,
          extra=1, split.cex=0.9, space=0.6, mtop=1.6)

# Main-text display: top four levels of the same within tree, for readability.
draw_tree(S1$within,
          "Within tree, top four levels: what moves solvency over time (within-firm deviations)",
          "panele2_within_diagram_top.png", w=2800, h=1080, tweak=1.05,
          diverging=TRUE,
          extra=1, split.cex=0.95, space=0.5, mtop=1.8, max_depth=4)
