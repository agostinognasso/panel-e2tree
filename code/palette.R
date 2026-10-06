# palette.R — the one place where the manuscript figures get their colours.
# Sourced by 03, 04, 07 and 18, so the seven figures cannot drift apart.
#
# Two hues, chosen as a warm/cool pair that reads as opposite:
#   teal       #0f8fb5   the cool pole, and the single-series accent
#   terracotta #d1552b   the warm pole
# Checked for colour-vision deficiency before adoption: worst adjacent pair
# dE 19.6 under protanopia and 27.9 under normal vision (OKLab x100, against
# thresholds of 8 and 15), and both clear 3:1 contrast on a white surface.
#
# Which ramp goes where follows the data, not taste:
#   signed values (within deviations, R2 above/below zero) -> DIVERGING,
#     terracotta for negative, teal for positive, a neutral grey at zero;
#   one-directional magnitude (firm-mean solvency, ICC of an attribution)
#     -> SEQUENTIAL, the teal ramp alone, light to dark;
#   a single series -> the teal accent, no ramp and no legend.
#
# The box ramps stop well short of the dark end on purpose: rpart.plot writes
# the node text in black over the fill, and it has to stay readable.

PAL <- list(
  teal        = "#0f8fb5",
  teal_dark   = "#0a6b8a",
  terracotta  = "#d1552b",
  terra_dark  = "#a84120",
  neutral     = "#f0efec",   # the zero midpoint of the diverging ramps
  grid        = "grey55",
  ink         = "grey20"
)

# Sequential teal, light to dark. Lightness is monotone by construction.
PAL$seq_teal <- c("#eaf5f9", "#d3ebf3", "#b7dfec", "#97d1e4", "#74c2db",
                  "#4bb0d0", "#2e9fc4", "#0f8fb5")

# The light ends of each hue, for tree boxes that carry black text.
PAL$box_teal  <- c("#eaf5f9", "#d3ebf3", "#b7dfec", "#97d1e4", "#74c2db")
PAL$box_terra <- c("#fbede7", "#f7dace", "#f2c5b2", "#ecae93", "#e59572")

# rpart.plot takes a vector of box colours but not a pair of custom ramps (its
# two-sided palettes only work with its own named ones), and mapping a single
# vector across the value range would put the neutral at the midpoint of that
# range rather than at zero. So the diverging boxes are computed per node:
# zero gets the neutral, and each arm is scaled by the largest absolute value,
# which keeps the two signs comparable.
diverging_box_col <- function(y, lo = PAL$box_terra, hi = PAL$box_teal,
                              mid = PAL$neutral) {
  m <- suppressWarnings(max(abs(y), na.rm = TRUE))
  if (!is.finite(m) || m == 0) return(rep(mid, length(y)))
  ramp_lo <- grDevices::colorRamp(c(mid, lo))
  ramp_hi <- grDevices::colorRamp(c(mid, hi))
  vapply(y, function(v) {
    if (is.na(v)) return(mid)
    z <- ramp_hi(abs(v) / m)
    if (v < 0) z <- ramp_lo(abs(v) / m)
    grDevices::rgb(z[1], z[2], z[3], maxColorValue = 255)
  }, character(1))
}
