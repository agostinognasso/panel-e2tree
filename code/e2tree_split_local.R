# ============================================================================
# Self-contained candidate-split enumeration for the S3 "bridge" diagnostic.
#
# Pure-R reimplementation of e2tree's INTERNAL split()/ordSplit()/catSplit()
# (dev e2tree 1.2.0.9000). It is bundled here so the bridge diagnostic in
# paper_reproduce.R reproduces the SAME candidate-split matrix -- and the same
# reported number (aug_top10_bw = 0.90) -- against the PUBLIC e2tree (CRAN
# 1.2.0), whose internal split() does NOT accept `max_thresholds` and would
# otherwise enumerate every threshold (changing the diagnostic to 0.70).
#
# The wrapper is renamed `e2_candidate_splits()` (NOT `split`) to avoid masking
# base::split. The impurity step in paper_reproduce.R still uses the stable,
# C++-backed e2tree:::eImpurity(), which is identical in the public and dev
# packages.
# ============================================================================

# numeric / ordered: candidate thresholds capped at `max_thresholds` empirical
# quantiles (type = 1) when a variable has more unique values than the cap.
.e2_ordSplit <- function(x, max_thresholds = 256L) {
  values <- sort(unique(x))
  if (length(values) > 1) {
    values <- values[-length(values)]

    if (is.finite(max_thresholds) && length(values) > max_thresholds) {
      qs <- stats::quantile(x,
                            probs = seq_len(max_thresholds) / (max_thresholds + 1),
                            type = 1, names = FALSE)
      qs <- sort(unique(qs))
      qs <- qs[qs < max(x)]          # a threshold at max(x) is a void split
      if (length(qs) > 0) values <- qs
    }

    Sx <- as.data.frame(outer(x, values, "<=") + 0L)
    colnames(Sx) <- paste("<=", values, sep = "")
    return(Sx)
  } else {
    warning("constant variable: no split produced", call. = FALSE)
    return()
  }
}

# categorical: exhaustive binary partitions for k <= max_cat, otherwise a
# frequency-ordered ordinal fallback (k-1 splits).
.e2_catSplit <- function(x, max_cat = 10) {
  x <- as.character(x)
  values <- sort(unique(x))
  k <- length(values)

  if (k < 2) {
    warning("constant variable: no split produced", call. = FALSE)
    return()
  }

  if (k == 2) {
    Sx <- data.frame(rep(0, length(x)))
    Sx[x == values[1], 1] <- 1
    names(Sx) <- paste("%in% c('", values[1], "')", collapse = "", sep = "")
    return(Sx)
  }

  if (k > max_cat) {
    freq_order <- names(sort(table(x), decreasing = TRUE))
    x_rank <- as.numeric(factor(x, levels = freq_order))
    thresholds <- seq_len(k - 1)
    Sx <- as.data.frame(outer(x_rank, thresholds, "<=") + 0L)
    lab <- vapply(thresholds, function(thr) {
      cats_left <- freq_order[1:thr]
      paste("%in% c('", paste(cats_left, collapse = "', '"), "')", sep = "")
    }, character(1))
    colnames(Sx) <- lab
    return(Sx)
  }

  p <- partitions::listParts(k, do.set = FALSE)
  s <- lapply(p, function(x) {
    lab <- as.character(x)
    if (length(lab) == 2) x <- eval(parse(text = lab[1])) else x <- NA
  })
  s <- s[lengths(s) > 1]
  S <- lapply(s, function(l) cbind(x %in% values[l]))
  lab <- unlist(lapply(s, function(x) {
    paste("%in% c('", paste(values[x], collapse = "', '"), "')", sep = "")
  }))
  Sx <- as.data.frame(do.call(cbind, S))
  colnames(Sx) <- lab
  return(Sx)
}

# wrapper mirroring e2tree:::split(): per-column candidate splits merged into a
# single data frame S (one column per candidate, named "<var> <condition>").
e2_candidate_splits <- function(X, max_cat = 10, max_thresholds = 256L) {
  S <- lapply(X, function(x) {
    type <- class(x)
    if (length(type) > 1) {
      type <- if ("ordered" %in% type) "ordered" else type[1]
    }
    if (type %in% c("character", "factor")) {
      .e2_catSplit(x, max_cat = max_cat)
    } else if (type %in% c("numeric", "ordered", "integer")) {
      .e2_ordSplit(x, max_thresholds = max_thresholds)
    }
  })
  lab <- rep(names(X), lengths(S))
  l <- paste(lab, unlist(lapply(S, names)))
  S <- do.call(cbind, S)
  names(S) <- l
  row.names(S) <- row.names(X)
  list(S = S, lab = lab)
}
