# ---- Internal: histogram-based dissimilarity for one TF -----------------

#' Compute MinMax dissimilarity score for a single TF across labels
#'
#' For each label, a density histogram of AUC values is computed over
#' \code{n_breaks} equal-width bins on \code{[0, 1]}. The dissimilarity is
#' defined as \code{sum(|max - min|) / n_breaks / n_non_zero_rows}, where
#' \code{max} and \code{min} are taken column-wise across labels and
#' \code{n_non_zero_rows} is the number of labels with at least one non-zero
#' bin.
#'
#' @param values_list named list of numeric vectors (one per label) containing
#'   the AUC scores for this TF.
#' @param n_breaks integer number of histogram bins (default 100).
#' @return numeric scalar dissimilarity score.
#' @keywords internal
.minmax_dissimilarity <- function(values_list, n_breaks = 100L) {

  breaks <- c(seq(0, 1, 1 / n_breaks))

  # Build matrix: rows = labels, cols = bins
  hist_mat <- do.call(rbind, lapply(values_list, function(vals) {
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0L) return(rep(NA_real_, n_breaks))
    h <- graphics::hist(vals, breaks = breaks, plot = FALSE)
    h$density
  }))

  # Keep only complete rows (labels with no NA bins)
  hist_mat <- hist_mat[stats::complete.cases(hist_mat), , drop = FALSE]

  if (nrow(hist_mat) == 0L || ncol(hist_mat) == 0L) return(0)

  minmax_diff <- apply(stats::na.omit(hist_mat), 2, max) -
                 apply(stats::na.omit(hist_mat), 2, min)
  variant <- sum(abs(minmax_diff)) / n_breaks

  # Normalize by number of non-zero rows (labels with signal)
  n_non_zero <- sum(rowSums(hist_mat) != 0)
  if (n_non_zero > 0) {
    variant <- variant / n_non_zero
  }

  variant
}

# ---- Public: calculate_dissimilarity ------------------------------------

#' Calculate MinMax dissimilarity scores across phenotype labels
#'
#' For every TF in the AUC data, computes a histogram-based dissimilarity
#' score that quantifies how differently the TF activity is distributed
#' across the specified phenotype labels.
#'
#' This mirrors \code{SimiCPipeline.calculate_dissimilarity} from the
#' Python package.
#'
#' @param x \code{SimiCvizExperiment} object with AUC data.
#' @param tf_names character vector of TFs to evaluate (default: all).
#' @param labels integer vector of labels (or character label names) to
#'   compare (default: all). Must contain at least two labels.
#' @param n_breaks integer number of histogram bins (default 100).
#' @param verbose logical; if \code{TRUE} (default), print progress info.
#'
#' @return A \code{data.frame} with columns \code{TF} and
#'   \code{MinMax_score}, sorted by descending dissimilarity score.
#'   Row names are set to TF names.
#' @export
calculate_dissimilarity <- function(x,
                                    tf_names = NULL,
                                    labels = NULL,
                                    n_breaks = 100L,
                                    verbose = TRUE) {
  if (!is.SimiCvizExperiment(x)) stop("x must be a SimiCvizExperiment")

  labels   <- .resolve_labels(x, labels)
  tf_names <- .resolve_tf_names(x, tf_names)

  if (length(labels) < 2L) {
    stop("At least two labels are required to calculate dissimilarity scores.")
  }

  if (verbose) {
    message("Calculating dissimilarity scores for ", length(tf_names),
            " TFs across ", length(labels), " labels...")
  }

  # Pre-extract AUC subsets per label (avoids repeated subsetting)
  auc_per_label <- stats::setNames(
    lapply(labels, function(lab) {
      tryCatch(.subset_auc_for_label(x, lab), error = function(e) NULL)
    }),
    as.character(labels)
  )

  scores <- vapply(tf_names, function(tf) {
    vals_list <- lapply(auc_per_label, function(sub) {
      if (is.null(sub) || !tf %in% colnames(sub)) return(numeric(0))
      stats::na.omit(sub[[tf]])
    })
    # Only keep labels that have data
    vals_list <- vals_list[vapply(vals_list, length, integer(1)) > 0L]
    if (length(vals_list) < 2L) return(NA_real_)
    .minmax_dissimilarity(vals_list, n_breaks = n_breaks)
  }, numeric(1))

  df <- data.frame(
    TF = tf_names,
    MinMax_score = scores,
    stringsAsFactors = FALSE
  )
  rownames(df) <- df$TF

  # Sort descending

  df <- df[order(df$MinMax_score, decreasing = TRUE), , drop = FALSE]

  if (verbose) {
    n_show <- min(10L, nrow(df))
    message("Top ", n_show, " TFs by MinMax dissimilarity score:")
    top <- utils::head(df, n_show)
    for (i in seq_len(nrow(top))) {
      message(sprintf("  %s: %.4f", top$TF[i], top$MinMax_score[i]))
    }
  }

  df
}
