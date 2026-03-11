#' Compute MinMax dissimilarity scores (internal)
#'
#' @param auc_list Named list of data.frames (one per label).
#' @param tf_names Character vector of TF column names.
#' @param labels Vector of label keys to compare.
#' @param n_breaks Integer number of histogram bins (default 100).
#' @return Numeric vector of dissimilarity scores.
#' @keywords internal
#' @noRd
.compute_dissimilarity_scores <- function(auc_list, tf_names, labels,
                                          n_breaks = 100L) {
  breaks <- seq(0, 1, length.out = n_breaks + 1L)
  scores <- numeric(length(tf_names))

  for (i in seq_along(tf_names)) {
    tf <- tf_names[i]
    hist_mat <- list()

    for (lab in labels) {
      df <- auc_list[[as.character(lab)]]
      if (is.null(df) || !(tf %in% colnames(df))) next
      vals <- df[[tf]]
      vals <- vals[!is.na(vals)]
      if (length(vals) == 0L) next
      h <- graphics::hist(vals, breaks = breaks, plot = FALSE)
      # density-normalised counts (sum * bin_width ≈ 1)
      hist_mat[[as.character(lab)]] <- h$density
    }

    if (length(hist_mat) < 2L) {
      scores[i] <- 0
      next
    }

    mat <- do.call(rbind, hist_mat)              # labels × bins
    # Remove columns that are all NA
    keep <- apply(mat, 2, function(col) !all(is.na(col)))
    mat <- mat[, keep, drop = FALSE]

    if (ncol(mat) == 0L) {
      scores[i] <- 0
      next
    }

    minmax_diff <- apply(mat, 2, max, na.rm = TRUE) -
      apply(mat, 2, min, na.rm = TRUE)
    variant <- sum(abs(minmax_diff)) / n_breaks
    non_zero_rows <- sum(rowSums(mat, na.rm = TRUE) != 0)
    if (non_zero_rows > 0) variant <- variant / non_zero_rows

    scores[i] <- variant
  }

  scores
}


#' Resolve label arguments to integer label keys
#'
#' Accepts integer keys, character label names, or a mix, and returns
#' the corresponding integer label keys from the experiment object.
#'
#' @param x A \code{SimiCvizExperiment} object.
#' @param labels Vector of integer keys or character label names, or NULL (all).
#' @return Integer vector of resolved label keys.
#' @keywords internal
#' @noRd
.resolve_labels_dissim <- function(x, labels) {
  cl <- x@cell_labels
  all_labels <- sort(unique(cl$label))

  if (is.null(labels)) return(all_labels)

  # Map label_names back to integer keys
  name_to_key <- stats::setNames(
    as.integer(names(x@label_names)),
    x@label_names
  )

  resolved <- vapply(labels, function(lab) {
    # Already an integer key
    if (is.numeric(lab) || is.integer(lab)) {
      return(as.integer(lab))
    }
    # Try matching as a label name
    if (is.character(lab) && lab %in% names(name_to_key)) {
      return(name_to_key[[lab]])
    }
    # Try coercing string to integer (e.g. "0", "1")
    as_int <- suppressWarnings(as.integer(lab))
    if (!is.na(as_int) && as_int %in% all_labels) {
      return(as_int)
    }
    stop(sprintf(
      "Label '%s' not found. Valid keys: [%s]. Valid names: [%s].",
      lab,
      paste(all_labels, collapse = ", "),
      paste(x@label_names, collapse = ", ")
    ))
  }, integer(1), USE.NAMES = FALSE)

  unique(resolved)
}


#' Calculate dissimilarity scores across phenotype labels
#'
#' Compares per-cell TF activity-score distributions between phenotype labels
#' using a histogram-based MinMax divergence.
#'
#' @param x A \code{SimiCvizExperiment} object.
#' @param labels Integer vector of label keys or character vector of label
#'   display names to compare (default: all labels). Mixing types is allowed.
#' @param cell_groups Optional named list mapping group names to character
#'   vectors of cell identifiers.
#' @param n_breaks Number of histogram bins (default 100).
#' @param verbose Logical; print progress messages (default \code{TRUE}).
#'
#' @return A \code{data.frame} with TFs as rows, sorted by score (descending).
#'
#' @rdname calculate_dissimilarity
#' @export
calculate_dissimilarity <- function(x,
                                    labels = NULL,
                                    cell_groups = NULL,
                                    n_breaks = 100L,
                                    verbose = TRUE) {
  if (!is.SimiCvizExperiment(x)) {
    stop("calculate_dissimilarity: 'x' must be a SimiCvizExperiment.")
  }

  # Need collected AUC + cell_labels
  auc_df <- x@auc$collected
  if (is.null(auc_df) || nrow(auc_df) == 0L) {
    stop("No collected AUC data in the experiment object.")
  }
  cl <- x@cell_labels
  if (nrow(cl) == 0L) {
    stop("No cell_labels in the experiment object.")
  }

  # Resolve labels (accepts integers, label names, or mix)
  labels <- .resolve_labels_dissim(x, labels)
  if (length(labels) < 2L) {
    stop("Need at least 2 labels to compute dissimilarity.")
  }

  # TF columns (exclude 'label' if present)
  tf_cols <- setdiff(colnames(auc_df), "label")

  # Build per-label AUC list (subsetting cells)
  .build_auc_list <- function(cell_subset = NULL) {
    out <- list()
    for (lab in labels) {
      cells_in_lab <- cl$cell[cl$label == lab]
      if (!is.null(cell_subset)) {
        cells_in_lab <- intersect(cells_in_lab, cell_subset)
      }
      idx <- which(rownames(auc_df) %in% cells_in_lab)
      if (length(idx) == 0L) next
      out[[as.character(lab)]] <- auc_df[idx, tf_cols, drop = FALSE]
    }
    out
  }

  # --- No cell groups: global dissimilarity ---
  if (is.null(cell_groups)) {
    if (verbose) message("Calculating dissimilarity scores (all cells)...")
    auc_list <- .build_auc_list()
    scores <- .compute_dissimilarity_scores(auc_list, tf_cols,
                                            as.character(labels), n_breaks)
    result <- data.frame(MinMax_score = scores, row.names = tf_cols,
                         stringsAsFactors = FALSE)
    result <- result[order(result$MinMax_score, decreasing = TRUE), , drop = FALSE]
    if (verbose) {
      message(sprintf("Top 10 TFs by MinMax dissimilarity:"))
      print(utils::head(result, 10))
    }
    return(result)
  }


  # --- Per-group dissimilarity ---
  if (verbose) {
    message(sprintf("Calculating dissimilarity for %d cell groups...",
                    length(cell_groups)))
  }

  group_scores <- list()
  for (gname in names(cell_groups)) {
    cell_ids <- cell_groups[[gname]]
    if (verbose) message(sprintf("  Group: %s (%d cells)", gname, length(cell_ids)))
    auc_list <- .build_auc_list(cell_subset = cell_ids)
    if (length(auc_list) < 2L) {
      if (verbose) message(sprintf("    Warning: <2 labels with cells for '%s', skipping.", gname))
      group_scores[[gname]] <- rep(0, length(tf_cols))
    } else {
      group_scores[[gname]] <- .compute_dissimilarity_scores(
        auc_list, tf_cols, names(auc_list), n_breaks
      )
    }
  }

  result <- as.data.frame(group_scores, row.names = tf_cols,
                          stringsAsFactors = FALSE)
  result$mean_score <- rowMeans(result, na.rm = TRUE)
  result <- result[order(result$mean_score, decreasing = TRUE), , drop = FALSE]

  if (verbose) {
    message("Top 10 TFs by mean dissimilarity across groups:")
    print(utils::head(result, 10))
  }

  result
}
