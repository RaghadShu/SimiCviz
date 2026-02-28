#' Plot adjusted R² distributions
#'
#' Plots histograms of adjusted R² values per label, similar to the Python
#' \code{SimiCVisualization$plot_r2_distribution}.
#'
#' @param adjusted_r_squared A \strong{named list} of numeric vectors, one per
#'   label.
#'   Names must be label identifiers (e.g. \code{"0"}, \code{"1"}, …).
#'   This metric is specific to SimiC's regression step; other methods may not
#'   produce it, so it is kept as an explicit input rather than extracted from
#'   the experiment object.
#' @param x Optional \code{SimiCvizExperiment} object used solely to resolve
#'   display names and colors for labels.
#'   If \code{NULL}, default \code{"Label i"} naming and a built-in palette are
#'   used.
#' @param labels Optional vector of labels to plot (subset of
#'   \code{names(adjusted_r_squared)}). Defaults to all.
#' @param threshold Numeric R² threshold line and summary-statistic cutoff
#'   (default \code{0.7}).
#' @param nrow,ncol Optional grid layout. If both \code{NULL}, one label per
#'   row.
#'
#' @return Called for side effects (plots). Returns \code{invisible(NULL)}.
#' @export
plot_r2_distribution <- function(adjusted_r_squared,
                                 x = NULL,
                                 labels = NULL,
                                 threshold = 0.7,
                                 grid = NULL) {

  # --- validate adjusted_r_squared ---
  if (missing(adjusted_r_squared) || is.null(adjusted_r_squared) || !is.list(adjusted_r_squared)) {
    stop("`adjusted_r_squared` must be a named list of numeric vectors (one per label).")
  }

  obj <- adjusted_r_squared

  all_labels <- names(obj)
  if (is.null(all_labels) || any(all_labels == "")) {
    stop("`adjusted_r_squared` must be a *named* list with label identifiers as names.")
  }

  # --- optional experiment object for display names / colors ---
  sim_obj <- NULL
  if (!is.null(x)) {
    if (!is.SimiCvizExperiment(x)) {
      stop("`x` must be a SimiCvizExperiment or NULL.")
    }
    sim_obj <- x
  }

  # --- label subsetting ---
  if (!is.null(labels)) {
    lab_chr <- as.character(labels)
    missing_labs <- setdiff(lab_chr, all_labels)
    if (length(missing_labs)) {
      warning("Some requested labels not found and will be ignored: ",
              paste(missing_labs, collapse = ", "))
    }
    use_labels <- intersect(lab_chr, all_labels)
    if (!length(use_labels)) {
      stop("None of the requested labels are present in adjusted_r_squared.")
    }
  } else {
    use_labels <- all_labels
  }

  n_labels <- length(use_labels)

  # --- grid layout ---
  if (is.null(grid)){
    nrow <- n_labels
    ncol <- 1L
  } else if (length(grid) != 2 || !is.numeric(grid) || any(grid < 0)) {
    stop("`grid` must be a numeric positive vector of length 2 (nrow, ncol).")
  } else{
    # take the min number between the provided grid and the number of labels
    nrow <- ifelse(as.integer(grid[1]) == 0, n_labels, min(as.integer(grid[1]), n_labels)) 
    ncol <- ifelse(as.integer(grid[2]) == 0, ceiling(n_labels / nrow), min(as.integer(grid[2]), ceiling(n_labels / nrow))) 
  }

  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(mfrow = c(nrow, ncol))

  # --- iterate over labels ---
  for (i in seq_len(n_labels)) {
    lab_name <- use_labels[[i]]
    r2_vals  <- obj[[lab_name]]

    if (is.null(r2_vals) || !length(r2_vals)) {
      plot.new()
      title(main = paste("No data for label", lab_name))
      next
    }

    sel        <- r2_vals > threshold
    n_selected <- sum(sel, na.rm = TRUE)
    mean_r2    <- if (n_selected > 0) mean(r2_vals[sel], na.rm = TRUE) else 0

    lab_id_num <- tryCatch(as.integer(lab_name),
                           warning = function(.) NA_integer_,
                           error   = function(.) NA_integer_)
    if (is.na(lab_id_num)) lab_id_num <- i

    lab_display <- .get_label_name(sim_obj, lab_id_num)
    lab_color   <- .get_label_color(sim_obj, lab_id_num, idx = i)

    n_tfs <- max(100, length(x@tf_ids))
    main_title <- sprintf(
      "%s\nTargets selected: %d, Mean R\u00b2: %.3f",
      lab_display, n_selected, mean_r2
    )

    r2_histogram(
      adjusted_r_squared = r2_vals,
      threshold          = threshold,
      n_tfs              = n_tfs,
      main               = main_title,
      col                = lab_color
    )
  }
  invisible(NULL)
}

# internal helper for a single histogram, similar to Python implementation
r2_histogram <- function(adjusted_r_squared,
                         threshold = 0.7,
                         n_tfs = 100,
                         main = "",
                         col = "grey") {
  hist(
    adjusted_r_squared,
    col      = col,
    breaks   = n_tfs,
    xlab     = "Adjusted R²",
    main     = main,
    border   = "black"
  )
  abline(v = threshold, col = "red", lwd = 2, lty = 2)
  grid(nx = NA, ny = NULL, col = "grey80", lty = "dotted")
}
