# ---- Internal helpers ----------------------------------------------------

#' Subset collected AUC for a specific label
#'
#' @param x SimiCvizExperiment
#' @param label integer label to subset
#' @return data.frame with cells in rows, TFs in columns (no label column)
#' @keywords internal
.subset_auc_for_label <- function(x, label) {
  if (!is.SimiCvizExperiment(x)) stop("x must be a SimiCvizExperiment")

  lab_int <- as.integer(label)
  lab_chr <- as.character(lab_int)
  if (nrow(x@cell_labels) == 0L){
    stop(sprintf("No cell labels found in the experiment; cannot subset AUC for label '%s'.", lab_chr))
  }
  # Priority 1: collected data.frame with info from  x@cell_lables
  auc_df <- tryCatch(auc_df <- x@auc$collected, error = function(e) NULL)
  if (!is.null(auc_df) && "label" %in% colnames(x@cell_labels)) {
    cell_idx <- x@cell_labels$cell[x@cell_labels$label == lab_int]
    sub <- auc_df[cell_idx, , drop = FALSE]
    return(sub)
  }

  # Priority 3: per_label list
  pl <- tryCatch(x@auc$per_label, error = function(e) NULL)

  if (!is.null(pl) && "label" %in% colnames(x@cell_labels)) {
    if (lab_chr %in% names(pl)) {
      auc_tmp <- as.data.frame(pl[[lab_chr]])
      cell_idx <- x@cell_labels$cell[x@cell_labels$label == lab_int]
      return(auc_tmp[cell_idx, , drop = FALSE])
    }
  }

  stop(sprintf("No AUC data found for label '%s'.", lab_chr))
}

#' Get TF names from collected AUC
#' @keywords internal
.auc_tf_names <- function(x) {
  # Try collected first
  auc_df <- x@auc$collected
  if (!is.null(auc_df)) {
    return(setdiff(colnames(auc_df), "label"))
  }
  # Fallback: per_label — union of all column names
  pl <- x@auc$per_label
  if (!is.null(pl) && length(pl) > 0L) {
    return(unique(unlist(lapply(pl, colnames))))
  }
  character()
}

#' Resolve labels to use for plotting
#' @keywords internal
.resolve_labels <- function(x, labels = NULL) {
  if (!is.null(labels)) {
     if (is.character(labels) & !is.null(x@label_names) & all(labels %in% x@label_names)){
      labels <- names(x@label_names)[which(x@label_names %in% labels)]
     } else if (is.character(labels) & !all(labels %in% x@label_names)) {
      show <- labels[!labels %in% x@label_names]
      stop("Some labels not found in SimiCvizExperiment@label_names: \n\t missing -> ", paste(show,collapse = ", "))
      }

  return(as.integer(labels))
  }
  auc_df <- x@auc$collected
  if (!is.null(auc_df) && "label" %in% colnames(auc_df)) {

    return(sort(unique(auc_df$label)))
  }
  if (nrow(x@cell_labels) > 0L) {

    return(sort(unique(x@cell_labels$label)))
  }

  as.integer(names(x@weights))
}

#' Resolve TF names for plotting
#' @keywords internal
.resolve_tf_names <- function(x, tf_names = NULL) {
  all_tfs <- .auc_tf_names(x)
  if (is.null(tf_names)) return(sort(all_tfs))
  tf_names <- as.character(tf_names)
  valid <- intersect(tf_names, all_tfs)
  if (length(valid) == 0L) stop("None of the specified TFs found in AUC data.")
  valid
}

# ---- ECDF metrics (mirrors Python _compute_ecdf_metrics) ----------------

#' Compute ECDF-based metrics for a numeric vector
#'
#' @param values numeric vector of AUC scores
#' @param x_lower lower integration bound (default 0)
#' @param x_upper upper integration bound (default 1)
#' @param percentile cumulative probability threshold (default 0.5)
#' @return named list with ecdf_auc, auc_at_percentile, x_at_percentile
#' @keywords internal
.compute_ecdf_metrics <- function(values,
                                  x_lower = 0,
                                  x_upper = 1,
                                  percentile = 0.5) {
  n <- length(values)
  sorted_vals <- sort(values)

  in_range <- sorted_vals[sorted_vals > x_lower & sorted_vals < x_upper]
  breakpoints <- c(x_lower, in_range, x_upper)

  total_area <- 0
  area_at_p <- NULL
  x_at_p <- NULL

 for (i in seq_len(length(breakpoints) - 1L)) {
    x_left  <- breakpoints[i]
    x_right <- breakpoints[i + 1L]
    width   <- x_right - x_left

    ecdf_left  <- sum(sorted_vals <= x_left) / n
    ecdf_right <- sum(sorted_vals <= x_right) / n

    if (is.null(area_at_p) && ecdf_right >= percentile) {
      if (ecdf_left >= percentile) {
        x_at_p <- x_left
        area_at_p <- total_area
      } else {
        total_area <- total_area + ecdf_left * width
        x_at_p <- x_right
        area_at_p <- total_area
        next
      }
    }
    total_area <- total_area + ecdf_left * width
  }

  if (is.null(area_at_p)) {
    x_at_p <- x_upper
    area_at_p <- total_area
  }

  range_width <- x_upper - x_lower
  if (range_width > 0) total_area <- total_area / range_width

  list(ecdf_auc = total_area,
       auc_at_percentile = area_at_p,
       x_at_percentile = x_at_p)
}

# ---- Public: calculate_ecdf_auc -----------------------------------------

#' Calculate ECDF-based metrics for TFs across labels
#'
#' Computes three metrics per TF per label:
#' \enumerate{
#'   \item \strong{ECDF-AUC}: area under the ECDF over \code{[0, 1]}
#'   \item \strong{AUC50}: area under the ECDF from 0 to the median activity
#'   \item \strong{x\@p50}: median activity score
#' }
#'
#' @param x \code{SimiCvizExperiment}
#' @param tf_names character vector of TFs (default: all)
#' @param labels integer vector of labels (default: all)
#' @param integration_range numeric length-2 vector (default \code{c(0, 1)})
#' @param percentile cumulative probability threshold (default 0.5)
#' @return data.frame with TFs in rows; per-label metric columns plus
#'   delta columns when >1 label.
#' @export
calculate_ecdf_auc <- function(x,
                               tf_names = NULL,
                               labels = NULL,
                               integration_range = c(0, 1),
                               percentile = 0.5) {
  if (!is.SimiCvizExperiment(x)) stop("x must be a SimiCvizExperiment")

  labels   <- .resolve_labels(x, labels)
  tf_names <- .resolve_tf_names(x, tf_names)
  x_lower  <- integration_range[1]
  x_upper  <- integration_range[2]

  rows <- list()
  for (tf in tf_names) {
    row <- list(TF = tf)
    for (lab in labels) {
      lab_name <- simic@label_names[as.character(lab)]
      sub <- tryCatch(.subset_auc_for_label(x, lab), error = function(e) NULL)
      if (is.null(sub) || !tf %in% colnames(sub)) {
        row[[paste0(lab_name, "_ecdf_auc")]]  <- NA_real_
        row[[paste0(lab_name, "_auc50")]]      <- NA_real_
        row[[paste0(lab_name, "_x_at_p50")]]   <- NA_real_
        next
      }
      vals <- stats::na.omit(sub[[tf]])
      if (length(vals) < 1L) {
        row[[paste0(lab_name, "_ecdf_auc")]]  <- NA_real_
        row[[paste0(lab_name, "_auc50")]]      <- NA_real_
        row[[paste0(lab_name, "_x_at_p50")]]   <- NA_real_
        next
      }
      m <- .compute_ecdf_metrics(vals, x_lower, x_upper, percentile)
      row[[paste0(lab_name, "_ecdf_auc")]]  <- m$ecdf_auc
      row[[paste0(lab_name, "_auc50")]]      <- m$auc_at_percentile
      row[[paste0(lab_name, "_x_at_p50")]]   <- m$x_at_percentile
    }
    rows[[length(rows) + 1L]] <- row
  }

  df <- do.call(rbind, lapply(rows, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
  rownames(df) <- df$TF
  df$TF <- NULL

  # Delta columns
  if (length(labels) > 1L) {
    lab_names  <- simic@label_names
    # lab_names <- vapply(labels, function(l) .get_label_name(x, l), character(1))
    for (suffix in c("_ecdf_auc", "_auc50", "_x_at_p50")) {
      cols <- paste0(lab_names, suffix)
      cols <- intersect(cols, colnames(df))
      if (length(cols) > 1L) {
        vals_mat <- df[, cols, drop = FALSE]
        df[[paste0("delta", suffix)]] <- apply(vals_mat, 1, max, na.rm = TRUE) -
          apply(vals_mat, 1, min, na.rm = TRUE)
      }
    }
  }
  df
}

# ---- Public: plot_auc_distributions --------------------------------------

#' Plot AUC density distributions per TF across phenotype labels
#'
#' Produces paginated PDF output (or on-screen) with density plots for each
#' TF, coloured by label. Mirrors
#' \code{SimiCVisualization.plot_auc_distributions} from the Python package.
#'
#' @param x \code{SimiCvizExperiment}
#' @param tf_names character vector of TFs to plot (default: all)
#' @param labels integer vector of labels  or character label_names (default: all)
#' @param fill logical; fill density curves (default TRUE)
#' @param alpha numeric transparency (default 0.5)
#' @param bw_adjust bandwidth adjustment passed to \code{\link[stats]{density}}
#'   (default 1). Lower = more detail, higher = smoother.
#' @param rug logical; add rug marks (default FALSE)
#' @param grid tuple-like integer vector \code{c(nrow, ncol)} per page.
#'   \code{NULL} puts everything on one page.
#' @param save logical; save to PDF
#' @param filename PDF filename (default auto-generated)
#' @param out_dir output directory for PDF (default \code{getwd()})
#' @param width,height page dimensions in inches
#'
#' @return Invisibly, a list of ggplot objects (one per page).
#' @export
plot_auc_distributions <- function(x,
                                   tf_names = NULL,
                                   labels = NULL,
                                   fill = TRUE,
                                   alpha = 0.5,
                                   bw_adjust = 1,
                                   rug = FALSE,
                                   grid = c(4L, 2L),
                                   save = FALSE,
                                   filename = NULL,
                                   out_dir = getwd(),
                                   width = 14,
                                   height = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
  if (!requireNamespace("gridExtra", quietly = TRUE)) stop("gridExtra required")

  labels   <- .resolve_labels(x, labels)
  lab_names <- x@label_names
  tf_names <- .resolve_tf_names(x, tf_names)
  n_tfs    <- length(tf_names)
  col_map  <-  x@colors
  col_map_named <- stats::setNames(unname(col_map), lab_names[names(col_map)])
  diss_score <- calculate_dissimilarity(x,tf_names = tf_names, labels = labels, verbose = FALSE)
  # Grid setup
  if (is.null(grid)) {
    grid_cols <- 2L
    grid_rows <- ceiling(n_tfs / grid_cols)
    plots_per_page <- n_tfs
  } else {
    grid_rows <- as.integer(grid[1])
    grid_cols <- as.integer(grid[2])
    plots_per_page <- grid_rows * grid_cols
  }
  if (is.null(height)) height <- 5 * grid_rows

  # Build long data per TF
  plot_list <- list()
  for (tf in tf_names) {
    long <- .build_density_long(x, tf, labels)
    if (nrow(long) == 0L) {
      p <- ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5,
                          label = paste("No data for", tf), size = 5) +
        ggplot2::theme_void() +
        ggplot2::ggtitle(tf)
      plot_list[[length(plot_list) + 1L]] <- p
      next
    }
    label_counts <- table(long$label_name)
    new_labels <- paste0(names(label_counts)," (n=", label_counts, ")")
    p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$value,
                                             fill = .data$label_name,
                                             colour = .data$label_name))

    if (fill) {
      p <- p + ggplot2::geom_density(alpha = alpha, adjust = bw_adjust)
    } else {
      p <- p + ggplot2::geom_density(alpha = 0, adjust = bw_adjust, linewidth = 1)
    }

    if (rug) {
      p <- p + ggplot2::geom_rug(alpha = 0.4, colour = "grey40")
    }

    p <- p +
      ggplot2::scale_fill_manual(values = col_map_named, labels = new_labels) +
      ggplot2::scale_colour_manual(values = col_map_named,labels = new_labels) +
      ggplot2::labs(x = "Activity Score", y = "Density", fill = NULL, colour = NULL) +
      ggplot2::xlim(0, 1) +
      ggplot2::ggtitle(tf, subtitle = paste("Dissimilarity score:", signif(diss_score[tf,"MinMax_score"],3)))+
      ggplot2::theme_classic() +
      ggplot2::theme(legend.position.inside = c(1, 1),
                     legend.justification = c(1, 1), # anchor to top-right
                     legend.direction = "vertical",
      plot.title = ggplot2::element_text(face = "bold", size = 12))

    plot_list[[length(plot_list) + 1L]] <- p
  }

  # Paginate & output
  pages <- .paginate_plots(plot_list, plots_per_page, grid_rows, grid_cols)

  if (save) {
    fname <- filename %||% "AUC_distributions.pdf"
    fpath <- file.path(out_dir, fname)
    dir.create(dirname(fpath), recursive = TRUE, showWarnings = FALSE)
    grDevices::pdf(fpath, width = width, height = height, onefile = TRUE)
    for (pg in pages) {
      grid::grid.newpage()
      grid::grid.draw(pg)
    }
    grDevices::dev.off()
    message("Saved AUC distributions to: ", fpath)
  } else {
    for (pg in pages) {
      grid::grid.newpage()
      grid::grid.draw(pg)
    }
  }

  invisible(pages)
}

# ---- Public: plot_auc_cumulative -----------------------------------------

#' Plot empirical cumulative distribution functions of AUC scores
#'
#' @inheritParams plot_auc_distributions
#' @param percentile probability threshold for summary metrics (default 0.5)
#'
#' @return Invisibly, a list of ggplot objects (one per page).
#' @export
plot_auc_cumulative <- function(x,
                                tf_names = NULL,
                                labels = NULL,
                                alpha = 0.8,
                                rug = FALSE,
                                grid = c(4L, 2L),
                                percentile = 0.5,
                                save = FALSE,
                                filename = NULL,
                                out_dir = getwd(),
                                width = 14,
                                height = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
  if (!requireNamespace("gridExtra", quietly = TRUE)) stop("gridExtra required")

  labels   <- .resolve_labels(x, labels)
  tf_names <- .resolve_tf_names(x, tf_names)
  n_tfs    <- length(tf_names)
  lab_names <- x@label_names
  col_map   <- x@colors
  col_map_named <- stats::setNames(unname(col_map), lab_names[names(col_map)])

  if (is.null(grid)) {
    grid_cols <- 2L
    grid_rows <- ceiling(n_tfs / grid_cols)
    plots_per_page <- n_tfs
  } else {
    grid_rows <- as.integer(grid[1])
    grid_cols <- as.integer(grid[2])
    plots_per_page <- grid_rows * grid_cols
  }
  if (is.null(height)) height <- 5 * grid_rows

  # Pre-compute ECDF-AUC metrics for subtitles
  ecdf_df <- tryCatch(
    calculate_ecdf_auc(x, tf_names = tf_names, labels = labels,
                       percentile = percentile),
    error = function(e) NULL
  )

  # col_map <- stats::setNames(
  #   vapply(labels, function(l) .get_label_color(x, l, idx = which(labels == l)),
  #          character(1)),
  #   vapply(labels, function(l) .get_label_name(x, l), character(1))
  # )
  plot_list <- list()
  for (tf in tf_names) {
    long <- .build_density_long(x, tf, labels)

    if (nrow(long) == 0L) {
      p <- ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5,
                          label = paste("No data for", tf), size = 5) +
        ggplot2::theme_void() + ggplot2::ggtitle(tf)
      plot_list[[length(plot_list) + 1L]] <- p
      next
    }
    label_counts <- table(long$label_name)
    new_labels <- paste0(names(label_counts)," (n=", label_counts, ")")

    p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$value,
                                             colour = .data$label_name)) +
      ggplot2::stat_ecdf(geom = "step", linewidth = 1, alpha = alpha) +
      ggplot2::scale_colour_manual(values = col_map_named, labels = new_labels) +
      ggplot2::labs(x = "Activity Score", 
                    y = "Cumulative Probability",
                    colour = NULL) +
      ggplot2::xlim(0, 1) +
      ggplot2::theme_classic() +
      ggplot2::theme(legend.position.inside = c(1, 0),
                     legend.justification = c(1, 0), # anchor to bottom-right
                     legend.direction = "vertical",
                     plot.title = ggplot2::element_text(face = "bold", size = 12))

    if (rug) {
      p <- p + ggplot2::geom_rug(alpha = 0.4, colour = "grey40",
                                  sides = "b")
    }

    # Build title with metrics
    title_parts <- tf
    if (!is.null(ecdf_df) && tf %in% rownames(ecdf_df)) {
      delta_cols <- grep("^delta_", colnames(ecdf_df), value = TRUE)
      sub_parts <- character()
      if ("delta_ecdf_auc" %in% delta_cols) {
        sub_parts <- c(sub_parts, sprintf("Delta*AUC: %.4f",
                                           ecdf_df[tf, "delta_ecdf_auc"]))
      }
      if ("delta_auc50" %in% delta_cols) {
        sub_parts <- c(sub_parts, sprintf("Delta*AUC50: %.4f",
                                           ecdf_df[tf, "delta_auc50"]))
      }
      if (length(sub_parts) > 0L) {
        subtitle_parts <- paste0(paste(sub_parts, collapse = "  |  "))
      }
    }
    p <- p + ggplot2::ggtitle(title_parts, subtitle_parts)

    plot_list[[length(plot_list) + 1L]] <- p
  }

  pages <- .paginate_plots(plot_list, plots_per_page, grid_rows, grid_cols)

  if (save) {
    fname <- filename %||% "AUC_cumulative.pdf"
    fpath <- file.path(out_dir, fname)
    dir.create(dirname(fpath), recursive = TRUE, showWarnings = FALSE)
    grDevices::pdf(fpath, width = width, height = height, onefile = TRUE)
    for (pg in pages) {
      grid::grid.newpage()
      grid::grid.draw(pg)
    }
    grDevices::dev.off()
    message("Saved AUC cumulative to: ", fpath)
  } else {
    for (pg in pages) {
      grid::grid.newpage()
      grid::grid.draw(pg)
    }
  }

  invisible(pages)
}

# ---- Public: plot_auc_summary_statistics ---------------------------------

#' Plot AUC summary statistics (boxplot, violin, mean bar, high-activity count)
#'
#' Creates a 2×2 panel mirroring
#' \code{SimiCVisualization.plot_auc_summary_statistics}.
#'
#' @param x \code{SimiCvizExperiment}
#' @param labels integer vector of labels (default: all)
#' @param high_threshold activity score threshold for "high activity" count
#'   (default 0.5)
#' @param save logical; save to PDF
#' @param filename PDF filename
#' @param out_dir output directory
#' @param width,height page dimensions
#'
#' @return Invisibly, the combined ggplot/grob object.
#' @export
plot_auc_summary_statistics <- function(x,
                                        labels = NULL,
                                        high_threshold = 0.5,
                                        save = FALSE,
                                        filename = NULL,
                                        out_dir = getwd(),
                                        width = 14,
                                        height = 10) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
  if (!requireNamespace("gridExtra", quietly = TRUE)) stop("gridExtra required")

  labels <- .resolve_labels(x, labels)

  # Build long-form data for all TFs and labels
  all_long <- do.call(rbind, lapply(labels, function(lab) {
    sub <- tryCatch(.subset_auc_for_label(x, lab), error = function(e) NULL)
    if (is.null(sub) || ncol(sub) == 0L) return(NULL)
    vals <- unlist(sub, use.names = FALSE)
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0L) return(NULL)
    data.frame(value = vals,
              #  label_name = .get_label_name(x, lab),
               label_name <- simic@label_names[as.character(lab)],
               stringsAsFactors = FALSE)
  }))

  if (is.null(all_long) || nrow(all_long) == 0L) {
    stop("No AUC data available for the requested labels.")
  }

  lab_order <- vapply(labels, function(l) .get_label_name(x, l), character(1))
  all_long$label_name <- factor(all_long$label_name, levels = lab_order)

  col_map <- x@colors[lab_order]

  # 1. Boxplot
  p_box <- ggplot2::ggplot(all_long, ggplot2::aes(x = .data$label_name,
                                                    y = .data$value,
                                                    fill = .data$label_name)) +
    ggplot2::geom_boxplot(alpha = 0.6, outlier.size = 0.5) +
    ggplot2::scale_fill_manual(values = col_map) +
    ggplot2::labs(x = NULL, y = "Activity Score",
                  title = "Activity Score Distribution (Boxplot)") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   plot.title = ggplot2::element_text(face = "bold"))

  # 2. Violin
  p_violin <- ggplot2::ggplot(all_long, ggplot2::aes(x = .data$label_name,
                                                      y = .data$value,
                                                      fill = .data$label_name)) +
    ggplot2::geom_violin(alpha = 0.6, draw_quantiles = c(0.25, 0.5, 0.75)) +
    ggplot2::scale_fill_manual(values = col_map) +
    ggplot2::labs(x = NULL, y = "Activity Score",
                  title = "Activity Score Distribution (Violin)") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   plot.title = ggplot2::element_text(face = "bold"))

  # 3. Mean bar
  summary_df <- do.call(rbind, lapply(labels, function(lab) {
    sub <- .subset_auc_for_label(x, lab)
    vals <- unlist(sub, use.names = FALSE)
    vals <- vals[!is.na(vals)]
    data.frame(label_name = simic@label_names[as.character(lab)],
               mean_val = mean(vals),
               sd_val = stats::sd(vals),
               stringsAsFactors = FALSE)
  }))
  summary_df$label_name <- factor(summary_df$label_name, levels = lab_order)

  p_mean <- ggplot2::ggplot(summary_df, ggplot2::aes(x = .data$label_name,
                                                      y = .data$mean_val,
                                                      fill = .data$label_name)) +
    ggplot2::geom_col(alpha = 0.7, colour = "black", linewidth = 0.5) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$mean_val - .data$sd_val,
                                         ymax = .data$mean_val + .data$sd_val),
                           width = 0.3) +
    ggplot2::scale_fill_manual(values = col_map) +
    ggplot2::labs(x = NULL, y = "Mean Activity Score",
                  title = "Mean Activity Score by Phenotype") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   plot.title = ggplot2::element_text(face = "bold"))

  # 4. High-activity TF count
  ha_df <- do.call(rbind, lapply(labels, function(lab) {
    sub <- .subset_auc_for_label(x, lab)
    mean_per_tf <- colMeans(sub, na.rm = TRUE)
    data.frame(label_name = simic@label_names[as.character(lab)],
               n_high = sum(mean_per_tf > high_threshold, na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
  ha_df$label_name <- factor(ha_df$label_name, levels = lab_order)

  p_high <- ggplot2::ggplot(ha_df, ggplot2::aes(x = .data$label_name,
                                                  y = .data$n_high,
                                                  fill = .data$label_name)) +
    ggplot2::geom_col(alpha = 0.7, colour = "black", linewidth = 0.5) +
    ggplot2::scale_fill_manual(values = col_map) +
    ggplot2::labs(x = NULL, y = "Number of TFs",
                  title = sprintf("TFs with High Activity (Mean AS > %.1f)",
                                  high_threshold)) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   plot.title = ggplot2::element_text(face = "bold"))

  combined <- gridExtra::arrangeGrob(p_box, p_violin, p_mean, p_high,
                                      ncol = 2, nrow = 2,
                                      top = grid::textGrob(
                                        "Activity Scores Summary Statistics",
                                        gp = grid::gpar(fontsize = 16,
                                                         fontface = "bold")))

  if (save) {
    fname <- filename %||% "AUC_summary_statistics.pdf"
    fpath <- file.path(out_dir, fname)
    dir.create(dirname(fpath), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(fpath, combined, width = width, height = height)
    message("Saved AUC summary statistics to: ", fpath)
  } else {
    grid::grid.newpage()
    grid::grid.draw(combined)
  }

  invisible(combined)
}

# ---- Public: plot_auc_heatmap (mean AUC per TF × label) -----------------

#' Plot heatmap of mean AUC per TF across phenotype labels
#'
#' @param x \code{SimiCvizExperiment}
#' @param tf_names character vector of TFs (default: all)
#' @param labels integer vector of labels (default: all)
#' @param top_n integer; show only top N TFs by max cross-label range
#'   (default NULL = all)
#' @param save logical
#' @param filename PDF filename
#' @param out_dir output directory
#' @param width,height page dimensions
#'
#' @return Invisibly, the ggplot object.
#' @export
plot_auc_heatmap <- function(x,
                             tf_names = NULL,
                             labels = NULL,
                             top_n = NULL,
                             save = FALSE,
                             filename = NULL,
                             out_dir = getwd(),
                             width = 10,
                             height = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")

  labels   <- .resolve_labels(x, labels)
  tf_names <- .resolve_tf_names(x, tf_names)

  # Compute mean AUC per TF per label
  rows <- list()
  for (tf in tf_names) {
    row <- list(tf = tf)
    for (lab in labels) {
      sub <- tryCatch(.subset_auc_for_label(x, lab), error = function(e) NULL)
      val <- if (!is.null(sub) && tf %in% colnames(sub)) {
        mean(sub[[tf]], na.rm = TRUE)
      } else {
        NA_real_
      }
      row[[simic@label_names[as.character(lab)]]] <- val
    }
    rows[[length(rows) + 1L]] <- as.data.frame(row, stringsAsFactors = FALSE)
  }
  mat_df <- do.call(rbind, rows)
  rownames(mat_df) <- mat_df$tf

  lab_cols <- setdiff(colnames(mat_df), "tf")
  num_mat <- mat_df[, lab_cols, drop = FALSE]

  # Optional top_n by range
  if (!is.null(top_n) && top_n < nrow(num_mat)) {
    rng <- apply(num_mat, 1, function(r) diff(range(r, na.rm = TRUE)))
    keep <- names(sort(rng, decreasing = TRUE))[seq_len(top_n)]
    num_mat <- num_mat[keep, , drop = FALSE]
  }

  if (is.null(height)) height <- max(6, nrow(num_mat) * 0.35 + 2)

  # Long form
  long <- data.frame(
    tf = rep(rownames(num_mat), times = ncol(num_mat)),
    condition = rep(colnames(num_mat), each = nrow(num_mat)),
    auc = unlist(num_mat, use.names = FALSE),
    stringsAsFactors = FALSE
  )
  long$tf <- factor(long$tf, levels = rev(rownames(num_mat)))
  long$condition <- factor(long$condition, levels = colnames(num_mat))

  p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$condition,
                                            y = .data$tf,
                                            fill = .data$auc)) +
    ggplot2::geom_tile(colour = "grey80") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", .data$auc)),
                       size = 3) +
    ggplot2::scale_fill_viridis_c(option = "C", na.value = "grey90") +
    ggplot2::labs(x = "Phenotype", y = "TF", fill = "Mean AUC",
                  title = "Mean Activity Score per TF") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(face = "bold", size = 14)
    )

  if (save) {
    fname <- filename %||% "AUC_heatmap.pdf"
    fpath <- file.path(out_dir, fname)
    dir.create(dirname(fpath), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(fpath, p, width = width, height = height)
    message("Saved AUC heatmap to: ", fpath)
  } else {
    print(p)
  }

  invisible(p)
}

# ---- Public: export_auc_pdfs (convenience wrapper) -----------------------

#' Export all AUC visualizations to PDF
#'
#' Convenience wrapper that calls \code{\link{plot_auc_distributions}},
#' \code{\link{plot_auc_cumulative}}, \code{\link{plot_auc_summary_statistics}},
#' and \code{\link{plot_auc_heatmap}}, saving each to a \code{plots/auc/}
#' subdirectory.
#'
#' @param x \code{SimiCvizExperiment}
#' @param out_dir root output directory
#' @param prefix filename prefix
#' @param labels integer vector of labels (default: all)
#' @param overwrite logical (currently unused; files are always overwritten)
#' @param ... additional arguments forwarded to individual plot functions
#'
#' @return Invisibly, a list of output file paths.
#' @export
export_auc_pdfs <- function(x,
                            out_dir,
                            prefix = "SimiCviz",
                            labels = NULL,
                            overwrite = FALSE,
                            ...) {
  plot_dir <- file.path(out_dir, "plots", "auc")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  res <- list()

  tryCatch({
    plot_auc_distributions(x, labels = labels, save = TRUE,
                           filename = paste0(prefix, "_auc_distributions.pdf"),
                           out_dir = plot_dir, ...)
    res$distributions <- file.path(plot_dir, paste0(prefix, "_auc_distributions.pdf"))
  }, error = function(e) message("Skipping distributions: ", conditionMessage(e)))

  tryCatch({
    plot_auc_cumulative(x, labels = labels, save = TRUE,
                        filename = paste0(prefix, "_auc_cumulative.pdf"),
                        out_dir = plot_dir, ...)
    res$cumulative <- file.path(plot_dir, paste0(prefix, "_auc_cumulative.pdf"))
  }, error = function(e) message("Skipping cumulative: ", conditionMessage(e)))

  tryCatch({
    plot_auc_summary_statistics(x, labels = labels, save = TRUE,
                                filename = paste0(prefix, "_auc_summary.pdf"),
                                out_dir = plot_dir, ...)
    res$summary <- file.path(plot_dir, paste0(prefix, "_auc_summary.pdf"))
  }, error = function(e) message("Skipping summary: ", conditionMessage(e)))

  tryCatch({
    plot_auc_heatmap(x, labels = labels, save = TRUE,
                     filename = paste0(prefix, "_auc_heatmap.pdf"),
                     out_dir = plot_dir, ...)
    res$heatmap <- file.path(plot_dir, paste0(prefix, "_auc_heatmap.pdf"))
  }, error = function(e) message("Skipping heatmap: ", conditionMessage(e)))

  invisible(res)
}

# ---- Internal: build long-form data for a single TF ---------------------

#' @keywords internal
.build_density_long <- function(x, tf, labels) {
  rows <- list()
  for (lab in labels) {
    sub <- tryCatch(.subset_auc_for_label(x, lab), error = function(e) NULL)
    if (is.null(sub) || !tf %in% colnames(sub)) next
    vals <- stats::na.omit(sub[[tf]])
    if (length(vals) < 1L) next
    rows[[length(rows) + 1L]] <- data.frame(
      value = vals,
      label = lab,
      label_name = unname(x@label_names[as.character(lab)]),
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0L) {
    return(data.frame(value = numeric(), label = integer(),  label_name = character(),
                      stringsAsFactors = FALSE))
  }
  long <- do.call(rbind, rows)
  # Factor preserving label order
  lab_order <- x@label_names
  # lab_order <- vapply(labels, function(l) .get_label_name(x, l), character(1))
  long$label_name <- factor(long$label_name,
                             levels = intersect(lab_order, unique(long$label_name)))
  long
}

# ---- Internal: paginate ggplots into gridExtra pages ---------------------

#' @keywords internal
.paginate_plots <- function(plot_list, per_page, nrow, ncol) {
  n <- length(plot_list)
  pages <- list()
  for (start in seq(1, n, by = per_page)) {
    end <- min(start + per_page - 1L, n)
    batch <- plot_list[start:end]
    # Pad with nullGrob if batch is smaller than grid
    while (length(batch) < nrow * ncol) {
      batch[[length(batch) + 1L]] <- grid::nullGrob()
    }
    pg <- gridExtra::arrangeGrob(grobs = batch, nrow = nrow, ncol = ncol)
    pages[[length(pages) + 1L]] <- pg
  }
  pages
}

# ---- Internal: null-coalescing operator ----------------------------------

#' @keywords internal
`%||%` <- function(a, b) if (!is.null(a)) a else b
