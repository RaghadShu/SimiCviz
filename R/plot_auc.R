#' Plot AUC distribution
#'
#' @param x [SimiCvizExperiment()] or data.frame with columns
#'   \code{tf} and \code{auc}.
#' @param condition optional condition filter if a \code{condition} column
#'   exists.
#'
#' @return ggplot object.
#' @export
plot_auc_distribution <- function(x, condition = NULL) {
  if (is.SimiCvizExperiment(x)) {
    df <- x$auc
  } else {
    df <- x
  }
  if (is.null(df)) stop("No AUC data available.")
  if (!all(c("tf", "auc") %in% colnames(df))) {
    stop("Input must contain columns: tf, auc")
  }

  if (!is.null(condition) && "condition" %in% colnames(df)) {
    df <- df[df$condition == condition, , drop = FALSE]
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$auc)) +
    ggplot2::geom_histogram(bins = 30, fill = "#2166AC", color = "white") +
    ggplot2::labs(x = "AUC", y = "Count") +
    ggplot2::theme_bw()

  p
}

#' Plot AUC heatmap (TF vs condition)
#'
#' Expects columns \code{tf}, \code{condition}, \code{auc}.
#'
#' @param x [SimiCvizExperiment()] or data.frame.
#'
#' @return ggplot object.
#' @export
plot_auc_heatmap <- function(x) {
  if (is.SimiCvizExperiment(x)) {
    df <- x$auc
  } else {
    df <- x
  }
  if (is.null(df)) stop("No AUC data available.")
  required <- c("tf", "condition", "auc")
  if (!all(required %in% colnames(df))) {
    stop("Input must contain columns: tf, condition, auc")
  }

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$condition, y = .data$tf, fill = .data$auc)
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(option = "C") +
    ggplot2::labs(x = "Condition", y = "TF", fill = "AUC") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  p
}

#' Export AUC plots to PDF
#'
#' Exports distribution and, if possible, heatmap plots to a standard layout:
#' \itemize{
#'   \item \code{plots/auc/} for all AUC-related plots.
#' }
#'
#' @param x [SimiCvizExperiment()] or AUC data.frame.
#' @param out_dir root directory.
#' @param prefix filename prefix.
#' @param condition optional condition filter for distribution plot.
#' @param width,height PDF size.
#' @param overwrite logical.
#'
#' @return Invisibly, a list of file paths.
#' @export
export_auc_pdfs <- function(x,
                            out_dir,
                            prefix = "SimiCviz",
                            condition = NULL,
                            width = 7,
                            height = 5,
                            overwrite = FALSE) {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  plot_dir <- file.path(out_dir, "plots", "auc")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  res <- list()

  # Distribution
  suffix <- if (!is.null(condition)) paste0("_", condition) else ""
  f_d <- file.path(plot_dir, paste0(prefix, "_auc_distribution", suffix, ".pdf"))
  if (!overwrite && file.exists(f_d)) {
    stop("File already exists and overwrite = FALSE: ", f_d)
  }
  grDevices::pdf(f_d, width = width, height = height)
  p_d <- plot_auc_distribution(x, condition = condition)
  print(p_d)
  grDevices::dev.off()
  res$distribution <- f_d

  # Heatmap, only if columns available
  df <- if (is.SimiCvizExperiment(x)) x$auc else x
  if (!is.null(df) && all(c("tf", "condition", "auc") %in% colnames(df))) {
    f_h <- file.path(plot_dir, paste0(prefix, "_auc_heatmap.pdf"))
    if (!overwrite && file.exists(f_h)) {
      stop("File already exists and overwrite = FALSE: ", f_h)
    }
    grDevices::pdf(f_h, width = width, height = height)
    p_h <- plot_auc_heatmap(x)
    print(p_h)
    grDevices::dev.off()
    res$heatmap <- f_h
  }

  invisible(res)
}
