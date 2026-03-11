#' Plot dissimilarity heatmap
#'
#' Displays a heatmap of TF regulatory dissimilarity scores, optionally
#' broken down by cell group. Uses \pkg{ggplot2} tile geometry.
#'
#' @param x A \code{SimiCvizExperiment} object.
#' @param labels Integer vector of label keys or character vector of label
#'   display names to compare (default: all). Mixing types is allowed.
#' @param top_n Integer; number of top TFs to display (default: all).
#' @param cell_groups Optional named list of cell-ID vectors (see
#'   \code{\link{calculate_dissimilarity}}).
#' @param sort_by Column name to sort TFs by.
#'   Default: \code{"MinMax_score"} (no groups) or \code{"mean_score"} (groups).
#' @param dissim_df Optional pre-computed dissimilarity data.frame. If
#'   \code{NULL}, \code{\link{calculate_dissimilarity}} is called internally.
#' @param cmap Colour palette specification. Can be:
#'   \itemize{
#'     \item A viridis palette name: \code{"viridis"}, \code{"magma"},
#'       \code{"plasma"}, \code{"inferno"}, or \code{"cividis"}
#'       (requires \pkg{viridisLite}).
#'     \item A single colour string (gradient from white to that colour).
#'     \item A character vector of 2+ colours for a custom gradient.
#'     \item \code{NULL} (default): built-in viridis-like gradient.
#'   }
#' @param show_values Logical; annotate cells with numeric values (default \code{TRUE}).
#' @param save Logical; save the plot to a PDF file (default \code{TRUE}).
#' @param out_dir Output directory for the PDF (default: working directory).
#' @param filename Custom filename (default: auto-generated).
#' @param width,height PDF dimensions in inches.
#' @param ... Additional arguments passed to \code{calculate_dissimilarity}.
#'
#' @return Invisibly, a list with \code{plot} (the ggplot object) and
#'   \code{data} (the dissimilarity data.frame).
#'
#' @examples
#' \dontrun{
#'   plot_dissimilarity_heatmap(simic, top_n = 20)
#'   plot_dissimilarity_heatmap(simic, top_n = 15, cmap = "magma", save = FALSE)
#'   plot_dissimilarity_heatmap(simic, cmap = c("white", "red", "darkred"))
#' }
#' @import ggplot2
#' @rdname plot_dissimilarity_heatmap
#' @export
plot_dissimilarity_heatmap <- function(x,
                                       labels = NULL,
                                       top_n = NULL,
                                       cell_groups = NULL,
                                       sort_by = NULL,
                                       dissim_df = NULL,
                                       cmap = NULL,
                                       show_values = TRUE,
                                       save = TRUE,
                                       out_dir = getwd(),
                                       filename = NULL,
                                       width = NULL,
                                       height = NULL,
                                       ...) {

  if (!is.SimiCvizExperiment(x)) {
    stop("plot_dissimilarity_heatmap: 'x' must be a SimiCvizExperiment.")
  }

  # Compute or use pre-computed scores
  if (is.null(dissim_df)) {
    dissim_df <- calculate_dissimilarity(x, labels = labels,
                                         cell_groups = cell_groups,
                                         verbose = FALSE, ...)
  }
  if (is.null(dissim_df) || nrow(dissim_df) == 0L) {
    stop("No dissimilarity scores to plot.")
  }

  # Sort
  if (is.null(sort_by)) {
    sort_by <- if ("mean_score" %in% colnames(dissim_df)) "mean_score" else colnames(dissim_df)[1]
  }
  if (!(sort_by %in% colnames(dissim_df))) {
    warning(sprintf("sort_by='%s' not found; using first column.", sort_by))
    sort_by <- colnames(dissim_df)[1]
  }
  dissim_df <- dissim_df[order(dissim_df[[sort_by]], decreasing = FALSE), , drop = FALSE]

  # Top N (after sorting descending, take top_n then re-sort for plot)
  if (!is.null(top_n) && top_n < nrow(dissim_df)) {
    # Take the top_n highest (they are at the bottom after ascending sort)
    dissim_df <- utils::tail(dissim_df, top_n)
  }

  n_rows <- nrow(dissim_df)
  n_cols <- ncol(dissim_df)

  # Prettify column names for display
  col_display <- colnames(dissim_df)
  col_display[col_display == "mean_score"]   <- "Mean"
  col_display[col_display == "MinMax_score"] <- "MinMax Score"
  col_map <- stats::setNames(col_display, colnames(dissim_df))

  # Convert to long format for ggplot
  dissim_df$TF <- factor(rownames(dissim_df), levels = rownames(dissim_df))
  long_df <- reshape2::melt(dissim_df, id.vars = "TF",
                            variable.name = "Metric", value.name = "Score")
  # Apply display names
  long_df$Metric <- factor(col_map[as.character(long_df$Metric)],
                           levels = col_display)

  # Build colour scale
  fill_scale <- .build_ggplot_fill_scale(cmap)

  # Auto dimensions
  if (is.null(width))  width  <- max(5, 1.5 + n_cols * 1.5)
  if (is.null(height)) height <- max(5, n_rows * 0.35 + 2)

  # Title
  title_txt <- if (!is.null(cell_groups)) {
    "Regulatory Dissimilarity Scores by Cell Group"
  } else {
    "Regulatory Dissimilarity Scores"
  }

  # Determine text colour threshold
  vmin <- min(long_df$Score, na.rm = TRUE)
  vmax <- max(long_df$Score, na.rm = TRUE)
  thresh <- (vmin + vmax) / 2

  p <- ggplot2::ggplot(long_df, ggplot2::aes(x = Metric, y = TF, fill = Score)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    fill_scale +
    ggplot2::labs(title = title_txt,
                  x = NULL, y = "Transcription Factors",
                  fill = "Dissimilarity\nScore") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y  = ggplot2::element_text(size = 8),
      plot.title    = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5),
      panel.grid    = ggplot2::element_blank(),
      legend.position = "right"
    )

  if (show_values) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.4f", Score),
                   colour = ggplot2::after_stat(ifelse(Score > thresh, "high", "low"))),
      size = 2.8, show.legend = FALSE
    )
    # Use manual colour so text is readable on both light and dark tiles
    p <- p + ggplot2::scale_colour_manual(
      values = c("high" = "white", "low" = "black"),
      guide = "none"
    )
    # Replace the geom_text above with a simpler approach that avoids after_stat
    # Remove the last two layers and redo
    p <- ggplot2::ggplot(long_df, ggplot2::aes(x = Metric, y = TF, fill = Score)) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
      fill_scale +
      ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%.4f", Score)),
        colour = ifelse(long_df$Score > thresh, "white", "black"),
        size = 2.8
      ) +
      ggplot2::labs(title = title_txt,
                    x = NULL, y = "Transcription Factors",
                    fill = "Dissimilarity\nScore") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y  = ggplot2::element_text(size = 8),
        plot.title    = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5),
        panel.grid    = ggplot2::element_blank(),
        legend.position = "right"
      )
  }

  if (save) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    fname <- if (is.null(filename)) "dissimilarity_heatmap.pdf" else filename
    fpath <- file.path(out_dir, fname)
    ggplot2::ggsave(fpath, plot = p, width = width, height = height)
    message(sprintf("Saved dissimilarity heatmap to %s", fpath))
  } else {
    print(p)
  }

  # Remove the TF column we added for melting
  dissim_df$TF <- NULL

  invisible(list(plot = p, data = dissim_df))
}


#' Build a ggplot2 fill scale from a cmap specification
#'
#' @param cmap Palette spec (NULL, viridis name, single colour, or vector).
#' @return A ggplot2 scale object.
#' @keywords internal
#' @noRd
.build_ggplot_fill_scale <- function(cmap) {
  viridis_names <- c("viridis", "magma", "plasma", "inferno", "cividis",
                     "rocket", "mako", "turbo")

  # Default: viridis-like

  if (is.null(cmap)) {
    return(ggplot2::scale_fill_gradientn(
      colours = c("#440154", "#31688e", "#35b779", "#fde725")
    ))
  }

  # Named viridis palette
  if (is.character(cmap) && length(cmap) == 1L &&
      tolower(cmap) %in% viridis_names) {
    palette_name <- tolower(cmap)
    if (requireNamespace("viridisLite", quietly = TRUE)) {
      cols <- viridisLite::viridis(256, option = palette_name)
      return(ggplot2::scale_fill_gradientn(colours = cols))
    }
    # Fallback approximations
    fallback <- list(
      viridis = c("#440154", "#31688e", "#35b779", "#fde725"),
      magma   = c("#000004", "#51127c", "#b73779", "#fc8961", "#fcfdbf"),
      plasma  = c("#0d0887", "#7e03a8", "#cc4778", "#f89540", "#f0f921"),
      inferno = c("#000004", "#420a68", "#932667", "#dd513a", "#fca50a", "#fcffa4"),
      cividis = c("#00224e", "#414d6b", "#7b7b78", "#b8a951", "#fdea45"),
      rocket  = c("#03051a", "#4c1d4e", "#a11a5b", "#e04f39", "#faebdd"),
      mako    = c("#0b0405", "#2a1858", "#245f8a", "#30b09a", "#def5e5"),
      turbo   = c("#30123b", "#4662d7", "#35abf8", "#1ae4b6", "#72fe5e",
                  "#c8ef34", "#faba39", "#f66b19", "#ca240e", "#7a0403")
    )
    anchors <- fallback[[palette_name]] %||% fallback[["viridis"]]
    message(sprintf("viridisLite not installed; using approximate '%s' palette.", palette_name))
    return(ggplot2::scale_fill_gradientn(colours = anchors))
  }

  # Single colour: white → colour
  if (is.character(cmap) && length(cmap) == 1L) {
    return(ggplot2::scale_fill_gradient(low = "white", high = cmap))
  }

  # Custom vector of colours
  if (is.character(cmap) && length(cmap) >= 2L) {
    return(ggplot2::scale_fill_gradientn(colours = cmap))
  }

  warning("Unrecognised `cmap`; using default viridis-like palette.")
  ggplot2::scale_fill_gradientn(
    colours = c("#440154", "#31688e", "#35b779", "#fde725")
  )
}
