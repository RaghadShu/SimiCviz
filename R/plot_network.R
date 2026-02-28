#' Plot regulatory network from weights
#'
#' Creates a simple TF-target regulatory network using \pkg{igraph}.
#'
#' @param x [SimiCvizExperiment()] or data.frame with columns
#'   \code{tf}, \code{target}, \code{weight}.
#' @param condition optional condition filter if a \code{condition} column
#'   exists.
#' @param top_n optional integer; keep only top \code{top_n} edges by
#'   absolute weight.
#' @param layout_fn igraph layout function (default: [igraph::layout_with_fr]).
#'
#' @return igraph plot (invisibly returns the igraph object).
#' @export
plot_network <- function(x,
                         condition = NULL,
                         top_n = 500,
                         layout_fn = igraph::layout_with_fr) {
  if (is.SimiCvizExperiment(x)) {
    df <- x@weights
  } else {
    df <- x
  }
  if (!all(c("tf", "target", "weight") %in% colnames(df))) {
    stop("Input must contain columns: tf, target, weight")
  }

  if (!is.null(condition) && "condition" %in% colnames(df)) {
    df <- df[df$condition == condition, , drop = FALSE]
  }

  if (!is.null(top_n) && nrow(df) > top_n) {
    ord <- order(abs(df$weight), decreasing = TRUE)
    df <- df[ord[seq_len(top_n)], , drop = FALSE]
  }

  g <- igraph::graph_from_data_frame(
    d = data.frame(
      from = df$tf,
      to = df$target,
      weight = df$weight,
      stringsAsFactors = FALSE
    ),
    directed = TRUE
  )

  w <- igraph::E(g)$weight
  e_col <- ifelse(w >= 0, "#2166AC", "#B2182B")
  e_width <- scales::rescale(abs(w), to = c(0.5, 4))

  plot(
    g,
    layout = layout_fn(g),
    vertex.label = igraph::V(g)$name,
    vertex.size = 5,
    vertex.label.cex = 0.6,
    edge.arrow.size = 0.4,
    edge.color = e_col,
    edge.width = e_width
  )

  invisible(g)
}

#' Export network plot to PDF
#'
#' @param x [SimiCvizExperiment()] or weights data.frame.
#' @param out_dir root output directory.
#' @param prefix filename prefix.
#' @param condition optional condition filter.
#' @param top_n integer; passed to [plot_network()].
#' @param width,height PDF dimensions (in inches).
#' @param overwrite logical; overwrite if file exists.
#'
#' @return Invisibly, the PDF filepath.
#' @export
export_network_pdf <- function(x,
                               out_dir,
                               prefix = "SimiCviz",
                               condition = NULL,
                               top_n = 500,
                               width = 7,
                               height = 7,
                               overwrite = FALSE) {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  plot_dir <- file.path(out_dir, "plots", "network")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  suffix <- if (!is.null(condition)) paste0("_", condition) else ""
  f <- file.path(plot_dir, paste0(prefix, "_network", suffix, ".pdf"))
  if (!overwrite && file.exists(f)) {
    stop("File already exists and overwrite = FALSE: ", f)
  }

  grDevices::pdf(f, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  plot_network(x, condition = condition, top_n = top_n)

  invisible(f)
}
