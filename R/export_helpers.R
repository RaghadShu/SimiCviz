#' Export all standard visualizations and tables
#'
#' Creates a canonical directory tree under \code{out_dir}:
#' \itemize{
#'   \item \code{data/weights}, \code{data/auc}
#'   \item \code{plots/network}, \code{plots/auc}
#' }
#'
#' @param x [SimiCvizExperiment()] object.
#' @param out_dir root directory for exports.
#' @param prefix filename prefix.
#' @param conditions optional vector of conditions to loop over. If \code{NULL},
#'   no condition stratification is applied.
#' @param overwrite logical; passed to underlying exporters.
#'
#' @return Invisibly, a list with data and plot file paths.
#' @export
export_SimiCviz_all <- function(x,
                                out_dir,
                                prefix = "SimiCviz",
                                conditions = NULL,
                                overwrite = FALSE) {
  if (!is.SimiCvizExperiment(x)) {
    stop("x must be a SimiCvizExperiment")
  }
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # Tables
  data_files <- export_SimiCviz_csv(
    x,
    out_dir = out_dir,
    prefix = prefix,
    overwrite = overwrite
  )

  # Network plots
  net_files <- list()
  if (is.null(conditions) || !"condition" %in% colnames(x$weights)) {
    net_files$all <- export_network_pdf(
      x,
      out_dir = out_dir,
      prefix = prefix,
      condition = NULL,
      overwrite = overwrite
    )
  } else {
    for (cond in conditions) {
      net_files[[cond]] <- export_network_pdf(
        x,
        out_dir = out_dir,
        prefix = prefix,
        condition = cond,
        overwrite = overwrite
      )
    }
  }

  # AUC plots
  auc_files <- list()
  if (!is.null(x$auc)) {
    if (is.null(conditions) || !"condition" %in% colnames(x$auc)) {
      auc_files$all <- export_auc_pdfs(
        x,
        out_dir = out_dir,
        prefix = prefix,
        condition = NULL,
        overwrite = overwrite
      )
    } else {
      for (cond in conditions) {
        auc_files[[cond]] <- export_auc_pdfs(
          x,
          out_dir = out_dir,
          prefix = prefix,
          condition = cond,
          overwrite = overwrite
        )
      }
    }
  }

  invisible(list(
    data = data_files,
    network = net_files,
    auc = auc_files
  ))
}
