.ensure_python <- function() {
  reticulate::py_require(
    python_version = ">=3.8",
    packages = c("numpy", "pandas", "anndata")
  )
}

# Silence R CMD check notes for NSE/data-mask variables used in tidy pipelines.
utils::globalVariables(c(
  "Metric", "Score", "TF", "default_colors", "label", "score",
  "target", "text_color", "tf", "weight"
))

#' @importFrom dplyr %>%
#' @importFrom graphics abline grid hist par plot.new title
#' @importFrom methods as new setClass setGeneric setMethod show
#' @importFrom BiocParallel bplapply MulticoreParam SnowParam SerialParam bpnworkers
#' @importFrom utils head read.csv
NULL