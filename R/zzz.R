# .onLoad <- function(libname, pkgname) {
#   reticulate::py_require(
#     python_version = ">=3.8",
#     packages = c("numpy", "pandas")  # whatever you need
#   )
# }
.ensure_python <- function() {
  reticulate::py_require(
    python_version = ">=3.8",
    packages = c("numpy", "pandas")
  )
}

