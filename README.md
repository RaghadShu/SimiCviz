# SimiCviz

SimiCviz — Visualization tools for SimiC and SimiCPipeline outputs.

A lightweight R/Bioconductor-oriented package to import, summarize, and visualize
single-cell gene regulatory network (GRN) outputs (weights, TF activity/AUC,
network summaries, and dissimilarity metrics) from SimiC/SimiCPipeline and
other GRN inference tools.

## Status
- Current: Install from GitHub (development).
- Planned: Submit to Bioconductor for formal distribution.

## Installation

Install from GitHub:
```r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("irenemaring/SimiCviz")
```

(When available on Bioconductor, installation instructions will be updated.)

## Quick start

```r
library(SimiCviz)

# Load SimiCvizExperiment directly from SimiCPipeline run output
simic <- load_SimiCPipeline(project_dir = "path/to/SimiCPipeline_run_dir",
                            run_name = "example1",
                            lambda1 = 0.01, lambda2 = 0.001)

# Visualize TF activity distributions for 4 TFs
plot_auc_distributions(
    simic,
    tf_names = simic@tf_ids[1:4],
    fill = TRUE,
    alpha = 0.6,
    bw_adjust = 1 / 8,
    rug = TRUE,
    grid = c(2, 2)
)
```

## Generic CSV + AUC example

```r
library(SimiCviz)

# Example files bundled with the package
weights_path <- system.file("extdata", "example_weights.csv", package = "SimiCviz")
cell_labels_path <- system.file("extdata", "inputFiles/treatment_annotation.csv", package = "SimiCviz")

# Load GRN weights and cell labels
weights_df <- read_weights_csv(weights_path)
cell_labels <- load_cell_labels(cell_labels_path, header = TRUE, sep = ",")

# Load expression matrix from your file (csv/pickle/h5ad/rds)
expression_path <- expression_mat_path <- system.file("extdata",file.path("inputFiles", "example1_expression.pickle"),package = "SimiCviz")

# Compute activity scores (AUC) with AUCProcessor
processor <- AUCProcessor(
    weights = weights_df,
    expression = expression_path,
    cell_labels = cell_labels,
    n_cores = 2,
    backend = "multisession"
)

processor <- compute_auc(
    processor,
    sort_by = "expression",
    percent_of_target = 1.0,
    verbose = TRUE
)

# Activity matrix in wide format (cells x TFs)
auc_wide <- get_auc(processor, format = "wide")
head(auc_wide)
```

## Plot examples

```r
# Build a SimiCvizExperiment from loaded weights + computed AUC
viz_obj <- SimiCvizExperiment(
    weights = weights_df,
    auc = auc_wide,
    cell_labels = cell_labels,
    label_names = c("control", "PD-L1", "DAC", "Combination"),
    colors = c("#e0e0e0", "#a8c8ff", "#ffb6b6", "#c1a9e0")
)

# 1) Weight barplots (top targets per TF)
plot_tf_weights(viz_obj, tf_names = viz_obj@tf_ids[1:4], top_n = 20)

# 2) Dissimilarity heatmap for top TFs
plot_dissimilarity_heatmap(viz_obj, top_n = 10, cmap = "viridis")

# 3) AUC summary heatmap
plot_auc_heatmap(viz_obj, top_n = 20)
```

## Main features
- Load and standardize GRN outputs from SimiCPipeline and generic CSV formats.
- Build and manage `SimiCvizExperiment` containers for weights, AUC/activity,
    cell labels, and metadata.
- Compute activity scores from expression + GRN weights using
    `calculate_activity_scores` or `AUCProcessor`.
- Perform TF-level dissimilarity analysis across labels and optional cell groups.
- Compute ECDF-based comparison metrics with `calculate_ecdf_auc`.
- Import/export tabular results for reproducible analysis workflows.

## Plotting functions
- Weight visualization: `plot_tf_weights`, `plot_target_weights`
- Network view: `plot_tf_network_heatmap`
- Model fit diagnostics (SimiC): `plot_r2_distribution`
- Dissimilarity: `plot_dissimilarity_heatmap`
- Activity distributions: `plot_auc_distributions`, `plot_auc_cumulative`
- Summary views: `plot_auc_heatmap`, `plot_auc_summary_statistics`

## Graphical abstract / Logo
(placeholder — add an image here when available)

![Graphical abstract](docs/figures/graphical_abstract.png)
<!-- Replace the path above with your generated image -->

## Documentation
Vignettes and function documentation are included in the package. After installing, open the main vignette:
```r
browseVignettes("SimiCviz")
```
## Contributing
Contributions, issues and feature requests are welcome. Please open issues or pull requests on the GitHub repository.

## Contact
Irene Marín-Goñi — imarin.4@alumni.unav.es

## License
MIT
