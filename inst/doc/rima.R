## ----echo=FALSE---------------------------------------------------------------
# Go around issue when a data.table is not printed in knitted doc
knit_print.data.table = function(x, ...) {
  res = paste(c("", "", knitr::kable(as.data.frame(x))), collapse = "\n")
  knitr::asis_output(res)
}

knit_print.data.frame = function(x, ...) {
  res = paste(c("", "", knitr::kable(x)), collapse = "\n")
  knitr::asis_output(res)
}

registerS3method(
  "knit_print", "data.frame", knit_print.data.frame,
  envir = asNamespace("knitr")
)

registerS3method(
  "knit_print", "data.table", knit_print.data.table,
  envir = asNamespace("knitr")
)

## -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(RIMA)
  library(SingleCellExperiment)
  library(miloR)
  library(data.table)
  library(ggplot2)
  library(scater)
})
set.seed(7)

# SingleCellExperiment objects containing the gastrulation data
sce_mouse <- RIMA::sce_mouse_gastrulation
sce_rabbit <- RIMA::sce_rabbit_gastrulation

# Color scheme for cell types
l_ct_colors <- lapply(list(mouse = sce_mouse, rabbit = sce_rabbit), function(x) {
  as.data.table(unique(colData(x)[, c("celltype", "color")]))
})
dt_ct_palette <- unique(rbindlist(l_ct_colors))

## ----message=FALSE------------------------------------------------------------
# PCA is a batch-corrected embedding, cells that are close to each other in this embedding will form neighbourhoods
define_neighbourhoods <- function(sce, prop_seeds, knn=10, reduced.dim="PCA"){
  n_components <- ncol(reducedDim(sce, reduced.dim))  # use all available PCs
  mi <- Milo(sce)
  mi <- miloR::buildGraph(mi, k = knn, d = n_components, reduced.dim = "PCA")
  mi <- miloR::makeNhoods(mi, prop = prop_seeds, k = knn, d=n_components, refined = TRUE)
  return(mi)
}

# prop_seeds controls the proportion of cells selected as neighbourhoods' seeds. Decrease to have less neighbourhoods.
mi_mouse <- define_neighbourhoods(sce_mouse, prop_seeds = 0.02)
mi_rabbit <- define_neighbourhoods(sce_rabbit, prop_seeds = 0.02)

## -----------------------------------------------------------------------------
milos <- preprocess_milos(mi_mouse, mi_rabbit)
milos

## -----------------------------------------------------------------------------
dt_sims <- RIMA::calculate_similarities(milos, method = "spearman")
head(dt_sims)

## ----message=FALSE,results='hide'---------------------------------------------
dt_sims_withSignif <- RIMA::calculate_nhoodnhood_significance(
  milos,
  dt_sims,
  n_scrambles = 10,
  col_scramble_label = "celltype",
  alpha = 0.05,
  direction = "b"
)

## -----------------------------------------------------------------------------
head(dt_sims_withSignif)

## ----fig.height=8, fig.width=10, message=FALSE--------------------------------
RIMA::plot_trimming_map(milos, dt_sims_withSignif, cols_label = c("celltype", "celltype")) + labs(x="Celltype rabbit", y="Celltype mouse")
RIMA::plot_trimming_map(milos, dt_sims_withSignif, cols_label = c("celltype", "stage")) + labs(x="Stage rabbit", y="Celltype mouse")

## -----------------------------------------------------------------------------
dt_match <- RIMA::match_nhoods(dt_sims_withSignif[is_significant == TRUE])
head(dt_match)

## -----------------------------------------------------------------------------
# UMAP 2D embedding for visualisation
milos <- lapply(milos, function(x) {
  scater::runUMAP(x, ncomponents = 2, dimred = "PCA")
})

## ----fig.width=6, fig.height=5------------------------------------------------
p <- RIMA::plot_matches_embed(milos, dt_match, dt_palette = dt_ct_palette, cols_color = c("celltype", "celltype")) + guides(fill="none")
p

## ----message=F, fig.height=8, fig.width=10------------------------------------
RIMA::plot_matches_map(milos, dt_match, cols_label = c("celltype", "celltype"), transform="log10")  + labs(x="Celltype rabbit", y="Celltype mouse")
RIMA::plot_matches_map(milos, dt_match, cols_label = c("celltype", "stage"), cluster_cols = FALSE)  + labs(x="Stage rabbit", y="Celltype mouse")

## ----warning=F----------------------------------------------------------------
# Set genes to NULL to calculate CoPE of all the genes
dt_copes <- RIMA::calculate_cope(milos, dt_match, genes = NULL)

# Most and least conserved gene expression across matched neighbourhoods
setorder(dt_copes, cope)
top_genes <- tail(dt_copes[!is.na(cope), gene], 3)
worst_genes <- head(dt_copes[!is.na(cope), gene], 3)

dt_copes[gene %in% c(top_genes, worst_genes)]

## ----warning=F----------------------------------------------------------------
RIMA::plot_paired_expression(milos, dt_match, c(top_genes, worst_genes)) + scale_x_log10() + scale_y_log10() + labs(x="Expression mouse [log10]", y="Expression rabbit [log10]")

