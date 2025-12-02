# RIMA: RIgorous Matching of single-cell transcriptomics Atlases

## Overview

RIMA (RIgorous Matching of single-cell transcriptomics Atlases) is an R package for matching single-cell transcriptomics data across two datasets at the level of **cell neighbourhoods**. It is particularly useful for cross-species comparisons, challenging integration scenarios, or comparing different experimental conditions where traditional integration methods may be difficult.

Unlike most integration-based approaches, RIMA identifies one-to-one mappings between cell neighbourhoods (small groups of similar cells) and retains quantitative gene expression information. This allows you to:

- **Compare cell states** across different datasets in a rigorous, statistically sound manner
- **Retain neighbourhood-level expression profiles** for downstream analysis
- **Identify conserved gene expressions** between matched cell populations
- **Analyse associations between neighbourhood's metadata**,  for example to align trajectories

## Key Features

- **Neighbourhood-based matching**: Works at the level of cell neighbourhoods rather than individual cells, improving robustness
- **Statistical significance testing**: Derives p-values for neighbourhood pairs similarities by comparing the to similarities against neighbourhoods where cell identities have been scrambled
- **Flexible gene mapping**: Supports feature mapping between datasets (e.g., for cross-species comparisons using orthologs)
- **Gene conservation scoring**: Identifies genes with conserved expression (CoPE score) across matched neighbourhoods
- **Visualise matching and metadata association**: Includes functions for creating informative heatmaps and embedding visualizations of RIMA's matching

## Installation

You can install RIMA directly from GitHub using the `devtools` package:

```r
# Install devtools if you don't have it
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

# Install RIMA from GitHub
devtools::install_github("majpark21/RIMA")
```

### Requirements

RIMA requires R ≥ 4.0 and the following packages:
- `data.table`
- `miloR` (Bioconductor)
- `SingleCellExperiment` (Bioconductor)
- `Matrix`
- `igraph`
- `ggplot2`
- `viridis`
- `rdist`

The required packages will be automatically installed with the command hereabove.

## Quick Start

Here's a minimal example to get you started:

```r
library(RIMA)
library(miloR)
library(SingleCellExperiment)

# Assuming you have two SingleCellExperiment objects, containing sc data for mouse and rabbit
# sce_mouse and sce_rabbit

# Step 0: Define the neighbourhoods (here with Milo's implementation, but could use others, e.g. metacells)
define_neighbourhoods <- function(sce, prop_seeds, knn=10, reduced.dim="PCA"){
  n_components <- ncol(reducedDim(sce, reduced.dim))  # use all available PCs
  mi <- Milo(sce)
  mi <- miloR::buildGraph(mi, k = knn, d = n_components, reduced.dim = "PCA")
  mi <- miloR::makeNhoods(mi, prop = prop_seeds, k = knn, d=n_components, refined = TRUE)
  return(mi)
}
mi_mouse <- define_neighbourhoods(sce_mouse, prop_seeds = 0.02)
mi_rabbit <- define_neighbourhoods(sce_rabbit, prop_seeds = 0.02)

# Step 1: Preprocess the Milo objects
milos <- preprocess_milos(mi_mouse, mi_rabbit)

# Step 2: Calculate neighbourhood similarities
dt_sims <- calculate_similarities(milos, method = "spearman")

# Step 3: Assess statistical significance of nhood-nhood similarity
dt_sims_sig <- calculate_nhoodnhood_significance(
  milos, dt_sims,
  n_scrambles = 10,
  col_scramble_label = "celltype",
  direction = "b"
)

# Step 4: Match significant nhood-nhood connections
dt_match <- match_nhoods(dt_sims_sig[is_significant == TRUE])

# Step 5: Visualize and analyze results
plot_matches_embed(milos, dt_match, cols_color = c("celltype", "celltype"))
plot_matches_map(milos, dt_match, cols_label = c("celltype", "stage"))

# Example downstream analysis: Find the 3 genes with the most conserved expression across matches
dt_cope <- calculate_cope(milos, dt_match, genes = NULL)
setorder(dt_cope, cope)
plot_paired_expression(milos, dt_match, genes = tail(dt_cope$gene, 3))
```

## Citation

If you use RIMA in your research, please cite:

```
Jacques, M.-A., et al. (2025). RIMA: Rigorous Matching of single-cell transcriptomics Atlases. 
GitHub: https://github.com/majpark21/RIMA
```

## Documentation

- **Vignette**: Run `browseVignettes("RIMA")` or check the [vignettes directory](vignettes/) for detailed usage examples
- **Function help**: Use `?function_name` (e.g., `?calculate_nhoodnhood_significance`) for detailed parameter descriptions
- **Package help**: Run `?RIMA` for an overview

## Example Data

RIMA includes example data from mouse and rabbit gastrulation atlases:
- `sce_mouse_gastrulation`: Mouse gastrulation SingleCellExperiment
- `sce_rabbit_gastrulation`: Rabbit gastrulation SingleCellExperiment

Load them with:
```r
RIMA::sce_mouse_gastrulation
RIMA::sce_rabbit_gastrulation
```

## License

RIMA is licensed under the GPL-3 License. See [LICENSE](LICENSE) for details.

## Contact

For questions or feedback, please contact the package maintainer:
- **Marc-Antoine Jacques** - jacques@ebi.ac.uk

## Related Packages

- [**miloR**](https://bioconductor.org/packages/miloR): Differential abundance testing in neighbourhoods
