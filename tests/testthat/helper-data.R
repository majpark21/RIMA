# library(SingleCellExperiment)
# library(miloR)
# library(S4Vectors)

make_test_milo <- function() {
  # Use the simulated data from MiloR package for test, see Milo vignette: https://bioconductor.org/packages/release/bioc/vignettes/miloR/inst/doc/milo_demo.html
  sim_trajectory <- readRDS("miloR_sim_trajectory.rds")

  ## Extract SingleCellExperiment object
  traj_sce <- sim_trajectory[['SCE']]
  ## Extract sample metadata to use for testing
  traj_meta <- sim_trajectory[["meta"]]
  ## Add metadata to colData slot
  SummarizedExperiment::colData(traj_sce) <- S4Vectors::DataFrame(traj_meta)
  colnames(traj_sce) <- SummarizedExperiment::colData(traj_sce)$cell_id
  redim <- SingleCellExperiment::reducedDim(traj_sce, "PCA")
  dimnames(redim) <- list(colnames(traj_sce), paste0("PC", c(1:50)))
  SingleCellExperiment::reducedDim(traj_sce, "PCA") <- redim

  SingleCellExperiment::logcounts(traj_sce) <- log(SingleCellExperiment::counts(traj_sce) + 1)

  traj_milo <- miloR::Milo(traj_sce)
  traj_milo <- miloR::buildGraph(traj_milo, k = 10, d = 30)
  traj_milo <- miloR::makeNhoods(traj_milo, prop = 0.1, k = 10, d=30, refined = TRUE)

  return(traj_milo)
}
