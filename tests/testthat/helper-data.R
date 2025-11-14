make_test_milo <- function() {
  # Use the simulated data from MiloR package for test, see Milo vignette: https://bioconductor.org/packages/release/bioc/vignettes/miloR/inst/doc/milo_demo.html
  # Has been turned into a SingleCellExperiment object, slimmed down to 100 genes, and removed PCA. Reduce size package.
  set.seed(7)

  traj_sce <- readRDS("miloR_sim_trajectory_slim.rds")

  traj_milo <- miloR::Milo(traj_sce)
  traj_milo <- miloR::buildGraph(traj_milo, k = 10, d = 30)
  traj_milo <- miloR::makeNhoods(traj_milo, prop = 0.1, k = 10, d=30, refined = TRUE)

  return(traj_milo)
}
