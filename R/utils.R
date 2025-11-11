library(dplyr)
library(S4Vectors)
library(data.table)
library(igraph)
library(scran)

#' Left-join colData or rowData with a DataFrame
#' 
#' This function is here to ensure that the order of cells (colData) or features
#' (rowData) is well-preserved after merging with dataframe or datatable. Also
#' insures that colnames/rownames are not lost after merging.
#' 
#' Arguments inherited from dplyr::left_join.
#' 
#' sce: sce object
#' df: dataframe or datatable to join
#' data: "col" or "row", whether to perform the join on colData or rowData
#' by: column name to merge on, should be present in both col/rowData(sce) and df
#' 
#' Example
#' colData(sce) <- left_join(colData(sce), dt_join, by="celltype")
#' 
left_join_sce <-
  function(sce,
           df,
           data = c("col", "row"),
           by = NULL,
           copy = FALSE,
           suffix = c(".x", ".y"),
           ...,
           keep = FALSE) {
    
    if(!data %in% c("col", "row")){
      stop("'data' must be one of c('col', 'row').")
    }
    
    if(data=="col"){
      v_names <- colnames(sce)
      colData(sce) <-
        S4Vectors::DataFrame(
          dplyr::left_join(
            as.data.frame(colData(sce)),
            df,
            by = by,
            copy = copy,
            suffix = suffix,
            keep = keep,
            ...
          )
        )
      colnames(sce) <- v_names
    } else if(data=="row"){
      v_names <- rownames(sce)
      rowData(sce) <-
        S4Vectors::DataFrame(
          dplyr::left_join(
            as.data.frame(rowData(sce)),
            df,
            by = by,
            copy = copy,
            suffix = suffix,
            keep = keep,
            ...
          )
        )
      rownames(sce) <- v_names
    }
    
    return(sce)
  }



# MNN matching ----------------------------------------------------------------------------

# Calculate pairwise nhood correlation and format in long data.table
getSims <- function(milo1, milo2, subset_row1=NULL, subset_row2=NULL, sim="spearman"){
  # Extract nhood expression matrix from Milo or use matrix directly
  if(is.matrix(milo1)){
    nhoods1 <- milo1
  } else if(class(milo1)=="Milo"){
    nhoods1 <- nhoodExpression(milo1)    
  } else {
    stop("'milo1' must be of class Milo or a matrix.")
  }
  if(is.matrix(milo2)){
    nhoods2 <- milo2
  } else if(class(milo2)=="Milo"){
    nhoods2 <- nhoodExpression(milo2)    
  } else {
    stop("'milo2' must be of class Milo or a matrix.")
  }
  
  if(!is.null(subset_row1)){
    nhoods1 <- nhoods1[subset_row1,]
  }
  if(!is.null(subset_row2)){
    nhoods2 <- nhoods2[subset_row2,]
  }
  out <- cor(
    nhoods1,
    nhoods2,
    method = sim
  )
  
  out <- as.data.table(reshape2::melt(as.matrix(out)))
  setnames(out, c("nhoods1", "nhoods2", "sim"))
  return(out)
}


# dt_sim is a long-format data table with 3 columns: nhood center 1, nhood center 2, similarity values; Usually obtained by converting a pairwise correlation matrix. e.g.: as.data.table(reshape2::melt(as.matrix(cor(x))))
# col_ref is the group where every cell is getting matched to the other column
# col_sim is the column holding similarities
# eps_diffMax is to return all hits that are comprised between [MaxSim-eps, MaxSim]
getMappings <- function(dt_sim, col_ref="Var1", col_sim="sim", eps_diffMax=0){
  nhood_sim <- copy(dt_sim)
  setnames(nhood_sim, col_sim, "sim")
  nhood_sim[, thresh := max(.SD[, sim]) - eps_diffMax, by=col_ref]
  nhood_sim <- nhood_sim[sim >= thresh]
  nhood_sim[, thresh := NULL]
  
  return(nhood_sim)
}


# dt_mapping is a 3 column data.table with IDs of nhoods in species 1, IDs of nhoods in species 2, similarity between nhoods
# When eps_diffMax=0 each nhood will appear exactly once. If > 0 some nhoods are "candidate matches" for several nhoods 
getMNNs <- function(dt_mapping, col_sim="sim", eps_diffMax = 0){
  cols_ref <- setdiff(colnames(dt_mapping), col_sim)
  map1 <- getMappings(dt_mapping, col_ref = cols_ref[1], eps_diffMax=eps_diffMax)
  map2 <- getMappings(dt_mapping, col_ref = cols_ref[2], eps_diffMax=eps_diffMax)
  
  mnn <- rbindlist(list(map1[, ..cols_ref], map2[, ..cols_ref]))
  mnn <- mnn[duplicated(mnn)]
  return(mnn)
}

# dt_mnn is a data table with 2 columns: first column comprises nhood names in the source species, target first column comprises nhood names in the target species
# dt_sim is a long-format data table with 3 columns: nhood center 1, nhood center 2, similarity values; Usually obtained by converting a pairwise correlation matrix. e.g.: as.data.table(reshape2::melt(as.matrix(cor(x))))
# Returns an igraph matching list if return_matchObject=T
getBipartiteMatch <- function(dt_mnn, dt_sim, return_matchObject=F){
  # "weight" is automatically recognized as edge attribute
  dt_mnn <- merge(dt_mnn, dt_sim, by=colnames(dt_mnn))
  setnames(dt_mnn, "sim", "weight")
  original_cols <- setdiff(colnames(dt_mnn), "weight")
  # Ensure names of nodes are unique for each species
  for(ii in seq_along(original_cols)){
    dt_mnn[[original_cols[ii]]] <- paste(paste0(".GROUP", ii), dt_mnn[[original_cols[ii]]], sep="__")
  }
  
  # vertex type, one species is coded as up, the other as down
  dt_types <- list(up=unique(dt_mnn[[1]]), down=unique(dt_mnn[[2]]))
  dt_types <- rbindlist(lapply(dt_types, as.data.table, keep.rownames=T), idcol = "type")
  setnames(dt_types, "V1", "name")
  setcolorder(dt_types, "name")
  dt_types[, type := type=="up"]  # convert to Bool for matching
  
  g <- graph_from_data_frame(dt_mnn, directed = F, vertices = dt_types)
  # NULL weights and types will read directly from graph's attributes
  match <- max_bipartite_match(g, weights = NULL, types = NULL)
  if(return_matchObject){
    return(match)
  }
  
  out <- as.data.table(match$matching, keep.rownames=T)
  # Remove cells for which there is no match left
  out <- out[(!is.na(V1)) & (!is.na(V2))]
  # Each match is listed twice, because it is the same in both direction. Keep only one
  out <- out[str_detect(V1, ".GROUP1")]
  out[, V1 := str_replace(V1, ".GROUP1__", "")]
  out[, V2 := str_replace(V2, ".GROUP2__", "")]
  setnames(out, original_cols)
  
  return(out)
}

# l_mnn_eps <- lapply(l_sims, function(x) getMNNs(x, eps_diffMax = 0.01))
# out <- getBipartiteMatch(l_mnn_eps[["mfascicularis"]], l_sims[["mfascicularis"]])


# Plot MNN matching -------------------------------------------------------------

# Extract dimred coords of representative cells from each nhood
extractCentroidNhoods <- function(milo, dimred="UMAP", use_indexCell=T){
  mat <- reducedDim(milo, dimred)
  
  if(use_indexCell){
    l_nhood_centers <- lapply(milo@nhoodIndex, function(x) colnames(milo)[x])
    names(l_nhood_centers) <- unlist(milo@nhoodIndex)
    l_nhood_centers <- lapply(l_nhood_centers, function(x) mat[x,])
  } else {
    l_nhood_cells <- nhoodCells(milo)
    l_nhood_coords <- lapply(l_nhood_cells, function(x) mat[x,])
    l_nhood_avg <- lapply(l_nhood_coords, function(x) matrix(colMeans(x), ncol = 2))
    
    l_nhood_centers <- lapply(names(l_nhood_avg), function(x) rdist::cdist(l_nhood_coords[[x]], l_nhood_avg[[x]], metric = "euclidean"))
    names(l_nhood_centers) <- names(l_nhood_avg)
    l_nhood_centers <- lapply(l_nhood_centers, which.min)
    l_nhood_centers <- lapply(names(l_nhood_centers), function(x) l_nhood_coords[[x]][l_nhood_centers[[x]],])
    names(l_nhood_centers) <- names(l_nhood_avg)
  }
  
  l_nhood_centers <- lapply(l_nhood_centers, function(x) data.table(matrix(x, ncol=2)))
  out <- rbindlist(l_nhood_centers, idcol = "nhood")
  out[, V1 := as.numeric(V1)]
  out[, V2 := as.numeric(V2)]
  return(out)
}

# Returns number of cells in each nhood
nhoodSize <- function(milo){
  out <- colSums(nhoods(milo))
  out <- as.data.table(out, keep.rownames = T)
  setnames(out, c("nhood", "Ncell"))
  return(out)
}

# Annotate nhoods with a column in colData() of the Milo
nhoodAnnot <- function(milo, cols_annot=c("celltype")){
  l_nhood_celltype <- list()
  nhoods_sce <- nhoods(milo)
  nhood_annot <- data.table(Nhood = 1:ncol(nhoods_sce), Nhood_center = colnames(nhoods_sce))
  l_annot <- lapply(cols_annot, function(x) miloR::annotateNhoods(milo, nhood_annot, coldata_col = x))
  merge_by <- function(x,y){merge(x,y,by=c("Nhood", "Nhood_center"))}
  dt_annot <- Reduce(merge_by, l_annot)
  # nhood_annot <- miloR::annotateNhoods(milo, nhood_annot, coldata_col = cols_annot)
  return(dt_annot)
}

# Get a list with cells composing each nhood in the milo object
nhoodCells <- function(milo){
  l_nhood_cells <- apply(nhoods(milo), 2, function(x){names(x)[x==1]})
  return(l_nhood_cells)
}

plotNhoods <- function(milo, color_by=NULL, alpha_points=0.5, fill_breaks_colors=NULL, ...){
  dt_plot <- extractCentroidNhoods(milo, ...)
  dt_plot <- merge(dt_plot, nhoodSize(milo))
  if(!is.null(color_by)){
    dt_annot <- nhoodAnnot(milo, color_by)
    cols_merge <- c("Nhood_center", color_by)
    dt_plot <- merge(dt_plot, dt_annot[, ..cols_merge], by.x="nhood", by.y="Nhood_center")
  }
  p <- ggplot(dt_plot, aes(x=V1, y=V2)) +
    geom_point(aes(fill=.data[[color_by]], size=Ncell), color="black", pch=21, alpha=alpha_points) +
    theme_minimal() +
    guides(fill="none", size="none") +
    labs(x="", y="") +
    theme(panel.grid = element_blank(), axis.text.x = element_blank(), axis.text.y=element_blank())
  
  if(!is.null(fill_breaks_colors)){
    p <- p + scale_fill_manual(breaks = fill_breaks_colors[[1]], values=fill_breaks_colors[[2]])
  }
  
  return(p)
}


# Create a set of coordinates where both manifolds are in the same coordinates range and add an offset
# dt_matches is a data table with 2 columns indicating, one for the IDs of each set. 1st column must refer to milo1, 2nd column must refer to milo2
plotManifoldMapping <- function(milo1, milo2, dt_matches, color_by=NULL, alpha_points=0.5, alpha_lines=0.8, hjust=1.2, vjust=0, fill_breaks_colors=NULL, draw_lines=T, color_unmapped=T, size_byNcell=T, ...) {
  prepare_dt_plot <- function(milo, color_by=NULL, ...) {
    dt_plot <- extractCentroidNhoods(milo, ...)
    dt_plot <- merge(dt_plot, nhoodSize(milo))
    if (!is.null(color_by)) {
      dt_annot <- nhoodAnnot(milo, color_by)
      cols_merge <- c("Nhood_center", color_by)
      dt_plot <- merge(dt_plot, dt_annot[, ..cols_merge], by.x = "nhood", by.y = "Nhood_center")
    }
    return(dt_plot)
  }
  
  dt_plot1 <- prepare_dt_plot(milo1, color_by = color_by, ...)
  dt_plot2 <- prepare_dt_plot(milo2, color_by = color_by, ...)
  
  # Rescale coords and offset one manifold
  dt_plot1[, c("V1", "V2") := list(scales::rescale(V1, to=c(0,1)), scales::rescale(V2, to=c(0,1)))]
  dt_plot2[, c("V1", "V2") := list(scales::rescale(V1, to=c(0+hjust,1+hjust)), scales::rescale(V2, to=c(0+vjust,1+vjust)))]
  if(!color_unmapped){
    dt_plot1[!nhood %in% dt_matches[[1]], (color_by) := NA]
    dt_plot2[!nhood %in% dt_matches[[2]], (color_by) := NA]
  }
  dt_points <- rbindlist(list(milo1=dt_plot1, milo2=dt_plot2), idcol = "milo")
  
  dt_lines <- copy(dt_matches)
  setnames(dt_lines, c("milo1", "milo2"))
  dt_lines <- as.data.table(lapply(dt_lines, as.character))
  
  dt_lines <- merge(dt_lines, dt_points[milo=="milo1", .(nhood, V1, V2)], by.x="milo1", by.y="nhood")
  setnames(dt_lines, c("V1", "V2"), c("xstart", "ystart"))
  dt_lines <- merge(dt_lines, dt_points[milo=="milo2", .(nhood, V1, V2)], by.x="milo2", by.y="nhood")
  setnames(dt_lines, c("V1", "V2"), c("xend", "yend"))
  
  if(size_byNcell){
    p <- ggplot(dt_points, aes(x=V1, y=V2)) +
      geom_point(aes(fill=.data[[color_by]], size=Ncell), color="black", pch=21, alpha=alpha_points) +
      theme_minimal() +
      guides(fill="none", size="none") +
      labs(x="", y="") +
      theme(panel.grid = element_blank(), axis.text.x = element_blank(), axis.text.y=element_blank())
  } else {
    p <- ggplot(dt_points, aes(x=V1, y=V2)) +
      geom_point(aes(fill=.data[[color_by]]), color="black", pch=21, alpha=alpha_points) +
      theme_minimal() +
      guides(fill="none", size="none") +
      labs(x="", y="") +
      theme(panel.grid = element_blank(), axis.text.x = element_blank(), axis.text.y=element_blank())
  }
  
  if(draw_lines){
    p <- p +
      geom_segment(data=dt_lines, aes(x=xstart, y=ystart, xend=xend, yend=yend), alpha = alpha_lines, color="grey50")
  }
  
  if(!is.null(fill_breaks_colors)){
    p <- p + scale_fill_manual(breaks = fill_breaks_colors[[1]], values=fill_breaks_colors[[2]], na.value = alpha("white", alpha=0))
  }
  
  return(p)
}



# Save MiloR to Anndata file with nhoods slot  ---------------------------------
# Overwrite zellkonverter H5AD writer so that it can pass on nhoods
.H5ADwriter_custom <-
  function (sce,
            file,
            X_name,
            skip_assays,
            compression,
            verbose = NULL,
            ...)
  {
    adata <-
      zellkonverter::SCE2AnnData(
        sce,
        X_name = X_name,
        skip_assays = skip_assays,
        verbose = verbose,
        ...
      )
    if (class(sce) == "Milo"){
      adata$obsm <- list(nhoods = nhoods(sce))
      adata$uns <- list(nhoodIndex = unlist(nhoodIndex(sce)))
    }
    zellkonverter:::.ui_step("Writing {.file { .trim_path(file)} }",
                             msg_done = "Wrote {.file { .trim_path(file)} }",
                             spinner = TRUE)
    if (!is.null(compression)) {
      zellkonverter:::.ui_info("Using {.field compression} compression")
    }
    adata$write_h5ad(file, compression = compression)
  }


writeH5AD_custom <-
  function (sce,
            file,
            X_name = NULL,
            skip_assays = FALSE,
            compression = c("none",
                            "gzip", "lzf"),
            version = NULL,
            verbose = NULL,
            ...)
  {
    require(SingleCellExperiment)
    require(zellkonverter)
    require(basilisk)
    compression <- match.arg(compression)
    if (compression == "none") {
      compression <- NULL
    }
    ass_list <- assays(sce)
    is_da <- logical(length(ass_list))
    for (a in seq_along(ass_list)) {
      if (is(ass_list[[a]], "DelayedMatrix") &&
          !is_sparse(ass_list[[a]])) {
        is_da[a] <- TRUE
        assay(sce, a, withDimnames = FALSE) <-
          zellkonverter:::.make_fake_mat(dim(sce))
      }
    }
    env <- zellkonverterAnnDataEnv(version)
    version <- gsub("zellkonverterAnnDataEnv-", "", slot(env,
                                                         "envname"))
    zellkonverter:::.ui_info("Using {.field anndata} version {.field {version}}")
    file <- path.expand(file)
    basiliskRun(
      env = env,
      fun = .H5ADwriter_custom,
      sce = sce,
      file = file,
      X_name = X_name,
      skip_assays = skip_assays,
      compression = compression,
      verbose = verbose,
      ...
    )
    if (any(is_da)) {
      for (p in which(is_da)) {
        if (p == 1L) {
          curp <- "X"
        }
        else {
          curp <- file.path("layers", assayNames(sce)[p])
        }
        rhdf5::h5delete(file, curp)
        mat <- ass_list[[p]]
        if (!is_sparse(mat)) {
          HDF5Array::writeHDF5Array(mat, filepath = file,
                                    name = curp)
        }
        else {
          zellkonverter:::.write_CSR_matrix(file, name = curp, mat = mat)
        }
      }
    }
    invisible(NULL)
  }