#' Plot matches in side-by-side embeddings
#'
#' @param milos List of 2 Milo objects.
#' @param dt_match A data.table with 3 columns containing the nhood-nhood matching and similarities. Typically generated with match_nhoods.
#' @param cols_color A vector of length 2 containing the names of the columns in colData used to color the neighbourhoods. If NULL, no coloring.
#' @param dt_palette NULL or a data.table with 2 columns, one column must be named after cols_color and contains the neighbourhoods labels, the second column (named 'color') contains a HEX color code.
#' @param dimred Which reduced dimensions to use?
#' @param args_process_coordinates A list of 2 lists with preprocessing instructions for the neighbourhoods' coordinates. The nested lists must contain an entry called 'angle' with a scalar, and an entry called 'shift' with a numeric vector of length 2. The shift vector regulates the amount of horizontal and vertical shift between the 2 atlases to avoid overlapping. The angle is useful to rotate the atlases and disentangle matches lines.
#' @param add_non_matched_category not implemented yet
#' @param sim_limits Numeric vector of length 2. Limits the colorscale of the similarity.
#' @param title Title of the plot.
#' @param linewd Numeric, width of the matching lines
#'
#' @returns A ggplot object.
#' @export
#'
#' @examples
plot_matches_embed <- function(milos,
                         dt_match,
                         cols_color = c("celltype", "celltype"),
                         dt_palette = NULL,
                         dimred = "UMAP",
                         args_process_coordinates = list(list(angle = 0, shift = c(0, 0)),
                                                          list(angle = 0, shift = c(10, 0))),
                         add_non_matched_category = c(),
                         sim_limits = NULL,
                         title = "",
                         linewd = 0.35) {
  # Rotate, scale and shift such that can plot lines in-between umap projections
  .process_coordinates <- function(df, angle, shift, s="scale"){

    .rotate_coordinates <- function(df, angle, coord_cols = c("V1", "V2")) {
      # Convert angle to radians
      theta <- angle * pi / 180

      # Create the rotation matrix
      rotation_matrix <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), nrow = 2, byrow = TRUE)

      # Rotate the coordinate columns
      coords <- as.matrix(df[, ..coord_cols])
      rotated_coords <- as.data.table(coords %*% rotation_matrix)
      colnames(rotated_coords) <- coord_cols

      # Combine the rotated coordinates with the other columns
      other_cols <- setdiff(colnames(df), coord_cols)
      rotated_df <- cbind(df[,..other_cols], rotated_coords)

      return(rotated_df)
    }

    df <- copy(df)
    df <- .rotate_coordinates(df, angle=angle)
    if(s=="scale"){
      df[, V1 := scale(V1) + shift[1]]
      df[, V2 := scale(V2)+ shift[2]]
    } else if (s=="01"){
      df[, V1 := (V1-min(V1))/(max(V1)-min(V1)) + shift[1]]
      df[, V2 := (V2-min(V2))/(max(V2)-min(V2)) + shift[2]]
    }

    return(df)
  }

  # Check arguments, set defaults
  if (!is.null(cols_color)) {
    if (length(cols_color) != 2) {
      stop("cols_color must either be NULL or of length 2.")
    } else if (cols_color[1] != cols_color[2]) {
      warning(
        sprintf(
          "cols_color[1] and cols_color[2] are different. This is not yet supported, both values will be set to: %s",
          cols_color[1]
        )
      )
      cols_color[2] <- cols_color[1]
    }
    annot1 <- annotate_nhoods(milos[[1]], cols_color[1])
    annot2 <- annotate_nhoods(milos[[2]], cols_color[2])
  } else {
    # If cols_color is NULL, use a placeholder to have all the points appearing in black
    cols_color <- c("PLACEHOLDER", "PLACEHOLDER")
    annot1 <- data.table(
      Nhood = 1:ncol(nhoods(milos[[1]])),
      Nhood_center = colnames(nhoods(milos[[1]])),
      PLACEHOLDER = rep("_", ncol(nhoods(milos[[1]])))
    )
    annot2 <- data.table(
      Nhood = 1:ncol(nhoods(milos[[2]])),
      Nhood_center = colnames(nhoods(milos[[2]])),
      PLACEHOLDER = rep("_", ncol(nhoods(milos[[2]])))
    )
  }

  # Get centroids
  name1 <- colnames(dt_match)[1]
  name2 <- colnames(dt_match)[2]
  dt_dimred1 <- extract_centroid_nhoods(milos[[1]], dimred, use_indexCell=TRUE)
  dt_dimred2 <- extract_centroid_nhoods(milos[[2]], dimred, use_indexCell=TRUE)

  # Add annotations to matches
  matchtable <- merge(dt_match, annot1, by.x=name1, by.y="Nhood_center")
  matchtable <- merge(matchtable, annot2, by.x=name2, by.y="Nhood_center", suffixes=c(paste0("_", name1), paste0("_", name2)))

  # Get coordinates of paired nhoods merged with corresponding color annotation
  dt_coords1 <- copy(dt_dimred1[nhood %in% matchtable[[name1]]])
  dt_coords2 <- copy(dt_dimred2[nhood %in% matchtable[[name2]]])

  cols <- c(name1, name2, "sim", paste0(cols_color[1], "_", name1), paste0(cols_color[2], "_", name2))
  dt_coords1 <- merge(matchtable[, ..cols], dt_coords1, by.x=name1, by.y="nhood")
  dt_coords2 <- merge(matchtable[, ..cols], dt_coords2, by.x=name2, by.y="nhood")
  dt_points <- rbindlist(list(nhoods_set1=dt_coords1, nhoods_set2=dt_coords2), idcol = "nhoods_set", use.names = TRUE)

  # if(length(add_non_matched_category)>0){
  #   tmp1 <- annot1[celltype %in% add_non_matched_category, .(Nhood_center, celltype)]
  #   setnames(tmp1, "Nhood_center", "nhood")
  #   tmp1[, nhoods_set :=  "nhoods1"]
  #   tmp1 <- merge(tmp1, dt_umap1, by="nhood")
  #
  #   tmp2 <- annot2[celltype %in% add_non_matched_category, .(Nhood_center, celltype)]
  #   setnames(tmp2, "Nhood_center", "nhood")
  #   tmp2[, nhoods_set :=  "nhoods2"]
  #   tmp2 <- merge(tmp2, dt_umap2, by="nhood")
  #
  #   tmp <- rbindlist(list(tmp1, tmp2))
  #   setnames(tmp, "celltype", "value")
  #   matched <- dt_points[, nhood]
  #   tmp <- tmp[!nhood %in% matched]
  #
  #   setcolorder(tmp, colnames(dt_points))
  #   dt_points <- rbindlist(list(dt_points, tmp), use.names = TRUE)
  # }

  # Process coordinates after optionally adding the non-matched categories
  dt_points <- split(dt_points, by="nhoods_set")
  args1 <- args_process_coordinates[[1]]
  args1[["df"]] <- dt_points[["nhoods_set1"]]
  dt_points[["nhoods_set1"]] <- do.call(.process_coordinates, args1)
  args2 <- args_process_coordinates[[2]]
  args2[["df"]] <- dt_points[["nhoods_set2"]]
  dt_points[["nhoods_set2"]] <- do.call(.process_coordinates, args2)
  dt_points <- rbindlist(dt_points)
  dt_points[, unified_color_col := ifelse(nhoods_set=="nhoods_set1", get(paste0(cols_color[1], "_nhoods1")), get(paste0(cols_color[2], "_nhoods2")))]

  # Get coordinates of segments connecting a pair of neighbourhoods
  cols <- c("nhoods_set", name1, name2, "V1", "V2", "sim")
  dt_lines <- split(dt_points[, ..cols], by = "nhoods_set")
  dt_lines <- merge(dt_lines[["nhoods_set1"]], dt_lines[["nhoods_set2"]], by=c(name1, name2, "sim"), suffixes=c("_start", "_end"))

  p <- ggplot(dt_points, aes(x=V1, y=V2)) +
    geom_segment(data=dt_lines, aes(x=V1_start, y=V2_start, xend=V1_end, yend=V2_end, color=sim), alpha=0.5, linewidth=linewd) +
    geom_point(aes(fill=unified_color_col), shape=21) +
    theme_classic() +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank()) +
    scale_x_continuous(limits=c()) +
    scale_y_continuous(limits=c()) +
    labs(subtitle = title, x=paste0(dimred, "1 scaled"), y=paste0(dimred, "2 scaled"))

  if(!is.null(dt_palette)){
    p <- p + scale_fill_manual(breaks = dt_palette[[cols_color[1]]], values=dt_palette$color)
  }
  if(is.null(sim_limits)){
    p <- p + viridis::scale_color_viridis(option = "B")
  } else {
    p <- p + viridis::scale_color_viridis(option = "B", limits=sim_limits)
  }
  return(p)
}
