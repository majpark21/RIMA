suppressPackageStartupMessages({
  library(argparser)
  library(data.table)
  library(miloR)
  library(SingleCellExperiment)
  library(stringr)
  library(igraph)
  source("./utils.R")
})

t0 <- proc.time()
parser <- argparser::arg_parser(
  description = "This script takes 2 Milo objects, with populated nhoods slot,
  and matches their neighborhoods. Neighborhoods are matched using a maximum bipartite
  matching procedure, based on pairwise neighborhood similarities (correlation).
  Matching is ensured to map relevant pairs of neighborhoods by first filtering
  insignificant nhood-nhood edges thanks to a permutation test against
  unspecific neighborhoods which are formed by scrambling cell identities.
  
  The milo1 object is used as the reference object for which each nhood is
  independently tested against scrambled versions of milo2's nhoods; a 
  statistical cutoff is drawn from this test and only relevant edges in the
  direction milo1 -> milo2 are passed as potential matches to the bipartite matching procedure.",
  name = "Cross-atlas neighborhood matching.",
  hide.opts = TRUE
)

parser <- add_argument(
  parser = parser,
  arg = "milo1",
  help = "Path to RDS file with 1st milo object. Its @nhoods slot must be populated.
  Rownames (i.e. features/gene names) must be specified to ensure that the right features
  are compared; feature mapping can be provided, see --featuretable.
  During filtering of insignificant nhood-nhood edges, the nhoods of milo1 are kept
  intact and for each of them, a statistical similarity cutoff is determined by
  comparing them against scrambled versions of nhoods in milo2.",
  default = NULL,
  type = "character",
  nargs = 1
)

parser <- add_argument(
  parser = parser,
  arg = "milo2",
  help = "Path to RDS file with 1st milo object. Its @nhoods slot must be populated.
  Rownames (i.e. features/gene names) must be specified to ensure that the right features
  are compared; feature mapping can be provided, see --featuretable.
  During filtering of insignificant nhood-nhood edges, the nhoods of milo1 are kept
  intact and for each of them, a statistical similarity cutoff is determined by
  comparing them against scrambled versions of nhoods in milo2.",
  default = NULL,
  type = "character",
  nargs = 1 
)

parser <- add_argument(
  parser = parser,
  arg = "--similarity",
  help = 'Method to calculate pairwise nhood similarity. Must be one of c("pearson", "kendall", "spearman").',
  default = "spearman",
  type = "character",
  nargs = 1
)

parser <- add_argument(
  parser = parser,
  arg = "--nscrambles",
  help = "Number of cell identity scrambling used to determine nhood-nhood significance.",
  default = 5,
  type = "numeric",
  nargs = 1
)


parser <- add_argument(
  parser = parser,
  arg = "--adjust",
  help = "Pvalue adjustement before calling nhood-nhood edges. Any valid method
  option from 'p.adjust'. Default to Holm for FWER correction.",
  default = "holm",
  type = "character",
  nargs = 1
)


parser <- add_argument(
  parser = parser,
  arg = "--adjustalpha",
  help = "Level of significance used to call statistically significant nhood-nhood edge.
  This is applied to the adjusted p-values derived from the permutation test.",
  default = 0.05,
  type = "numeric",
  nargs = 1
)

parser <- add_argument(
  parser = parser,
  arg = "--labelscramble",
  help = "Column in 'colData(milo2)' used to weight cell resampling in nhood scrambling.
  This is to ensure that scrambled nhoods are unspecific do not look like the
  most abundant label. If 'false', permute cells without resampling.",
  default = "celltype",
  type = "character",
  nargs = 1
)

parser <- add_argument(
  parser = parser,
  arg = "--assay",
  help = "Milos' assay used to perform nhoods comparison.",
  default = "logcounts",
  type = "character",
  nargs = 1
)

parser <- add_argument(
  parser = parser,
  arg = "--featuretable",
  help = "Path to a csv file (or other formats recognized by data.table::fread)
  with 2 columns that maps features across both milos in case the features are
  not fully overlapping.
  The 1st (resp. 2nd) column refers to rownames in milo1 (resp. milo2).
  This is useful for example to compare atlases of different species.
  This option can also be used to perform the matching on a subset of features.
  If not specified, rownames must fully match between both milos.",
  default = NA,
  type = "character",
  nargs = 1
)

parser <- add_argument(
  parser = parser,
  arg = "--signiftable",
  help = "Path to a CSV file, such as the one produced with the --nomatching option.
  Used to provide edge-edge significance table and go straight to matching.",
  default = NA,
  type = "character",
  nargs=1
)

parser <- add_argument(
  parser = parser,
  arg = "--outdir",
  help = "Path of the directory where to write the matching table.",
  default = '.',
  type = "character",
  nargs = 1
)

# parser <- add_argument(
#   parser = parser,
#   arg = "--ndigits",
#   help = "Round up similarity and p-values to avoid large file.",
#   default = 4,
#   type = "numeric",
#   nargs = 1
# )

parser <- add_argument(
  parser = parser,
  arg = "--pval",
  help = "Flag, whether to return p-values associated with each nhood-nhood edge. [default: Unless this flag is specified, p-values are not returned.]",
  flag = TRUE
)

parser <- add_argument(
  parser = parser,
  arg = "--nomatching",
  help = "Flag, whether to stop before the matching step and return the edge-edge similarity values along with statistical significance. [default: Unless this flag is specified, matching will be performed.]",
  flag = TRUE
)

parser <- add_argument(
  parser = parser,
  arg = "--direction",
  help= "Direction of the matching.
  By default ('lr'; short for left-right), nhoods in milo1 (i.e. the left object) are kept intact and nhoods in milo2 (i.e. the right object) are scrambled to derive the significance threshold before the 1:1 matching.
  If 'rl' (right-left), this is inverted.
  If 'b' (bidirectional), the nhood-nhood insignificant edges trimming is performed in both direction and the intersection of significant nhood-nhood edges are passed for matching.
  
  The 'b' option is the most conservative and slower but accounts for the difference of heterogeneity between the 2 milos.",
  default = 'lr',
  type = "character",
  nargs = 1
)

# Parsing parameters ----
cat("Parsing arguments...\n")
argv <- parse_args(parser)

fi_milo1 <- argv$milo1
fi_milo2 <- argv$milo2
fi_features <- argv$featuretable
out_dir <- argv$outdir

sel_assay <- argv$assay
sim_measure <- argv$similarity
n_scrambles <- argv$nscrambles
method_adjust <- argv$adjust
alpha_adjust <- argv$adjustalpha
col_scramble_label <- argv$labelscramble

if(!sim_measure %in% c("pearson", "kendall", "spearman")){
  stop("'--similarity' must be one of: 'pearson', 'kendall', 'spearman'.")
}
if((alpha_adjust < 0) | (alpha_adjust > 1)){
  stop("'--alphaadjust' must be between 0 and 1.")
}
if(!argv$direction %in% c('lr', 'rl', 'b')){
  stop("--direction must be one of: 'lr', 'rl', 'b'.")
}

# Load data ----
cat("Loading data...\n")
# fi_genes <- "./genes_geneBasis.rds"

name1 <- str_replace(basename(fi_milo1), "\\.rds$", "")
name2 <- str_replace(basename(fi_milo2), "\\.rds$", "")

if(!dir.exists(out_dir)){dir.create(out_dir)}
out_fi <- sprintf("wbgm__%s__%s%s__%sscrambles__%s__%s.csv", sim_measure, str_to_upper(method_adjust), str_replace(alpha_adjust, "\\.", "_"), n_scrambles, name1, name2)
out_fi <- file.path(out_dir, out_fi)

milo1 <- readRDS(fi_milo1)
milo2 <- readRDS(fi_milo2)

# Subset to/reorder genes of interest ----
cat("Subset to/reorder genes of interest...\n")
if(is.na(fi_features)){
  if(any(sort(unique(rownames(milo1))) != sort(unique(rownames(milo2))))){
    stop("milo1 and milo2 must have the same rownames, or --featuretable must be specified.")
  }
  # Reorder features to same order
  milo2 <- milo2[rownames(milo1),]
} else {
  dt_features <- fread(fi_features)
  milo1 <- milo1[dt_features[[1]], ]
  milo2 <- milo2[dt_features[[2]], ]
  # Rename features in milo2 to their milo1's counterpart
  rownames(milo2) <- dt_features[[1]]
}

# Calculate nhoods similarities ----
if(is.na(argv$signiftable)){  # not needed if already got pvalues
  cat("Calculate nhoods similarities...\n")
  if(all(dim(milo1@nhoodExpression)==c(1,1))){
    milo1 <- miloR::calcNhoodExpression(milo1, sel_assay, subset.row = NULL)
  }
  if(all(dim(milo2@nhoodExpression)==c(1,1))){
    milo2 <- miloR::calcNhoodExpression(milo2, sel_assay, subset.row = NULL)
  }
  
  dt_sims <- getSims(
    milo1,
    milo2,
    subset_row1 = NULL,
    subset_row2 = NULL,
    sim = sim_measure
  )
  setnames(dt_sims, c(name1, name2, "sim"))
}

# Filter out insignificant edges ----
## Function definition ----

# milo: a milo object with filled nhoods slot
# n: number of permutations
# replace: is the resampling to be done with replacement?
# weight_sample: a vector of same length as number of cells that defines how likely a cell is to be picked. If specified, sampling is always done with replacement.
# return milo object with permuted cells
scramble_nhood_avg <- function(milo, n=5, replace=FALSE, weight_sample=NULL, assay="logcounts"){
  # Permute observations between cells and keep nhood structure as is to scramble them
  mat <- assay(milo, assay)
  all_cells <- colnames(milo)
  if(!is.null(weight_sample)){replace <- TRUE}
  
  l_out <- list()
  for(ii in 1:n){
    permut_cells <- sample(
      all_cells,
      size = length(all_cells),
      replace = replace,
      prob = weight_sample
    )
    permut_mat <- mat[, permut_cells]
    colnames(permut_mat) <- all_cells
    assay(milo, paste0("permut_", assay)) <- permut_mat
    milo <- calcNhoodExpression(milo, assay = paste0("permut_", assay))
    l_out[[ii]] <- nhoodExpression(milo)
  }
  return(l_out)
}

# Take 2 dts of matches with 3 columns: one corresponds to nhoods names in target species, one corresponds to nhood names in source species, last one corresponds to similarity values
# true is the authentic similarity values
# scrambled forms the null population against which true is compared
# alpha is level of significance
# adjust is passed to p.adjust to correct multiple testing. If NULL does not do pvalue correction
# col_group indicates at which group level the statistical tests should be performed; typically used to indicate that a different stats-significant cutoff is calculated for each nhood (i.e. group level); typically takes the column of target nhood to filter significant edges at nhood level. If NULL, significance is calculated from pooled.
# return_all==T means to return the original dt_sim_true with 2 additional columns indicating whether the hit is sgnificant and its associated pvalue
filter_signifMatches <- function(dt_sim_true, dt_sim_scrambled, alpha=0.05, adjust="fdr", col_sim="sim", col_group=NULL, return_notSignif=TRUE){
  # Basic function applied to each group from which to extract pvalues for each edge, or the significance cutoff
  .pvals_signifMatches <- function(v_sim_true, v_sim_scrambled, alpha=0.05, adjust="fdr", return_cutoff=FALSE){
    # Special case speed-up: no need to calculate costly ecdf
    if((return_cutoff) & (is.null(adjust))){
      cutoff <- quantile(v_sim_scrambled, probs=1-alpha)
      return(cutoff)
    }
    # All other cases need to calculate ecdf
    fn_ecdf <- ecdf(v_sim_scrambled)  # Estimate p-values function using scrambled distribution
    pvals <- 1 - fn_ecdf(v_sim_true)  # Assign p-values to authentic distribution
    if(!is.null(adjust)){
      pvals_adj <- p.adjust(pvals, method=adjust)
    } else {
      pvals_adj <- pvals
    }
    if(return_cutoff){
      signif_sim <- which(pvals <= alpha)
      cutoff <- min(v_sim_true[signif_sim])
      return(cutoff)
    } else {
      return(list(pvals=pvals, pvals_adj=pvals_adj))
    }
  }
  cols_original <- colnames(dt_sim_true)
  dt_sim_true <- copy(dt_sim_true)  # Avoid modifying original object
  dt_sim_scrambled <- copy(dt_sim_scrambled)
  
  # col_out <- ifelse(is.null(adjust), "pval", "pval_adjusted")
  col_pval <- "pval"
  col_pval_adj <- "pval_adjusted"
  col_nhoods <- setdiff(colnames(dt_sim_true), col_sim)
  # Make sure nhoods are not integer because confuses grouping
  for(col in col_nhoods){
    dt_sim_true[[col]] <- as.character(dt_sim_true[[col]])
    dt_sim_scrambled[[col]] <- as.character(dt_sim_scrambled[[col]])
  }
  
  if(is.null(col_group)){
    tmp <- .pvals_signifMatches(
      v_sim_true = dt_sim_true[[col_sim]],
      v_sim_scrambled = dt_sim_scrambled[[col_sim]],
      alpha = alpha,
      adjust = adjust,
      return_cutoff = FALSE
    )
    dt_sim_true[[col_pval]] <- tmp$pvals
    dt_sim_true[[col_pval_adj]] <- tmp$pvals_adj
    # dt_sim_true[[col_out]] <- .pvals_signifMatches(
    #   v_sim_true = dt_sim_true[[col_sim]],
    #   v_sim_scrambled = dt_sim_scrambled[[col_sim]],
    #   alpha = alpha,
    #   adjust = adjust,
    #   return_cutoff = FALSE
    # )
  } else {
    setkeyv(dt_sim_true, col_group)
    setkeyv(dt_sim_scrambled, col_group)
    # dt_sim_true[, (col_out) := -1.0]
    groups <- unique(dt_sim_true[[col_group]])
    pb <- txtProgressBar(min = 0, max = length(groups), style = 3)
    for(idx_groups in seq(length(groups))){
      curr_group <- groups[idx_groups]
      tmp <- .pvals_signifMatches(
        v_sim_true = dt_sim_true[curr_group, get(col_sim)],
        v_sim_scrambled = dt_sim_scrambled[curr_group, get(col_sim)],
        alpha = alpha,
        adjust = adjust,
        return_cutoff = FALSE
      )
      dt_sim_true[curr_group, (col_pval) := tmp$pvals]
      dt_sim_true[curr_group, (col_pval_adj) := tmp$pvals_adj]
      # dt_sim_true[curr_group, (col_out) := .pvals_signifMatches(
      #   v_sim_true = dt_sim_true[curr_group, get(col_sim)],
      #   v_sim_scrambled = dt_sim_scrambled[curr_group, get(col_sim)],
      #   alpha = alpha,
      #   adjust = adjust,
      #   return_cutoff = FALSE
      # )]
      setTxtProgressBar(pb, idx_groups)
    }
    close(pb)
  }
  
  dt_sim_true[, is_significant := pval_adjusted <= alpha]
  if(!return_notSignif){
    dt_sim_true <- dt_sim_true[is_significant == TRUE]
  }
  
  return(dt_sim_true)
}

getWeightScrambling <- function(milo, col_scramble){
  if(col_scramble=="false"){
    weight_cells <- NULL
  } else {
    weight_cells <- 1/table(colData(milo)[[col_scramble]])
    weight_cells <- weight_cells / sum(weight_cells)
    weight_cells <- weight_cells[milo[[col_scramble]]]
  }
  return(weight_cells)
}

# Simes pvalue combination special case with 2 values, speed up and easier in data.table
combine_2pvals_sime <- function(pval1, pval2){
  pvals <- c(pval1, pval2)
  ordered_pvals <- sort(pvals)
  adj_pvals <- 2 *(ordered_pvals/(1:2)) 
  return(min(1, adj_pvals))
}

pipelineEdgeEdgeSignif <- function(milo_intact, milo_scramble, name_intact, name_scramble, col_scramble_label, n_scrambles, sel_assay, sim_measure, dt_sims, adjust, alpha_adjust){
  ## Generate null (scrambled) distribution of similarities, upsample rare celltypes ----
  cat("Generate scrambled neighborhoods...\n")
  l_sims_scrambled <- list()
  
  # Weight is inversely proportional to the number of cells of the corresponding celltype
  weight_cells <- getWeightScrambling(milo = milo_scramble, col_scramble = col_scramble_label)
  # Generate scrambled nhood average expressions
  l_scrambled_assay_nhoods <- scramble_nhood_avg(
    milo_scramble,  # already subset here for speed up
    replace = TRUE,
    n = n_scrambles,
    weight_sample = weight_cells,
    assay=sel_assay
  )
  
  ## Calculate similarities against scrambled nhood average expression ----
  cat("Calculate similarities against scrambled nhood average expression...\n")
  for(ii in 1:length(l_scrambled_assay_nhoods)){
    dt_sim <- getSims(
      milo1 = milo_intact,
      milo2 = l_scrambled_assay_nhoods[[ii]],
      subset_row1 = NULL,
      subset_row2 = NULL,
      sim = sim_measure
    )
    
    setnames(dt_sim, c(name_intact, name_scramble, "sim"))
    l_sims_scrambled[[ii]] <- dt_sim
  }
  l_sims_scrambled <- rbindlist(l_sims_scrambled, idcol="iteration")
  
  ## Filter insignificant nhood-nhood edges ----
  cat("Filter insignificant nhood-nhood edges...\n")
  dt_sims_withSignif <- filter_signifMatches(
    dt_sim_true = dt_sims,
    dt_sim_scrambled = l_sims_scrambled,
    col_group = name_intact,
    alpha = alpha_adjust,
    adjust = adjust,
    col_sim = "sim",
    return_notSignif = TRUE
  )
  
  return(dt_sims_withSignif)
}

## Scramble nhoods and calculate edge-egde pval ----
if(is.na(argv$signiftable)){
  ls_args_common <- list(
    col_scramble_label = col_scramble_label,
    n_scrambles = n_scrambles,
    sel_assay = sel_assay,
    sim_measure =  sim_measure,
    adjust = method_adjust,
    alpha_adjust = alpha_adjust
  )
  if(argv$direction=='lr'){
    cols_sim <- c(name1, name2, "sim")
    ls_args <- list(
      milo_intact = milo1,
      milo_scramble = milo2,
      name_intact = name1,
      name_scramble = name2,
      dt_sims = dt_sims[, ..cols_sim]
    )
    ls_args <- append(ls_args, ls_args_common)
    dt_sims_withSignif <- do.call(pipelineEdgeEdgeSignif, ls_args)
  } else if(argv$direction=='rl'){
    cols_sim <- c(name2, name1, "sim")
    ls_args <- list(
      milo_intact = milo2,
      milo_scramble = milo1,
      name_intact = name2,
      name_scramble = name1,
      dt_sims = dt_sims[, ..cols_sim]
    )
    ls_args <- append(ls_args, ls_args_common)
    dt_sims_withSignif <- do.call(pipelineEdgeEdgeSignif, ls_args)
  } else if(argv$direction=='b'){
    ls_pvals <- list()
    # left -> right
    cols_sim <- c(name1, name2, "sim")
    ls_args <- list(
      milo_intact = milo1,
      milo_scramble = milo2,
      name_intact = name1,
      name_scramble = name2,
      dt_sims = dt_sims[, ..cols_sim]
    )
    ls_args <- append(ls_args, ls_args_common)
    ls_pvals[[1]] <- do.call(pipelineEdgeEdgeSignif, ls_args)
    # right -> left
    cols_sim <- c(name2, name1, "sim")
    ls_args <- list(
      milo_intact = milo2,
      milo_scramble = milo1,
      name_intact = name2,
      name_scramble = name1,
      dt_sims = dt_sims[, ..cols_sim]
    )
    ls_args <- append(ls_args, ls_args_common)
    ls_pvals[[2]] <- do.call(pipelineEdgeEdgeSignif, ls_args)
    
    
    # Return p-values from both directions; combine with Simes ,ethod and adjust the combined Pvalues
    # Discard individually adjusted pvalues, because adjustement is made on combined
    cols <- c(name1, name2, "pval")
    dt_sims_withSignif <- merge(ls_pvals[[1]], ls_pvals[[2]], by=c(name1, name2), suffixes=c("_lr", "_rl"))
    cat("Combining left-right and right-left pvalues...")
    dt_sims_withSignif[, pval_combined := combine_2pvals_sime(pval_lr, pval_rl), by=1:nrow(dt_sims_withSignif)]
    dt_sims_withSignif[, pval_combined_adjusted := p.adjust(pval_combined, method=method_adjust)]
    dt_sims_withSignif[, is_significant := pval_combined_adjusted <= alpha_adjust]
    
    
    # Consistent column names and order with unidirectional
    dt_sims_withSignif[, sim := sim_lr]
    setcolorder(dt_sims_withSignif, c(name1, name2, "sim"))
    
    # # Keep the intersection of significant edges -  Low power
    # cols <- c(name1, name2, "sim")
    # dt_sims_withSignif <- rbindlist(list(
    #   ls_pvals[[1]][is_significant==T, ..cols],
    #   ls_pvals[[2]][is_significant==T, ..cols]
    # ))
    # dt_sims_withSignif <- dt_sims_withSignif[duplicated(dt_sims_withSignif)]
    # # TODO: look into combining both directions p-val
    # dt_sims_withSignif[, pval_adjusted := NA]
    # dt_sims_withSignif[, is_significant := T]
  }
} else {
  dt_sims_withSignif <- fread(argv$signiftable)
  # Make sure that the columns containing the nhood names were not read as numerical
  cols <- colnames(dt_sims_withSignif)[1:2]
  dt_sims_withSignif[, (cols) := lapply(.SD, as.character), .SDcols = cols]
}

# Perform bipartite matching with significant edges only ----
if(!argv$nomatching){
  cat("Perform bipartite matching with significant edges...\n")
  cols_names<- c(name1, name2)
  cols_names_sim <- c(cols_names, "sim")
  dt_match <- getBipartiteMatch(
    dt_mnn = dt_sims_withSignif[is_significant == T, ..cols_names],
    dt_sim = dt_sims_withSignif[is_significant == T, ..cols_names_sim],
    return_matchObject = F
  )
  
  if(argv$pval){
    dt_match <- merge(
      dt_match,
      dt_sims_withSignif[is_significant == T],
      by = cols_names
    )
    # dt_match[, sim := round(sim, argv$ndigits)]
    # dt_match[, pval_adjusted:= round(pval_adjusted, argv$ndigits)]
  }
  
  cat(paste0("Exporting matching results to:", out_fi))
  fwrite(dt_match, out_fi)
} else {
  cat(paste0("Exporting edge-significance results to:", out_fi))
  # dt_sims_withSignif[, sim := round(sim, argv$ndigits)]
  # dt_sims_withSignif[, pval_adjusted:= round(pval_adjusted, argv$ndigits)]
  fwrite(dt_sims_withSignif, out_fi)
}

data.table::timetaken(t0)
