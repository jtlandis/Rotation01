#' Useful Functions to run simulations on data
#' @name .__module__.
#' @description 
#' helpful functions to perform simulations of rnaseq and evaluate
#' tree climbinb algorithms.
#' 
#' @examples 
#' 
#' set.seed(123)
#' clust <- rnbinom(100, size = 100, mu = sample(2:18, 100, T)) |> dist() |> hclust()
#' phylo <- ape::as.phylo(clust)
#' tse <- TreeSummarizedExperiment(
#'   assays = list(counts = matrix(0, nrow = 100, ncol = 50,
#'                                 dimnames = list(phylo$tip.label, sprintf('cell%02i',1:50)))),
#'   rowTree = phylo)
#' 
#' # add simulation means
#' # creates 'mu' assay and initializes all values with
#' # 10
#' tse <- sim_mu(tse, delta = 10)
#' # add random variation
#' tse <- sim_mu_random(tse, iterations = 100,
#'                      delta_range = c(-1.5, 1.5))
#' # add delta to row nodes
#' #
#' tse <- sim_mu_node(tse, affected_nodes = c(180, 181),
#'                    affected_cols = 26:50,
#'                    delta = 2)
#' 
#' tse <- sim_counts(tse)
#' ntse <- tse_by_nodes(tse)
#' 
#' # Note, before using DESeq2, some colData field
#' # should be set to perform a Differential expression
#' # analysis
#' 
#' ntse$group1 <- assign_group_members(
#'   ncol(ntse),
#'   ngroups = 2
#' ) |> as.factor()
#' 
#' # we can also make random groups
#' # assigning groups may be helpful
#' # earlier when choosing which columns
#' # should be perturbed
#' ntse$group2 <- assign_group_members(
#'   ncol(ntse), ngroups = 4, random = TRUE
#' )
#' 
#' out <- deseq(ntse, ~ group1, size_factors = 1L)
#' plot_results(out[[2]])
#' climb(out[[2]],
#'       unlist(descendants_of_nodes(c(180, 181),
#'                                   clust)),
#'       method = "ord", merge_alpha = 0.05) |>
#'   summarise_climb()
#' 
#' @export
box::use(
  ./sim_utils[...],
  ./sim_data[...]
)