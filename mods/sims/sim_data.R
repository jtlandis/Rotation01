

box::use(mods/gr,
         mods/similarity[cluster_data])

#'@export
gr_sub <- gr$GR[grep('LINC01128', gr$GR$gene_names)]

#' @export
.hclust <- cluster_data(gr_sub)

#' @export
tse_sim <- function(ncols = 10){
  TreeSummarizedExperiment::TreeSummarizedExperiment(
    assays = list(counts = matrix(0L, nrow = length(gr_sub), ncol = ncols, dimnames = list(names(gr_sub), paste0("C", 1:ncols)))),
    rowRanges = gr_sub,
    rowTree = ape::as.phylo(.hclust)
  )
}
