

box::use(
  ./sims[...],
  ./tree[...]
)

dtse <- shift_mu(tse, delta = 2.5)
dtse <- sim_counts(dtse)

split_counts <- function(se) {
  
  sims <- cnts <- counts(se)
  sims[] <- rbinom(prod(dim(se)), size = cnts, prob = 0.5)
  assay(se, "train") <- sims
  assay(se, "test") <- cnts - sims
  se
}
dtse <- split_counts(dtse)
train <- dtse
assay(train, "counts") <- assay(train, "train")
train <- tse_by_nodes(train)
train_res <- deseq(train, ~ groups, size_factors = 1L)
train_res <- train_res[[2]]
train_res_climbed <- climb(train_res, "order")
plot_results(climb(train_res, "order"))
qcolSum <- function(x, na.rm = FALSE) {
  n <- nrow(x)
  dn <- ncol(x)
  .Internal(colSums(x, n, dn, na.rm))
}
tse_by_climb <- function(tse, climb_res) {
  
  
  keep <- subset(climb_res, keep)
  n_col <- ncol(tse)
  col_nms <- colnames(tse)
  dimnms <- list(keep$nodes, col_nms)
  lst_mat <- vector('list', nrow(keep))
  Assays <- assays(tse)
  NewAssays <- vector('list', length(Assays))
  i_seq <- seq_along(lst_mat)
  desc <- keep$descendants
  node_name <- keep$nodes
  for (j in seq_along(NewAssays)) {
    assay_j <- Assays[[j]]
    for (i in i_seq) {
      lst_mat[[i]] <- matrix(
        data = qcolSum(assay_j[desc[[i]],]),
        nrow = 1L,
        ncol = n_col,
        dimnames = list(node_name[i], col_nms)
      )
    }
    NewAssays[[j]] <- do.call('rbind', lst_mat)
  }
  names(NewAssays) <- names(Assays)
  
  SummarizedExperiment::SummarizedExperiment(
    assays = NewAssays,
    rowData = methods::as(keep, "DataFrame"),
    colData = SummarizedExperiment::colData(tse),
    metadata = list(
      TreeSummarizedExperiment = tse
    )
  )
  
  
}