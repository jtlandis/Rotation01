# Simulating Reads


## Setup ---- 

### Packages ----
box::use(mods/gtf[...],
         mods/rnaseq[...],
         mods/similarity[...],
         tibble[tibble, as_tibble],
         mods/sims/sim_utils[...],
         mods/stamps[...],
         mods/tree[find_common_parent_node_merge_index, descendants_of_merge_index,
                   tse_by_nodes],
         ggplot2[...])
# library(TreeSummarizedExperiment)
# library(ggplot2)

### Data ----
gtf <- read_gtf("mods/gencode.v44.annotation.gtf", nrows = 1000)
tbl <- subset(gtf, n_trans==44)$trans[[1]]
tbl$trans_name <- extract_name(tbl$attribute, "transcript")


nIsoform <- nrow(tbl)
nSamps <- 10
.hclust <- cluster_data(tbl, similarity_method = "so", name = "trans")
phylo <- ape::as.phylo(.hclust)
mat <- matrix(0, nrow = nIsoform, ncol = nSamps,
              dimnames = list(tbl$trans_name,
                              sprintf(paste0("C%0",floor(log10(nSamps)) + 1L ,"i"),
                                      1:nSamps)))
tse <- TreeSummarizedExperiment(
  assays = list(counts = mat),
  rowRanges = GenomicRanges::GRanges(tbl$seqname,
                      ranges = IRanges::IRanges(start = tbl$start,
                                       end = tbl$end,
                                       name = tbl$trans_name),
                      strand = tbl$strand),
  rowTree = phylo
)


### Paramters ----

# These parameters are very simple conditions

#### RowData sim ----

# The 'affected_range' is manually chosen according to the visual clustering
# 
# In the future, maybe we choose random nodes to up/down regulate
row_mod <- descendants_of_nodes(77, .hclust)[[1]]


#### ColData sim ----

# Number of groups to split samples into
ngroups <- 2
colData(tse)$groups <- assign_group_members(nSamps = nSamps, ngroups = ngroups)

tse <- sim_mu(tse, delta = 10)

### Functions ----
shift_fun <- function(affected_rows, affected_cols) {
  force(affected_rows)
  force(affected_cols)
  function(tse, delta) {
    sim_mu(tse, affected_rows, affected_cols, delta = delta)
  }
}
shift_mu <- shift_fun(affected_rows = row_mod, affected_cols = 6:10)


simCountsEvalDeseq <- function(tse, ...) {

  tse <- sim_counts(tse)
  ntse <- tse_by_nodes(tse)
  deseq(ntse, ~ groups, size_factors = 1L, ...)
}

do_sim <- function(tse, delta, n_sims = 10, ...) {
  
  tses <- lapply(delta, shift_mu, tse = tse)
  delta <- rep(delta, each = n_sims)
  tses <- rep(tses, each = n_sims)
  
  out <- purrr::map(tses, simCountsEvalDeseq, .progress = TRUE, ...)
  out <- purrr::map(out, `[[`, 2L)
  out <- purrr::map2(out, delta, \(x, y, z) {x[[z]] <- y; x}, z = 'delta')
  out
  
}

do_climb <- function(.data, method = c("strict", "ordered_pval"), note = "none", ..., affected_rows = row_mod) {
  method <- match.arg(method, c("strict", "ordered_pval"))
  args <- list(...)
  if (!is.list(.data))
    stop("`.data` should either be alist of `dds_node_res` or a `dds_node_res` object")
  
  #check if it is
  if (is.object(.data)) {
    #it is classed
    if (!inherits(.data, "dds_node_res"))
      stop("`.data` is not a `dds_node_res` object")
    .data <- list(.data)
  }
  
  .data <- lapply(.data, climb, method = method, ...)
  out <- lapply(.data, summarise_climb, affected_rows = affected_rows, delta = unique(delta))
  out <- dplyr::bind_rows(out)
  out$note <- note
  out$.data <- .data
  out$.method <- method
  out$.args <- list(args)
  out$.sim_number <- 1:nrow(out)
  out
} 

plot_summary <- function(climb_res, annotate_missing = TRUE) {
  p <- ggplot(climb_res, aes(x = FDR, y = sens, group = .sim_number)) +
    geom_point(data = ~subset(.x, !is.nan(FDR)),
               color = "steelblue", fill = "steelblue") +
    facet_wrap(~delta) +
    theme_bw() 
  
  if (annotate_missing) {
    p <- p + geom_label(
      data = ~subset(.x, is.nan(FDR)) |>
        dplyr::summarise(nexclude = dplyr::n(), .by = c(delta, note)),
      aes(x = 0, y = -.05, label = paste0('excluded: ',round(nexclude, 0)),
          group = delta), inherit.aes = F,
      hjust = 0, vjust = 1) +
      ylim(c(-.2, 1))
  }
  p
  
}

## MAIN ----

# out <- do_sim(tse, seq(.5, 6, by = .5), n_sims = 50)
# saveRDS(out, 'rds_cache/sim-50-out.rds')
out <- readRDS("rds_cache/sim-50-out.rds")

method1 <- do_climb(out, method = "strict", note = "Bottom-up") 
method2 <- do_climb(out, method = "ord", note = "Top-Down")
method3 <- do_climb(out, method = "ord", merge_alpha = 0.1, note = "Top-Down: merge_alpha=0.1") 
# method4 <- do_climb(out, method = "ord", merge_alpha = 0.05, note = "Top-Down: merge_alpha=0.05") 

plot_summary(method1)
plot_summary(method2)
.p <- plot_summary(method3, annotate_missing = F) +
  labs(y = "Sensitivity",
       title = expression("50 Similations of " ~ delta))
ggsave("figures/Sensitivity-FDR-delta-sims.png", plot = .p, dpi = 600)
# .p1 <- plot_pvals(dplyr::bind_rows(method2$.data), sig_nodes_only = T) +
#   facet_wrap(~delta, scales = "free_y") +
#   labs(title = "Only Retained Nodes")
# 
# .p2 <- plot_pvals(dplyr::bind_rows(method2$.data), sig_nodes_only = F) +
#   facet_wrap(~delta, scales = "free_y") +
#   labs(title = "All tested Nodes")
# 
# ggsave("/figures/pvals-hist-topdown-retainedNodes.png", .p1)
# ggsave("/figures/pvals-hist-topdown-AllNodes.png", .p2)

#### gganimate ----

# anim <- dplyr::bind_rows(method1, method2, method3) |>
#   plot_summary() +
#   gganimate::transition_states(note)  +
#   ggtitle("{closest_state}")  +
#   gganimate::enter_grow() +
#   gganimate::exit_recolor(color = "red") +
#   gganimate::exit_fade() +
#   gganimate::exit_shrink()
# 
# .anim <- gganimate::animate(
#   anim,
#   renderer = gganimate::gifski_renderer(file = "gifs/sim-50.gif"))

### Split Counts ----

# delta >2 seemed like a high sensitivity and low FDR threshold

dtse <- shift_mu(tse, delta = 2.5)
dtse <- sim_counts(dtse)

split_counts <- function(se) {
  
  sims <- cnts <- assay(se, 'counts')
  sims[] <- rbinom(prod(dim(se)), size = cnts, prob = 0.5)
  test <- train <- se
  assay(train, 'counts') <- sims
  assay(test, 'counts') <- cnts - sims
  test@metadata$original <- se
  list(train = train, test = test)
}

do_splitCounts <- function(tse, delta, n_sims = 10, ..., method = c("ordered_pval", "strict")) {
  box::use(purrr[map, map2])
  method <- match.arg(method, choices = c("ordered_pval", "strict"))
  tses <- lapply(delta, shift_mu, tse = tse)
  delta <- rep(delta, each = n_sims)
  tses <- rep(tses, each = n_sims)
  
  #sim counts
  tses <- map(tses, sim_counts)
  splits <- map(tses, split_counts)
  tests <- map(splits, `[[`, 2L)
  train <- map(splits, `[[`, 1L)
  trained <- map2(train, 
                  delta,
                 function(x, d, ...) {
                   ntse <- tse_by_nodes(x)
                   out <- deseq(ntse, ~ groups, size_factors = 1L, ...)[[2L]]
                   out[['delta']] <- d
                   out
                 }, ..., .progress = "train-data")
  
  trained <- do_climb(trained, method = method, note = "split_counts")
  
  trans_test <- map2(tests, trained$.data, tse_by_climb, .progress = "transform-test-data")
  trans_tested <- map2(trans_test,
                       delta,
                       purrr::safely(
                         function(x, d, ...) {
                           suppressMessages({
                             dds <- DESeqDataSet(x, design = ~ groups)
                             sizeFactors(dds) <- 1L
                             dds <- DESeq(dds, quiet = T, ...)
                             res_names <- resultsNames(dds)
                             out <- lapply(res_names, function(x) results(object = dds, name = x))
                           })
                           out <- out[[2]]
                           out[['delta']] <- d
                           out
                         }, otherwise = "FAILED"), ..., .progress = "test-data")
  
  list(trained_res = trained, tested = trans_tested)
  
  
  
  
}


# tse_split <- split_counts(dtse)
# train <- tse_by_nodes(tse_split$train)
# train_res <- deseq(train, ~ groups, size_factors = 1L)
# train_res <- train_res[[2]]
# train_res_climbed <- climb(train_res, "order")
# plot_results(climb(train_res, "order"))
qcolSum <- function(x, na.rm = FALSE) {
  n <- nrow(x)
  dn <- ncol(x)
  .Internal(colSums(x, n, prod(dn), na.rm))
}
tse_by_climb <- function(tse, climb_res) {
  

  keep <- subset(climb_res, keep)
  n_col <- ncol(tse)
  col_nms <- colnames(tse)
  dimnms <- list(keep$nodes, col_nms)
  lst_mat <- vector('list', nrow(keep))
  Assays <- SummarizedExperiment::assays(tse)
  NewAssays <- vector('list', length(Assays))
  i_seq <- seq_along(lst_mat)
  desc <- keep$descendants
  node_name <- keep$nodes
  for (j in seq_along(NewAssays)) {
    assay_j <- Assays[[j]]
    for (i in i_seq) {
      lst_mat[[i]] <- matrix(
        data = qcolSum(assay_j[desc[[i]],, drop = FALSE]),
        nrow = 1L,
        ncol = n_col,
        dimnames = list(node_name[i], col_nms)
      )
    }
    NewAssays[[j]] <- do.call('rbind', lst_mat)
  }
  names(NewAssays) <- names(Assays)
  md <- c(list(
    TreeSummarizedExperiment = tse
  ), tse@metadata)
  md <- md[!duplicated(names(md))]
  SummarizedExperiment::SummarizedExperiment(
    assays = NewAssays,
    rowData = methods::as(keep, "DataFrame"),
    colData = SummarizedExperiment::colData(tse),
    metadata = md
  )
  
  
}


# tse_test <- tse_by_climb(dtse, train_res_climbed)
# tse_test

# plot_pvals(dplyr::bind_rows(method2$.data)) +
#   facet_wrap(~delta, scales = "free_y")

### VST ----
# .delta <-  seq(.5, 6, by = .5)
# tses <- lapply(.delta, shift_mu, tse = tse)
# .delta <- rep(.delta, each = 50)
# tses <- rep(tses, each = 50)
# 
# #sim counts
# tses <- purrr::map(tses, sim_counts, .progress = "sim-counts")
# tses <- purrr::map(tses,
#                    \(x) {
#                      suppressMessages({
#                        mat <- varianceStabilizingTransformation(assay(x))
#                      })
#                      hc <- hclust(dist(mat))
#                      rowTree(x) <- ape::as.phylo(hc)
#                      x
#                      
#                    }, .progress = "vst")
# ntses <- purrr::map(tses,
#                     tse_by_nodes, .progress = "expand-nodes")
# ddss <- purrr::map2(ntses, .delta,
#                    \(x, d) {
#                      suppressMessages({
#                        dds <- deseq(x, ~ groups, size_factors = 1L)[[2L]]
#                      })
#                      dds[['delta']] <- d
#                      dds
#                    }, .progress = "dds")
# # ddss_clmb <- purrr::map(ddss, climb, method = "ord", .progress = "climbing")
# # saveRDS(ddss_clmb, "rsd_cache/vst-dds-climb-res.rds")
# ddss_clmb <- readRDS("rsd_cache/vst-dds-climb-res.rds")
# plot_pvals(dplyr::bind_rows(ddss_clmb), sig_nodes_only = F) +
#   facet_wrap(~delta)

#reset

# ntse <- tse_by_nodes(tse)
# deseq(ntse, ~ groups, size_factors = 1L, ...)
# 
# splits <- purrr::map(tses, split_counts)
# tests <- purrr::map(splits, `[[`, 2L)
# train <- purrr::map(splits, `[[`, 1L)
# trained <- purrr::map2(train, 
#                 delta,
#                 function(x, d, ...) {
#                   ntse <- tse_by_nodes(x)
#                   out <- deseq(ntse, ~ groups, size_factors = 1L, ...)[[2L]]
#                   out[['delta']] <- delta
#                   out
#                 }, ..., .progress = "train-data")
d
## Error Control Count Split ----
### H_0 null hypothesis error control experiment ----

# tse_null <- do_sim(tse = tse, delta = 0, n_sims = 500)
# saveRDS(tse_null, "rds_cache/sim-null.rds")
tse_null <- readRDS("rds_cache/sim-null.rds")
tse_null_res <- do_climb(tse_null, method = "order", note = "null", affected_rows = integer())

plot_pvals(dplyr::bind_rows(tse_null_res$.data))

plot_pvals2 <- function(data) {
  dplyr::bind_rows(data) |> 
    subset(keep) |>
    dplyr::mutate(internalNode = factor(ifelse(n_children>1, "Internal Node", "Leaf"),c("Leaf", "Internal Node"))) |>
    ggplot(aes(x = pvalue)) +
    geom_histogram(aes(y = after_stat(count),
                       fill = internalNode), bins = 15, position = "stack") +
    # geom_density(aes(y = after_stat(density)), color = "red") +
    # facet_grid(cols = vars(internalNode)) +
    theme_classic() +
    theme(legend.position = "top") +
    scale_fill_manual(values = c("steelblue", "darkred")) +
    labs(fill = "")
}



tse_null_split <- do_splitCounts(tse = tse, delta = 0, n_sims = 500)
plot_pvals2(tse_null_split$trained_res$.data)
test_res <- lapply(tse_null_split$tested, function(x) x$result)
failed_tests <- which(vapply(test_res, \(x) identical(x,"FAILED"), T ))
test_null_split <- test_res[-failed_tests] |>
  lapply(\(x) {
    x$internalNode <- factor(ifelse(as.integer(sub('node_','',rownames(x)))>44L, "Internal Node", "Leaf"),c("Leaf", "Internal Node"))
    x$keep <- TRUE
    tibble::as_tibble(x)
  }) |>
  dplyr::bind_rows() 
test_null_split |>
  ggplot(aes(x = pvalue)) +
  geom_histogram(aes(y = after_stat(count/n*1000),
                     fill = internalNode), bins = 15, position = "stack") +
  geom_density(aes(y = after_stat(density)), color = "red") +
  # facet_grid(cols = vars(internalNode)) +
  theme_classic() +
  theme(legend.position = "top") +
  scale_fill_manual(values = c("steelblue", "darkred")) +
  labs(fill = "")

### simulate split tree with deltas

# tse_split <- do_splitCounts(tse = tse, delta = seq(.5, 6, by = .5), n_sims = 50)
# saveRDS(tse_split, "rds_cache/splitCounts-sim-deltas.rds")
tse_split <- readRDS("rds_cache/splitCounts-sim-deltas.rds")
plot_pvals2(tse_split$trained_res$.data)
test_split_res <- lapply(tse_split$tested, function(x) x$result)
failed_tests_sims <- which(vapply(test_split_res, \(x) identical(x,"FAILED"), T ))
test_split <- test_split_res |>
  purrr::imap(\(x, y ) {
    if (y %in% failed_tests_sims) return(NULL)
    x$internalNode <- factor(ifelse(as.integer(sub('node_','',rownames(x)))>44L, "Internal Node", "Leaf"),c("Leaf", "Internal Node"))
    x$keep <- TRUE
    x$sim_n <- y
    tibble::as_tibble(x, rownames = "node")
  }) |>
  dplyr::bind_rows() 
test_split |>
  ggplot(aes(x = pvalue)) +
  geom_histogram(aes(y = after_stat(count/n*1000),
                     fill = internalNode), bins = 15, position = "stack") +
  geom_density(aes(y = after_stat(density)), color = "red") +
  # facet_grid(cols = vars(internalNode)) +
  theme_classic() +
  theme(legend.position = "top") +
  scale_fill_manual(values = c("steelblue", "darkred")) +
  labs(fill = "")

#cheating here
metadata <- method1$.data[[1]][,c("nodes", "descendants", "n_children")]
. <- left_join(test_split, metadata, by = c("node"="nodes")) |>
  mutate(is_affected = vapply(descendants, \(x, y) all(x %in% y), logical(1), y = row_mod),
         is_affected2 = vapply(descendants, \(x, y) any(x %in% y), logical(1), y = row_mod))


#### Investiagting loglikelihood ----

# nbin_ll <- function(x, size, prob) {
#   n <- length(x)
#   sum(log(gamma(x + size)/(gamma(size) * factorial(x)))) +
#     (n * size * log(prob)) +
#     (sum(x) * log(1 - prob))
# }
# 
# df <- tidyr::expand_grid(size = 60:140, prob = seq(0, 1, .01)[-c(1, 101)])
# 
# vec <- as.vector(assay(dtse)[1:44,1:5])
# vec2 <- as.vector(assay(dtse)[row_mod, 6:10])
# 
# df$ll <- purrr::map2_dbl(df$size, df$prob, nbin_ll, x = vec)
# 
# ggplot(df, aes(x = size, y = prob, color = ll)) +
#   geom_point() +
#   geom_point(
#     data = ~dplyr::slice(.x,which.max(ll),.by = size),
#     color = "red"
#   ) +
#   scale_color_gradientn(values = scales::rescale(c(0,-5,-10,-25,-50, -150, -300, -500, -1000, -25000)),
#                         colors = scales::viridis_pal()(9))
