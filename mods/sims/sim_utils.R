
box::use(
  ../tree[...],
  ../stamps[...],
  stats[as.hclust, rnbinom, runif],
  dplyr[mutate, select, summarise, filter, reframe, left_join,
        arrange, bind_rows, n, slice],
  tidyr[unnest],
  tibble[as_tibble, tibble],
  SummarizedExperiment[assay, `assay<-`],
  methods[is],
  # S7[...],
  ggplot2[...]
)

#' @export
box::use(
  ../rnaseq[...],
  TreeSummarizedExperiment[TreeSummarizedExperiment, rowTree, `rowTree<-`],
  ../tree[descendants_of_nodes]
)
## Sim Data ----



### Functions ----

#' @description
#' potential utility function to make random groups
mk_group_members <- function(nSamps, ngroups) {
  n_per_group <- rep_len(nSamps %/% ngroups, ngroups)
  r <- nSamps %% ngroups
  if (r>0) {
    n_per_group[1:r] <- n_per_group[1:r] + 1L
  }
  n_per_group
}

#' @description
#' Utility function to assign group memebers
#' @export
assign_group_members <- function(nSamps, ngroups, random = F) {
  s <- 1:nSamps
  out <- vector('list', ngroups)
  grp_members <- mk_group_members(nSamps, ngroups)
  if (random) {
    for (i in seq_len(ngroups)) {
      s_ <- sample(s, grp_members[i], replace = F)
      s <- setdiff(s, s_)
      out[[i]] <- s_
    }
    out <- lapply(out, function(x) x[order(x)])
  } else {
    start <- 1L
    csum <- cumsum(grp_members)
    for (i in seq_len(ngroups)) {
      end <- csum[i]
      out[[i]] <- start:end
      start <- end + 1L
    }
  }
  ans <- integer(nSamps)
  for (i in seq_len(ngroups)) {
    ans[out[[i]]] <- i
  }
  ans
  
}

#' @description
#' Meant to be called on a results object and the Summarized experiment object
#' That DESeq2 was run on
#' @export
collect_data <- function(res, se, clust) {
  res <- as_tibble(res, rownames = "nodes")
  dendro_d <- dendro_data(model = clust)
  n_orig <- length(clust$order)
  n <- (2L*n_orig) - 1L
  sprintfmt <- paste0("node_%0", floor(log10(n))+1L, "i")
  
  rdata <- rowData(se) |> as_tibble(rownames = "nodes")
  
  dendro_d$nodes <- rep(sprintf(sprintfmt, as.integer((n_orig + 1L):n)), each = 3)
  
  .data <- dendro_d |>
    as_tibble() |>
    slice(n(), .by = nodes)  |>
    reframe(
      node_x = (x + xend)/2,
      node_y = y,
      nodes = nodes
    ) |>
    left_join(
      y = tidyr::nest(dendro_d, segments = x:yend, .by = nodes),
      by = "nodes"
    ) |>
    bind_rows(
      tibble::tibble(
        node_x = 1:n_orig,
        node_y = 0,
        nodes = sprintf(sprintfmt, 1:n_orig)
      )
    ) |>
    arrange(nodes) |>
    left_join(res, by = "nodes")
  
  .data <- .data |>
    left_join(y = rdata, by = "nodes")
  
  attr(.data, "hclust") <- clust
  class(.data) <- c("dds_node_res", class(.data))
  .data
}


#' @description
#' Plot a dendrogram with information regarding DESeq2
#' after creating observations for each node
#' @export
plot_results <- function(collected_results, alpha = NULL) {
  
  p <- ggplot(collected_results, aes(x = node_x, y = node_y)) +
    geom_segment(
      aes(x = x, y = y, xend = xend, yend = yend, color = pvalue),
      data = \(.x) tidyr::unnest(collected_results, segments)
    ) 
  
  p <- p +
    scale_color_gradientn(
      colors = rev(scales::viridis_pal(option="D")(5)),
      trans = "log10",
      limits = \(.x) c(10^(floor(log10(.x[1]))), 1),
      breaks = \(.x) {.x <- log10(.x)
      round(10^(seq(.x[1], .x[2], length.out = 5)), abs(.x[1]))}
    ) + theme_classic() +
    theme(axis.line.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "top",
          legend.title = element_text(face = "bold", vjust = .8),
          legend.key.width = unit(x = .8, units = "in")) +
    labs(y = "Height")
  
  if (is.null(alpha)) {
    p <- p + 
      geom_point(
        aes(color = pvalue)
      ) +
      geom_text(
        aes(label = sub("node_", "", nodes)),
        color = "grey40",
        angle = 90,
        hjust = 1.2
      )
  } else {
    
    p <- p +
      geom_point(
        data = \(.x) filter(.x, pvalue <= .env$alpha),
        color = "red"
      ) +
      geom_point(
        data = \(.x) filter(.x, pvalue > .env$alpha | is.na(pvalue)),
        color = "grey70"
      ) +
      geom_text(
        aes(label = sub("node_", "", nodes)),
        data = \(.x) filter(.x, pvalue <= .env$alpha),
        color = "red",
        angle = 90,
        hjust = 1.2
      ) +
      geom_text(
        aes(label = sub("node_", "", nodes)),
        data = \(.x) filter(.x, pvalue > .env$alpha | is.na(pvalue)),
        color = "grey70",
        angle = 90,
        hjust = 1.2
      )
    
  }
  
  if ("keep" %in% names(collected_results)) {
    p <- p +
      geom_point(
        data = dplyr::filter(collected_results, keep),
        color = "red", size = .5
      )
  }
  
  
  p
}

#' @param .data a result object
#' @param ... unused
climb_pvalue_arrange <- function(.data, merge_alpha = 1, ...) {
  # res <- out$.data[[120]]
  # 
  n_leafs <- sum(.data$n_children==1)
  
  # leafs <- filter(res, n_children==1) |>
  #   tidyr::unnest(descendants)
  desc <- numeric()
  ord <- arrange(.data, pvalue)
  ord$keep <- FALSE
  i <- 0L
  theDesc <- ord$descendants
  pvals <- ord$pvalue
  nchild <- ord$n_children
  while(length(desc)<n_leafs) {
    i <- i + 1L
    desc2 <- theDesc[[i]]
    simi <- any(desc2 %in% desc)
    if (!simi && ( nchild[i]==1L || pvals[i] < merge_alpha)) {
      desc <- c(desc, desc2)
      ord$keep[i] <- TRUE
    }
  }
  arrange(ord, node)
}


#' @export
climb <- function(.data, method = c("strict", "ordered_pval"), ...) {
  
  method <- match.arg(method, c("strict", "ordered_pval"))
  f <- switch(method,
              strict = climb_strict,
              ordered_pval = climb_pvalue_arrange)
  
  f(.data, ...)
}

#' @description
#' climb an expanded results data with a tree.
#' @returns A data.frame with new logical column `$keep` 
climb_strict <- function(.data, ...) {
  .clust <- attr(.data, "hclust")
  m <- .clust$merge
  o <- .clust$order
  n <- nrow(m) + 1
  .data$keep <- NA
  pvals <- .data$pvalue
  is_na <- is.na(pvals)
  pvals[is_na] <- 1
  for (i in seq_len(nrow(m))) {
    
    merge_i <- m[i,]
    m1 <- merge_i[1]
    m2 <- merge_i[2]
    
    node1 <- if (m1 < 0) 
      which(o==-m1)
    else 
      m1 + n
    
    node2 <- if (m2 < 0) 
      which(o==-m2)
    else
      m2 + n
    
    parent_node <- i + n
    
    nodes <- c(node1, node2)
    #Before we check p-values, check if we need to keep
    # climbing. Specifically if the parent is FALSE
    # we do not need to compare anything.
    dont_use <- vapply(.data$keep[nodes], isFALSE, logical(1))
    if (any(dont_use)) {
      .data$keep[parent_node] <- FALSE
      # some data may be newly merged and havent been tested
      # we don't need to compare p-values for these if we know
      # this merge index wont be used
      .data$keep[nodes[!dont_use]] <- TRUE
      next
    }
    
    node_p <- pvals[nodes]
    if (all(node_p > pvals[parent_node])) {
      .data$keep[nodes] <- FALSE
      .data$keep[parent_node] <- TRUE
    } else {
      .data$keep[nodes] <- TRUE
      .data$keep[parent_node] <- FALSE
    }
    
    
  }
  
  # .data$pvalue[is_na] <- NA_real_
  .data
}

#' @description
#' Summarises information returned from `climb()` by
#' providing insight in which nodes were merged. 
#' Assuming we know which rows were actually affected we
#' have sensitivity and persision data returned
#' 
#' True Positive  <--> merged_affected
#' False Negative <--> unmerged_affected
#' False Positive <--> merged_unaffected
#' True Negative  <--> unmerged_unaffected
#' 
#' @export
summarise_climb <- function(.data, affected_rows, ...) {
  
  .data$n_affected <-
    vapply(.data$descendants,
           function(x, y) sum(x %in% y), y = affected_rows,
           FUN.VALUE = integer(1))
  kept <- filter(.data, keep) |>
    mutate(are_merged = n_children > 1)
  
  summarise(
    kept,
    ...,
    merge_from = sum(n_children[are_merged]),
    merge_to = sum(are_merged),
    unmerged = sum(!are_merged),
    TP = sum(n_affected[are_merged]),
    FN = sum(n_affected[!are_merged]),
    FP = merge_from - TP,
    TN = unmerged - FN,
    total_affected = sum(TP, FN),
    sens = TP/(TP + FN),
    spec = TN/(TN + FP),
    FDR = FP/(TP + FP),
    FPR = FP/(TN + FP)
  )
  
}


#' @name sim_mu
#' Modify assay matrix of the 
#' `mu` parameters for `rnbinom()`
#' @param tse A object inheriting from 
#' SummarizedExperiment class, such as a
#' TreeSummarizedExperiment.
#' @param affected_rows rows in which to add delta
#' @param affected_cols columns in which to add delta
#' @param delta The actual value to shift by
#' @export
sim_mu <- function(se,
                   affected_rows = 1:nrow(se),
                   affected_cols = 1:ncol(se),
                   delta = 0) {
  
  mat_mu <- tryCatch(assay(se, "mu"),
                     error = function(cnd) {
                       matrix(0, nrow = nrow(se),
                              ncol = ncol(se), dimnames = dimnames(se))
                       })
  mat_mu[affected_rows, affected_cols] <- mat_mu[affected_rows, affected_cols] + delta
  assay(se, 'mu') <- mat_mu
  se
  
  
}

#' @describeIn sim_mu add variation randomly in an object
#' @param iterations number of times to add variation
#' @param delta_range vecter of numbers, whose min and max
#' values are uniformly sampled from for each
#' @export
sim_mu_random <- function(se,
                          iterations = 5,
                          delta_range = c(-1, 1)) {
  nrows <- nrow(se)
  r_seq <- 1:nrows
  ncols <- ncol(se)
  c_seq <- 1:ncols
  unif_min <- min(delta_range)
  unif_max <- max(delta_range)
  for (i in iterations) {
    nrow_mod <- sample(r_seq, 1L)
    ncol_mod <- sample(c_seq, 1L)
    rows_mod <- sample(r_seq, nrow_mod, replace = F)
    cols_mod <- sample(c_seq, ncol_mod, replace = F)
    delta <- runif(1, unif_min, unif_max)
    se <- sim_mu(se, affected_rows = rows_mod, affected_cols = cols_mod, delta = delta)
  }
  se
}

#' @describeIn sim_mu adds delta to decendants of a specific
#' node index. For a TreeSummarizedExperiment object of N rows, leaf nodes
#' are indexed 1:N by their hclust ordering (ie se[1,] is not necessarily node 1.). 
#' nodes (N + 1):(2*N -1) are indexed by their heights.
#' @param affected_nodes
#' @export
sim_mu_node <- function(tse,
                        affected_nodes = 2*nrow(tse) - 1,
                        affected_cols = 1:ncol(tse),
                        delta = 0) {
  
  nrows <- nrow(tse)
  nnodes <- (2*nrows)-1L
  if (any(affected_nodes > nnodes)) {
    stop("Specifying nodes that do not exist", call. = F)
  }
  
  clust <- rowTree(tse) |>
    as.hclust()
  
  affected_rows <- descendants_of_nodes(affected_nodes, clust)
  affected_rows <- unique(unlist(affected_rows))
  sim_mu(tse, affected_rows, affected_cols, delta)
  
  
  
  
}

#' @describeIn sim_mu simulates counts on a SummarizedExperiment object using `rnbinom(..., size = 100)`
#' @param size size parameter in `rnbinom()`
#' @export
sim_counts <- function(se, size = 100) {
  
  
  mat_mu <- tryCatch(assay(se, 'mu'),
                     error = function(cnd) {
                       stop("Attempted to call `sim_counts()` when no `mu` assay exists. Use `sim_mu()` first", call. = F)
                     })
  assay(se, 'counts')[] <- rnbinom(prod(dim(se)), size = size, mu = mat_mu)
  se
  
}


#' @export
deseq <- function(node_se, design, ..., size_factors = NULL) {
  
  stopifnot("node_se is likely not in correct format. please use `tse_by_node()`"=is(node_se,"SummarizedExperiment")
            &&('TreeSummarizedExperiment' %in% names(node_se@metadata)))
  suppressMessages({
    dds <- DESeqDataSet(node_se, design = design)
    if (!is.null(size_factors))
      sizeFactors(dds) <- size_factors
    dds <- DESeq(dds, quiet = T, ...)
    res_names <- resultsNames(dds)
    out <- lapply(res_names, function(x) results(object = dds, name = x))
  })
  out <- lapply(out, collect_data, se = node_se, clust = node_se@metadata$row_hclust)
  names(out) <- res_names
  out
}



#' @description
#' A short description...
#' @export
simulate_climb <- function(tse, affected_rows, affected_cols, mu, mu_delta) {
  
  time_start()
  clust <- as.hclust(rowTree(tse))
  tse <- sim_counts(tse)
  tse_nodes <- tse_by_nodes(tse)
  
  rowData(tse_nodes)$desc_affected <-
    vapply(rowData(tse_nodes)$descendants,
           function(x, y) sum(x %in% y), y = affected_rows,
           FUN.VALUE = integer(1))
  ##Assays sim ----
  ### DESeq2 ----
  time_stamp()
  colData(tse_nodes)$groups <- factor(tse_nodes$groups)
  suppressMessages({
    dds_1 <- DESeqDataSet(tse_nodes, design = ~ groups)
    sizeFactors(dds_1) <- rep(1L, ncol(dds_1))
    dds_1 <- DESeq(dds_1, quiet = TRUE)
    res <- results(dds_1, name = "groups_2_vs_1")
  })
  
  
  time_stamp()
  res_data <- collect_data(res, tse_nodes, clust) |>
    climb(clust)
  
  out <- summarise_climb(res_data)
  
  out$.params <- list(tibble(rows = list(affected_rows),
                      cols = list(affected_cols),
                      mu = mu,
                      mu_delta = mu_delta))
  #out$p <- list(plot_results(res_data))
  out$.data <- list(res_data)
  time_end()
  out
  
}

#' @export
plot_pvals <- function(climb_res, sig_nodes_only = TRUE) {
  
  plot_data <- climb_res
  if (sig_nodes_only)
    plot_data <- subset(plot_data, keep)
  
  ggplot(plot_data, aes(x = pvalue)) +
    geom_histogram(aes(y = after_stat(density)),
                   fill = "steelblue", color = "steelblue",bins = 15) +
    geom_density(color = "red") +
    theme_classic()
  
  
}
