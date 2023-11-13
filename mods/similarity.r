

box::use(
  S7[...],
  mods/tree[dendro_data],
  mods/gtf[exon_df, gene_df, transcript_df, extract_name, extract_attribute],
  mods/util[vec_rep],
  ggplot2[...],
  SummarizedExperiment[rowData],
  stats[dist, hclust]
)

box::use(ggside[geom_ysidesegment]) |>
suppressMessages()

change <- S7::new_external_generic("S7", "convert", dispatch_args = c("from", "to"))

## classes ----

GR <- methods::getClass("GenomicRanges")
GRL <- methods::getClass("GenomicRangesList")
SE <- methods::getClass("SummarizedExperiment")
RSE <- methods::getClass("RangedSummarizedExperiment")


## Generics ----

#' generates hclust from data
#' @param data an object
#' @param similarity_method one of 'sorensen', 'basepair',
#' or 'counts'.
#' @param dist_method method to provide to `stats::dist`
#' @param clust_method method to provide to `stats:hclust`
#' @export
cluster_data <- new_generic(
  "cluster_data",
  "data",
  fun = function(data, 
                 similarity_method = c("sorensen", "basepair", "counts"),
                 dist_method = "euclidean",
                 clust_method = "complete", name = NULL) {
    S7_dispatch()
})
  
#' @export
calc_Sorensen_coef <- new_generic(
  "calc_Sorensen_coef",
  "data"
)

#' @export
plot_cluster <- new_generic(
  "cluster_data",
  "data",
  fun = function(data, 
                 similarity_method = c("sorensen", "basepair", "counts"),
                 dist_method = "euclidean",
                 clust_method = "complete", name = NULL) {
    S7_dispatch()
  }
)


## Methods ----

#### cluster_data ----

# browser()
method(cluster_data, transcript_df) <-  function(data, 
                                   similarity_method = c("sorensen", "basepair", "counts"),
                                   dist_method = "euclidean",
                                   clust_method = "complete", name = NULL) {
         similarity_method <- match.arg(similarity_method, choices = c("sorensen", "basepair", "counts"))
         if (similarity_method=="basepair") {
           simi_mat <- bp_similarity_(data$start, data$end)
         } else if(similarity_method == "sorensen"){
           stopifnot("`exons` must be in the column names of data for sorensen coef."="exons" %in% names(data))
           simi_mat <- calc_Sorensen_coef(data)
         } else {
           stop("`counts` similarity method is not available")
         }
         if (!is.null(name)) {
           labels <- extract_name(data$attribute, type = name)
           rownames(simi_mat) <- colnames(simi_mat) <- labels
         }
         dist_mat <- dist(simi_mat, method = dist_method)
         clust <- hclust(dist_mat, method = clust_method)
         clust
       }

GR_cluster_data <- function(data,
                            similarity_method = c("sorensen", "basepair", "counts"),
                            dist_method = "euclidean",
                            clust_method = "complete", name = NULL) {
  # check if this object is compatible
  if (!('exons' %in% names(data@elementMetadata) && methods::isClass("GenomicRangesList", data$exons))) {
    stop("data should be generated with `read_gtf2GR`")
  }
  similarity_method <- match.arg(similarity_method, choices = c("sorensen", "basepair", "counts"))
  simi_fun <- switch(similarity_method,
                     sorensen = gr_sorensen_simi,
                     basepair = function(chrom_i) {
                       chrom_exons <- chrom_i$exons
                       chrom_exons <- unlist(chrom_exons)
                       gr_bp_similarity(chrom_exons)
                     },
                     counts = stop("method not implemented"))
  by_chrom <- GenomicRanges::split(data, data@seqnames)
  by_chrom <- by_chrom[vapply(by_chrom, length, integer(1))>0]
  n_chrom <- length(by_chrom)
  n_trans <- length(data)
  nms <- names(data)
  mat <- Matrix::sparseMatrix(dims = c(n_trans, n_trans),i = {}, j = {}, dimnames = list(nms, nms), x = numeric())
  sizes <- vapply(by_chrom, length, integer(1))
  index_end <- cumsum(sizes)
  index_start <- c(1L, index_end[-n_chrom] + 1L)
  for (i in cli::cli_progress_along(by_chrom)) {
    i_seq <- index_start[i]:index_end[i]
    chrom_i <- by_chrom[[i]]
    mat[i_seq, i_seq] <- simi_fun(chrom_i)
  }
  dist_mat <- dist(mat, method = dist_method)
  clust <- hclust(dist_mat, method = clust_method)
  clust
}


change  <- function(from, to) {
  GenomicR <- GenomicRanges::ranges(from)
  GenomicR@elementMetadata <- rowData(from)
  GenomicR
}

method(cluster_data, GR) <- GR_cluster_data



method(cluster_data, RSE) <-
  function(data,
           similarity_method = c("sorensen", "basepair", "counts"),
           dist_method = "euclidean",
           clust_method = "complete", name = NULL) {
         GR_cluster_data(data@rowRanges, similarity_method, dist_method, clust_method, name)
  }

#### Calc_Soren ----

method(calc_Sorensen_coef, transcript_df) <- function(data) {
  trans_data <- data
  trans_data$trans_name <- extract_name(trans_data$attribute, "trans")
  n <- nrow(trans_data)
  exons <- dplyr::select(trans_data, trans_name, exons) |>
    tidyr::unnest(exons)
  exons$id <- seq_len(nrow(exons))
  emat <- bp_similarity_(exons$start, exons$end)
  trans_mat <- matrix(0, nrow = n, ncol = n)
  seq_lst <- .Internal(split(exons$id, factor(exons$trans_name, levels = trans_data$trans_name)))
  for (i in 1:n) {
    i_seq <- seq_lst[[i]]
    vec <- numeric(n + 1L - i)
    for (j in i:n) {
      j_seq <- seq_lst[[j]]
      mat <- emat[i_seq, j_seq, drop = FALSE]
      vec[j - i + 1L] <- 2*sum(mat) / (sum(dim(mat)))
    }
    i_seq2 <- i:n
    trans_mat[i_seq2, i] <- vec
    trans_mat[i, i_seq2] <- vec
  }
  trans_mat
}

encode <- function(s) {
  n <- length(s)
  seq_i <- 1:n
  out <- integer(sum(seq_i))
  k <- 0L
  for (i in seq_i) {
    k <- k + 1L
    k_ <- k + n - i
    out[k:k_] <- i:n
    k <- k_
  }
  out
}

gr_sorensen_simi <- function(data) {
  
  chrom_exons <- data$exons
  chrom_exons <- unlist(chrom_exons)
  exon_mat <- gr_bp_similarity(chrom_exons)
  n_exons <- length(chrom_exons)
  
  mat <- Matrix::sparseMatrix(1:n_exons, vec_rep(1:length(data), vapply(data$exons, length, integer(1))), x = 1L)
  trans_mat <- Matrix::t(mat) %*% exon_mat %*% mat
  diag_vec <- diag2(trans_mat)
  
  logic_mat <- trans_mat>0
  n_per_col <- Matrix::colSums(logic_mat)
  
  mat_dims <- trans_mat
  mat_dims[logic_mat] <- vec_rep(diag_vec, n_per_col)
  
  trans_mat[logic_mat] <- 2 * trans_mat[logic_mat] / ((mat_dims + Matrix::t(mat_dims))[logic_mat])
  # chrom_exons_seq_i <- seq_len(length(chrom_exons))
  # n_trans_data <- length(data)
  # # trans_data_mat <- matrix(0, nrow = n_trans_data, ncol = n_trans_data)
  # trans_names <- names(data)
  # # fctr <- vec_rep(1:n_trans_data, vapply(data$exons, length, integer(1)))
  # # levels(fctr) <- trans_names
  # # class(fctr) <- 'factor'
  # fctr <- factor(vec_rep(trans_names, vapply(data$exons, length, integer(1))), trans_names)
  # browser()
  # seq_lst <- .Internal(split(chrom_exons_seq_i, fctr))
  # i <- vec_rep(1:n_trans_data, n_trans_data:1)
  # j <- encode(1:n_trans_data)
  # out <- numeric(length = length(i))
  # is <- integer(length(i))
  # js <- integer(length(j))
  # k_ <- 0L
  # for (k in cli::cli_progress_along(i)) {
  #   i_ <- i[k]
  #   j_ <- j[k]
  #   i_seq <- seq_lst[[i_]]
  #   j_seq <- seq_lst[[j_]]
  #   mat_sub <- exon_mat[i_seq, j_seq, drop = FALSE]
  #   v <- 2*sum(mat_sub)/(sum(dim(mat_sub)))
  #   if (v > 0) {
  #     k_ <- k_ + 1L
  #     is[k_] <- i_
  #     js[k_] <- j_
  #     out[k_] <- v
  #   }
  #   
  # }
  # cli::cli_progress_done()
  # k_seq <- 1:k_
  # Matrix::sparseMatrix(
  #   i = is[k_seq], j = js[k_seq], x = out[k_seq],
  #   dims = c(n_trans_data, n_trans_data), symmetric = T
  # )
  # for (k in 1:n_trans_data) {
  #   k_seq <- seq_lst[[k]]
  #   vec <- numeric(n_trans_data + 1L - k)
  #   for (j in k:n_trans_data) {
  #     j_seq <- seq_lst[[j]]
  #     mat_subset <- exon_mat[k_seq, j_seq, drop = FALSE]
  #     vec[j - k + 1L] <- 2*sum(mat_subset) / (sum(dim(mat_subset)))
  #   }
  #   k_seq2 <- k:n_trans_data
  #   trans_data_mat[k_seq2, k] <- vec
  #   trans_data_mat[k, k_seq2] <- vec
  # }
  trans_mat
}

method(calc_Sorensen_coef, GR) <- gr_sorensen_simi

method(calc_Sorensen_coef, RSE) <- function(data) {
  gr_sorensen_simi(data@rowRanges)
}

#### plot_cluster ----


method(plot_cluster, transcript_df) <- function(data, similarity_method = c("sorensen", "basepair", "counts"),
                         dist_method = "euclidean", clust_method = "complete", name = NULL) {
  similarity_method <- match.arg(similarity_method, c("sorensen", "basepair", "counts"))
  data_type <- data_type(data)
  if (data_type$type == "exons") {
    stop("`plot_cluster` is meant to be used on a set of transcripts, not exons.")
  }
  clust <- cluster_data(data, similarity_method, dist_method, clust_method, name)
  data_clust <- dendro_data(clust)
  if (!is.null(name)) {
    if (anyDuplicated(clust$labels)) {
      clust$labels <- make.unique(clust$labels, sep = "_")
    }
    data$cluster_order <- factor(x = clust$labels, levels = clust$labels[clust$order])
  } else {
    data <- data[clust$order,]
    data$cluster_order <- 1:nrow(data)
    data <- data[order(clust$order),]
  }
  p <- ggplot(data, aes(x = start, y = cluster_order)) +
    geom_segment(aes(xend = end, yend = cluster_order)) +
    geom_ysidesegment(
      aes(x = y, y = x, xend = yend, yend = xend),
      data = data_clust,
      inherit.aes = F
    ) +
    theme_bw() +
    labs(x = "Genome Position",
         y = "Cluster Label",
         subtitle = sprintf("Similarity: %s | dist: %s | clust: %s", 
                            match.arg(similarity_method, c("sorensen", "basepair", "counts")),
                            dist_method, clust_method))
  if ("exons" %in% colnames(data)) {
    exon_data <- dplyr::select(data, cluster_order, exons) |>
      tidyr::unnest(exons)
    p <- p +
      geom_rect(
        aes(xmin = start, ymin = as.integer(cluster_order) - .25,
            xmax = end, ymax = as.integer(cluster_order) + .25),
        data = exon_data,
        fill = "blue", inherit.aes = F
      ) +
      labs(
        title = sprintf("transcripts from %s", data_type$parent)
      )
  }
  p
}


method(plot_cluster, GR) <-
  GR_plot_cluster <- function(data, similarity_method = c("sorensen", "basepair", "counts"),
                                     dist_method = "euclidean", clust_method = "complete", name = NULL){
  
  
  similarity_method <- match.arg(similarity_method, c("sorensen", "basepair", "counts"))
  if (!('exons' %in% names(data@elementMetadata) && methods::isClass("GenomicRangesList", data$exons))) {
    stop("data should be generated with `read_gtf2GR`")
  }
  
  
  clust <- cluster_data(data, similarity_method, dist_method, clust_method, name)
  data_clust <- dendro_data(clust)
  
  plot_data <- tibble::tibble(
    start = data@ranges@start,
    end = start + data@ranges@width - 1L,
    seqnames = as.character(data@seqnames)
  )
  if (!is.null(name)) {
    if (anyDuplicated(clust$labels)) {
      clust$labels <- make.unique(clust$labels, sep = "_")
    }
    plot_data$cluster_order <- factor(x = clust$labels, levels = clust$labels[clust$order])
  } else {
    plot_data <- plot_data[clust$order,]
    plot_data$cluster_order <- 1:nrow(plot_data)
    plot_data <- plot_data[order(clust$order),]
  }
  
  p <- ggplot(plot_data, aes(x = start, y = cluster_order)) +
    geom_segment(aes(xend = end, yend = cluster_order)) +
    geom_ysidesegment(
      aes(x = y, y = x, xend = yend, yend = xend),
      data = data_clust,
      inherit.aes = F
    ) +
    theme_bw() +
    labs(x = "Genome Position",
         y = "Cluster Label",
         subtitle = sprintf("Simularity: %s | dist: %s | clust: %s", 
                            "sorensen",
                            dist_method, clust_method))
  if ('exons' %in% names(data@elementMetadata)) {
    exon_ <- unlist(data$exons)
    exon_data <- tibble::tibble(
      start = exon_@ranges@start,
      end = start + exon_@ranges@width - 1L,
      seqnames = as.character(exon_@seqnames),
      cluster_order = vec_rep(plot_data$cluster_order, vapply(data$exons, length, 1L))
    )
    p <- p +
      geom_rect(
        aes(xmin = start, ymin = as.integer(cluster_order) - .25,
            xmax = end, ymax = as.integer(cluster_order) + .25),
        data = exon_data,
        fill = "blue", inherit.aes = F
      ) 
  }
  p
} 

method(plot_cluster, RSE) <- function(data, similarity_method = c("sorensen", "basepair", "counts"),
                              dist_method = "euclidean", clust_method = "complete", name = NULL) {
  similarity_method <- match.arg(similarity_method, c("sorensen", "basepair", "counts"))
    GR_plot_cluster(data@rowRanges, similarity_method, dist_method, clust_method, name)
}


#### Other ----

gr_intersect_bp <- function(start, stop, ref_start, ref_stop) {
  to_compare <- (ref_start <= stop & ref_stop >= start) #| (ref_start > start & ref_stop > stop)# only look at positions that will overlap
  #res <- numeric(length(start))
  indx <- which(to_compare)
  uni_start <- int_start <- start[indx]
  int_start[int_start < ref_start] <- ref_start
  uni_stop <- int_stop <- stop[to_compare]
  int_stop[int_stop > ref_stop] <- ref_stop
  uni_start[uni_start > ref_start] <- ref_start
  uni_stop[uni_stop < ref_stop] <- ref_stop
  intersection <- int_stop-int_start
  unionv <- uni_stop-uni_start
  res <- intersection/unionv
  if (any(is_Nan <- is.nan(res)) && any(is_0 <- intersection[is_Nan]==0)) {
    res[is_Nan & is_0] <- 1
  }
  list(indx = indx, res = res)
}

intersect_bp <- function(start, stop, ref_start, ref_stop) {
  # browser()
  to_compare <- (ref_start <= stop & ref_stop >= start) #| (ref_start > start & ref_stop > stop)# only look at positions that will overlap
  res <- numeric(length(start))
  uni_start <- int_start <- start[to_compare]
  int_start[int_start < ref_start] <- ref_start
  uni_stop <- int_stop <- stop[to_compare]
  int_stop[int_stop > ref_stop] <- ref_stop
  uni_start[uni_start > ref_start] <- ref_start
  uni_stop[uni_stop < ref_stop] <- ref_stop
  
  res[to_compare] <- (int_stop-int_start)/(uni_stop-uni_start)
  # if(any(res<0))
  #   browser()
  # res
  res
  
}


bp_similarity_ <- function(start, stop) {
  n <- length(start)
  mat <- matrix(nrow = n, ncol = n)
  for (i in seq_len(n)) {
    ref_start <- start[i]
    ref_stop <- stop[i]
    seq_i <- i:n
    vec <- intersect_bp(start[seq_i], stop[seq_i], ref_start = ref_start, ref_stop = ref_stop)
    mat[seq_i,i] <- vec
    mat[i,seq_i] <- vec
  }
  mat
}

gr_bp_similarity <- function(gr) {
  r <- gr@ranges
  start <- r@start
  stop <- start + r@width - 1L
  # number of transcripts
  n <- length(gr)
  # mat <- matrix(nrow = n, ncol = n)
  is <- vector('list', n)
  js <- vector('list', n)
  xs <- vector('list', n)
  i_shift <- 0L
  for (i in seq_len(n)) {
    ref_start <- start[i]
    ref_stop <- stop[i]
    seq_i <- i:n
    vec <- gr_intersect_bp(start[seq_i], stop[seq_i], ref_start = ref_start, ref_stop = ref_stop)
    indx <- vec$indx
    js[[i]] <- rep_len(i, length(indx))
    xs[[i]] <- vec$res
    is[[i]] <- indx + i_shift
    i_shift <- i_shift + 1L
  }
  Matrix::sparseMatrix(
    i = unlist(is),
    j = unlist(js),
    x = unlist(xs),
    dims = c(n, n),
    symmetric = T
  )
}






data_type <- function(data) {
  
  pattern <- c("gene_id", "transcript_id")
  this_attr <- data$attribute[1]
  if ("trans" %in% names(data)) {
    list(type = "genes",
         name = extract_attribute(this_attr, "gene_id"),
         parent = "none")
  } else if ("exons" %in% names(data)) {
    list(type = "transcripts",
         name = extract_attribute(this_attr, "transcript_name"),
         parent = extract_attribute(this_attr, "gene_id"))
  } else {
    list(type = "exons",
         name = extract_attribute(this_attr, "exon_id"),
         parent = extract_attribute(this_attr, "transcript_name"))
  } 
  
}




diag2 <- function(mat) {
  l <- nrow(mat)
  mat[(1:l) *l - l + 1:l]
}


