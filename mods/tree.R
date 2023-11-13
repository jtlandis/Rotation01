
box::use(methods[as],
         S4Vectors[DataFrame],
         SummarizedExperiment[...],
         stats[as.hclust])

suppressMessages(box::use(
  TreeSummarizedExperiment[...]))


#' @export
dendro_data <- function(model) {
  merge <- model$merge
  height <- model$height
  ord <- model$order
  x_start <- x_end <- y_start <- y_end <- double(3 * length(height))
  seq_ <- -2:0
  for (i in seq_len(length(height))) {
    
    h <- height[i]
    ii <- (i * 3) + seq_
    
    if ((m1 <- merge[i,1])<0) {
      x1 <- which(ord==abs(m1))
      y1 <- 0
    } else {
      sel <- m1*3
      x1 <- (x_start[sel] + x_end[sel])/2
      y1 <- y_start[sel]
    }
    
    if ((m2 <- merge[i, 2])<0) {
      x2 <- which(ord==abs(m2))
      y2 <- 0
    } else {
      sel <- m2*3 
      x2 <- (x_start[sel] + x_end[sel])/2
      y2 <- y_start[sel]
    }
    
    x_start[ii] <- c(x1, x2, x1)
    x_end[ii] <-  c(x1, x2, x2)
    y_start[ii] <- c(y1, y2, h)
    y_end[ii] <- h
    
    
  }
  
  out <- data.frame(x = x_start, y = y_start, xend = x_end, yend = y_end)
  attr(out, 'label') <- model$labels[ord]
  out
  
}

#' @name Find Common Parent
#' @description
#' Given two indexes `leaf1` and `leaf2`, return the index
#' of a `hclust` object's `$merge` element such that it is
#' the most recent common ancestor
#' 
#' @param .data an object of class `hclust`
#' @param leaf1 an scalar integer
#' @param leaf2 an scalar integer
#' @return a scalar integer
#' @export
find_common_parent_node_merge_index <- function(.data, leaf1, leaf2) {
  # browser()
  m <- .data$merge
  n_size <- nrow(m)
  
  target_leaf <- leaf1
  
  merge_index <- which(m==-target_leaf)
  if (merge_index>n_size) {
    side_check <- 1L
    merge_index <- merge_index - n_size
  } else {
    side_check <- 2L
  }
  
  #check other side
  node <- m[merge_index, side_check]
  if (node<0) {
    # this is a leaf
    if (abs(node)==leaf2)
      return(merge_index)
  } else {
    # positive node is a merge point
    if (leaf_in_index(.data, index = node, leaf = leaf2)) 
      return(merge_index)
  }
  
  # otherwise continue upwards
  find_common_parent_node_merge_index(.data, leaf1 = -merge_index, leaf2)
}

#' @name leaf_in_index
#' @description
#' Test if a leaf is a descendant of a merge index
#' @param .data an object of class `hclust`
#' @param index merge index (scalar integer)
#' @param leaf scalar integer
#' @export
leaf_in_index <- function(.data, index, leaf) {
  # browser()
  m <- .data$merge
  leafs <- m[index,]
  if (-leaf %in% leafs) {
    return(TRUE)
  } 
  pos_leafs <- leafs[leafs>0]
  for (l in pos_leafs) {
    new_res <- leaf_in_index(.data, l, leaf)
    if(new_res) return(TRUE)
  }
  
  FALSE
  
}

#' Helper function for finding descendants of 
#' nodes.
#' @export
descendants_of_nodes <- function(nodes, clust) {
  m <- clust$merge
  n_internal_nodes <- nrow(m)
  n_leafs <- n_internal_nodes + 1L
  n_nodes <- n_leafs + n_internal_nodes
  
  out <- vector('list', length = length(nodes))
  internal_nodes <- nodes > n_leafs
  out[!internal_nodes] <- clust$order[nodes[!internal_nodes]]
  out[internal_nodes] <- lapply(nodes[internal_nodes] - n_leafs, descendants_of_merge_index, .data = clust)
  out
}

#' @name descendants_of_merge_index
#' @description
#' Provide a list of merge indexs and return all leafs that
#' are children of the merges.
#' @param .data hclust object
#' @param index vector of indexes
#' @return vector of leaf indexes (in relation to clustered data)
#' @export
descendants_of_merge_index <- function(.data, index) {
  m <- .data$merge
  out <- vector('list', nrow(m))
  i <- 1L
  while(length(index)>0) {
    v <- as.vector(m[index,])
    pos <- v>0
    out[[i]] <- abs(v[!pos])
    index <- v[pos]
    i <- i + 1L
  }
  out <- out[1:(i-1L)]
  out <- unlist(out)
  out[order(out)]
}

#' @name gen_sub_grps
#' @description
#' Walk the tree and generate meta data regarding
#' the possible sub nodes and their descendants. 
#' Given a clustering of N objects, the first
#' N objects are the leafs and the next (N+1):(2N - 1)
#' are the internal nodes, ordered by the heights.
#' 
#' @param .data hclust object
#' @export
gen_node_data <- function(.data) {
  m <- .data$merge
  o <- .data$order
  l <- .data$labels
  # names(o) <- l
  h <- .data$height
  n_nodes <- nrow(m)
  out <- vector("list", nrow(m))
  
  for (i in seq_along(out)) {
    o_ind <- m[i,]
    
    child1 <- if ((o1 <- o_ind[1]) < 0) {
      -o1
    } else {
      out[[o1]]
    }
    
    child2 <- if ((o2 <- o_ind[2]) < 0) {
      -o2
    } else {
      out[[o2]]
    }
    descendants <- c(child1, child2)
    out[[i]] <- descendants[order(descendants)]
  }
  
  dplyr::bind_rows(
    tibble::tibble(node = 1:length(o),
                   height = 0,
                   descendants = lapply(o, identity),
                   n_children = 1L),
    tibble::tibble(node = (n_nodes + 2L):(2L*n_nodes + 1L),
                   height = h,
                   descendants = out,
                   n_children = vapply(out, length, integer(1)))
  )
  
  
}


#' @name tse_by_nodes
#' @description
#' Adapting a TreeSummarizedExperiment object to summarize counts across
#' features (genes) such that each internal node is a new feature. This should 
#' increase the size of the features available. The original data is stored
#' under the `@metadata` slot as "TreeSummarizedExperiment"
#' @param tse A TreeSummarizedExperiment object
#' @param summarise_fun a function for summarizing a matrix (columnwise), 
#' For a N by M tse object, the function should return a vector of size M.
#' The default is `base::colSums`
#' @return A SummarizedExperiment object
#' @export
tse_by_nodes <- function(tse, summarise_fun = base::colSums) {
  hclust_data <- TreeSummarizedExperiment::rowTree(tse) |> as.hclust()
  node_data <- gen_node_data(hclust_data)
  
  n_cols <- ncol(tse)
  col_names <- colnames(tse)
  
  out <- vector("list", nrow(node_data)) 
  col_fn <- summarise_fun
  sprint_fmt <- paste0("node_%0", floor(log10(nrow(node_data)))+1, "i")
  summarise_assay <- function(x, ind, rname) {
    matrix(
      data = col_fn(x[ind, ,drop = F]),
      nrow = 1, ncol = n_cols, dimnames = list(rname, col_names)
    )
  }
  
  Assays <- assays(tse)
  
  for (i in seq_len(nrow(node_data))) {
    
    node_name <- sprintf(sprint_fmt, i)
    ind <- node_data$descendants[[i]]
    # sub_tse <- tse[ind,]
    out[[i]] <- lapply(Assays, summarise_assay, ind = ind, rname = node_name)
    
  }
  new_Assays <- vector("list", length(Assays))
  
  for (i in seq_along(Assays)) {
    new_Assays[[i]] <- do.call(rbind, lapply(out, `[[`, i))
  }
  names(new_Assays) <- names(Assays)
  
  SummarizedExperiment(
    assays = new_Assays,
    rowData = as(node_data, "DataFrame"),
    colData = colData(tse), 
    metadata = list(
      row_hclust = hclust_data,
      TreeSummarizedExperiment = tse
    )
  )
  
}