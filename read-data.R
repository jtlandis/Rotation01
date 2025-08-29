
box::use(mods/gtf,
         mods/similarity[...],
         mods/util[vec_rep],
         SummarizedExperiment[...],
         mods/gr
         )


sf_files <- list.files(box::file("data/PRJNA883409-Taquila-seq"), pattern = ".sf$", full.names = T)
sf_data <- readr::read_delim(sf_files, id = 'file', col_names = T, show_col_types = F)
meta <- readr::read_csv(box::file("data/PRJNA883409-Taquila-seq/meta_PRJA883409-filtered.csv")) 
meta <- meta[1:12,]

sf_data <- dplyr::mutate(
  sf_data,
  out = stringr::str_match_all(file, "q/([^_]+)_TEQUILA-seq_?(.*h)?_(.).sf"),
  sample = vapply(out, `[`, FUN.VALUE = character(1), ,2) |> {function(x) ifelse(x=="Brain", "brain","neuroblastoma")}(),
  time_points = vapply(out, `[`, FUN.VALUE = character(1), ,3) |> sub('h','Hr', x = _),
  reps = vapply(out, `[`, FUN.VALUE = character(1), ,4) |> as.integer(),
  id = sprintf("%s%s_rep%i", sample,ifelse(is.na(time_points), "", time_points), as.integer(reps)),
  out = NULL) 


sf_data$file <- NULL
sf_data <- sf_data |>
  dplyr::left_join(y = meta, by = c("sample","time_points","reps"))

GR <- gr$GR

makeSEfromTidy <- function(sf_data, assay_data, rowData, colData, GR) {
  

  box::use(dplyr[...], tidyr[...], methods[as])
  rowData_ <- select(sf_data, Name, {{rowData}}) |> distinct()
  if (!missing(GR)) {
    .n <- names(GR)
    if (!all(purrr::map2_lgl(.n, rowData_$Name, ~grepl(.x, .y)))) {
      stop("names(GR) do not match distinct rowData")
    }
    GR@elementMetadata[,names(rowData_)] <- rowData_
    sf_data <- left_join(sf_data, data.frame(Name = rowData_$Name, .name = .n, stringsAsFactors = F), by = "Name")
    sf_data$Name <- sf_data$.name
    sf_data$.name <- NULL
  }
  
  colData_ <- select(sf_data, {{colData}}) |> distinct()
  
  indexes_cols <- vapply(colData_, function(x) anyDuplicated(x)==0, logical(1))
  
  col_names_from <- if (any(indexes_cols)) {
    as.name(names(colData_)[which(indexes_cols)[1L]])
  } else {
    sf_data <- left_join(
      sf_data,
      mutate(
        colData_,
        col_names_from = sprintf("C%i", 1:n())
      )
    )
    as.name("col_names_from")
  }
  Assays_pre <- select(sf_data, Name, {{assay_data}}, !!col_names_from)
  Assays <- vector('list', ncol(Assays_pre)-2L)
  value_nms <- select(Assays_pre, -c(Name, !!col_names_from)) |> names()
  for (i in seq_along(Assays)) {
    v <- as.name(value_nms[i])
    a <- pivot_wider(Assays_pre, id_cols = Name, names_from = !!col_names_from, values_from = !!v) 
    rnames <- a$Name
    a$Name <- NULL
    a <- as.matrix(a)
    rownames(a) <- rnames
    Assays[[i]] <- a
  }
  
  names(Assays) <- value_nms
  
  
  se <- if(!missing(GR)) {
    SummarizedExperiment::SummarizedExperiment(
      assays = Assays,
      colData = as(colData_, "DataFrame"),
      rowRanges = GR
    )
  } else {
    SummarizedExperiment::SummarizedExperiment(
      assays = Assays,
      rowData = as(rowData_, "DataFrame"),
      colData = as(colData_, "DataFrame")
    )
  }
  

  se
 
}


se <- makeSEfromTidy(sf_data, c(TPM, NumReads), rowData = c(Length, EffectiveLength), colData = sample:source_name, GR = GR)

mutate(se, counts = round(NumReads,0))
assay(se, 'counts') <- round(assay(se, 2), 0)
assays(se) <- assays(se)[c(3,1,2)]
se_sub <- se[rowSums(assay(se, 'counts'))>0,]@rowRanges$gene_names
se_samp <- sample(se_sub, 10)
se_ <- se[se@rowRanges$gene_names %in% se_sub,]
se_2 <- se_[c(31407, 31408, 31409, 31410, 31411, 31412, 31413, 31414, 31415, 31416, 31417, 31418, 31419, 31420, 31421, 31422, 31423, 31424, 31425, 31426, 31427, 31428, 31429, 31430),]
tse_ <- tse(se_)
tse_expand <- tse_by_nodes2(tse_)
se_small <- se[grep('LINC01128', names(se@rowRanges)),]
clst <- cluster_data(se_sub)


se_groupings <- function(se) {
  tbl <- tibble(
    seqnames = as.factor(se@rowRanges@seqnames),
    gene_names = factor(se@rowRanges$gene_names, levels = unique(se@rowRanges$gene_names))
  ) |>
    group_by(seqnames, gene_names)
  
  tbl_grp <- group_data(tbl)
  tbl_grp
  
}

gen_tse <- function(se, ...) {
  f <- S7::method(cluster_data, getClass("RangedSummarizedExperiment"))
  hclust_obj <- f(se, ...)
  tse <- as(se, "TreeSummarizedExperiment")
  TreeSummarizedExperiment::rowTree(tse) <- ape::as.phylo(hclust_obj)
  tse
}

tse <- function(se, ...) {
  # browser()
  tbl_grp <- se_groupings(se)
  f <-  S7::method(cluster_data, getClass("RangedSummarizedExperiment"))
  n <- nrow(tbl_grp)
  rows_ <- tbl_grp$.rows
  sizes <- vapply(rows_, length, integer(1L))
  nse <- nrow(se)
  treeName <- sprintf("%s|%s", tbl_grp$seqnames, tbl_grp$gene_names)
  
  lnk_data <- LinkDataFrame(nodeLab = rownames(se),
                            nodeLab_alias = "node_1",
                            nodeNum = integer(nse),
                            isLeaf = TRUE,
                            whichTree = "phylo_null")
  lnk_data_sub <- lnk_data[,c("nodeLab_alias", "nodeNum", "whichTree")]
  phylo_empty <- structure(list(),
                           class = c("phylo_null","phylo"))
  trees <- vector("list", n)
  names(trees) <- treeName
  to_cluster <- which(sizes > 1L)
  try_fetch({
    for (i in cli::cli_progress_along(to_cluster)) {
      i_ <- to_cluster[i]
      i_seq <- rows_[[i_]]
      i_n <- length(i_seq)
      # if (i_n==1L) next
      clust <- f(se[i_seq], ...)
      trees[[i_]] <- ape::as.phylo.hclust(clust)
      tname <- treeName[i_]
      o <- clust$order
      node_num <- (1:i_n)[order(o)]
      alias <- sprintf(
        paste0("node_%0", floor(log10(2*i_n - 1L)) + 1L, "i"),
        node_num
      )
      #note - in previous algorithms, nodes were ordered according to
      # hclust --> hclust$label[hclust$order], then indexed by height...
      # the phylo class is an inigma to me. I do not know how to encode the nodes
      # yet.
      lnk_data_sub[i_seq,] <- list(alias, node_num, tname)
    }
  }, error = function(cnd) {
    abort(c(sprintf("process failed on %ith step", i),
            "i" = sprintf("i = %i", i_),
            "i" = sprintf("indexes = %s", paste(i_seq, collapse = ", ")),
            "i" = sprintf("chrom = %s; gene = %s",tbl_grp$seqnames[i_], tbl_grp$gene_names[i_])), parent = cnd)
  })
  
  lnk_data[,c("nodeLab_alias", "nodeNum", "whichTree")] <- lnk_data_sub
  trees_final <- c(trees[vapply(trees, Negate(is.null), logical(1L))], list(phylo_null = phylo_empty))
  tse <- as(se, "TreeSummarizedExperiment")
  tse@rowTree <- trees_final
  tse@rowLinks <- lnk_data
  tse
}

ob <- tse(se_sub)
ob1 <- ob[ob@rowLinks$whichTree=="chr1|ENSG00000230021",]


expand_node <- function(tse) {
  #we expect this function to be called
  # on a tse object with 1 tree in it
  # not checking is done...
  tree_name <- names(tse@rowTree)
  if (tree_name=="phylo_null") {
    gr <- tse@rowRanges
    tse@rowRanges@elementMetadata <- gr@elementMetadata[,"gene_names", drop = FALSE]
    tse@rowRanges$descendants <- list(1L)
    tse@rowRanges$exons <- list(gr$exons)
    return(tse)
  }
  if (length(tree_name)!=1L) {
    stop("expand_node() called on unexpected data")
  } 
  clust <- rowTree(tse) |> as.hclust()
  node_data <- gen_node_data(clust)
  
  n_cols <- ncol(tse)
  col_names <- colnames(tse)
  n_tse <- nrow(tse)
  n_rows <- nrow(node_data)
  col_fun <- base::colSums
  n_digits <- floor(log10(nrow(node_data)))+1
  sprintf_fmt <- paste0(tree_name,"|node_%0", n_digits, "i")
  old_assays <- assays(tse)
  n_assays <- length(old_assays)
  NewAssays <- list(
    matrix(0, nrow = n_rows,
           ncol = n_cols,
           dimnames = list(
             sprintf(sprintf_fmt, 1:n_rows),
             col_names
           ))
   ) |> rep_len(length.out = length(assays(tse)))
  
  old_row_seq <- 1:n_tse
  node_order <- order(tse@rowLinks$nodeNum)
  assay_seq <- seq_along(old_assays)
  node_seq <- (n_tse + 1L):n_rows
  descen <- node_data$descendants
  for (i in assay_seq) {
    old_assay <- old_assays[[i]]
    NewAssay <- NewAssays[[i]]
    NewAssay[old_row_seq,] <- old_assay[node_order,]
    for (j in node_seq) {
      indx <- descen[[j]]
      NewAssay[j,] <- col_fun(old_assay[indx,]) 
    }
    NewAssays[[i]] <- NewAssay
  }
  node_num <- 1:n_rows
  n_internal <- n_rows - n_tse
  lnks <- LinkDataFrame(
    nodeLab = sprintf(sprintf_fmt, node_num),
    nodeLab_alias = sprintf(sprintf("node_%%0%ii", n_digits), node_num),
    nodeNum = node_num,
    isLeaf = c(rep_len(TRUE, n_tse), rep_len(FALSE, n_internal)),
    whichTree = tree_name
  )
  out <- lapply(node_data$descendants,
                             function(i, desc) {
                               d <- desc[i]
                               node_ord <- order(d)
                               d[node_ord]},
                             desc = tse@rowLinks$nodeNum)
  #lnks$descendants <- out
  
  orig_GR <- tse@rowRanges[node_order]
  gr_seqname <- as.character(unique(tse@rowRanges@seqnames))
  gr_genename <- as.character(unique(tse@rowRanges$gene_names))
  gr_start <- integer(n_internal)
  gr_end <- integer(n_internal)
  gr_name <- character(n_internal)
  gr_strand <- character(n_internal)
  gr_strand[] <- "*"
  for (i in seq_along(node_seq)) {
    
    indx <- out[[node_seq[i]]]
    gr <- orig_GR[indx]
    gr_start[i] <- min(gr@ranges@start)
    gr_end[i] <- max(gr@ranges@width - 1L + gr@ranges@start)
    gr_strand_ <- unique(gr@strand)
    if (length(gr_strand_)==1) {
      gr_strand[i] <- as.character(gr_strand_)
    }
  }
  
  
  GRNodes <- GenomicRanges::GRanges(gr_seqname,
                               ranges = IRanges::IRanges(start = gr_start,
                                                         end = gr_end,
                                                         name = gr_name),
                               strand = gr_strand)
  gr_exons <- lapply(out, function(i, GR) {GR[i]$exons}, GR = orig_GR)
  orig_GR@elementMetadata <- orig_GR@elementMetadata[,"gene_names", drop = FALSE]
  newGR <- c(orig_GR, GRNodes)
  newGR$gene_names <- gr_genename
  newGR$descendants <- out
  newGR$exons <- gr_exons
  names(newGR) <- sprintf(sprintf_fmt, 1:n_rows)
  names(NewAssays) <- names(old_assays)
  
  out_tse <- SummarizedExperiment(
    assays = NewAssays,
    rowRanges = newGR
  ) |> as("TreeSummarizedExperiment")
  out_tse@rowTree <- tse@rowTree
  out_tse@rowLinks <- lnks

  out_tse
  
}

tse_by_nodes2 <- function(tse, BPPARAM = BiocParallel::MulticoreParam()) {
  tbl_grp <- se_groupings(tse)
  n_grps <- nrow(tbl_grp)
  out <- vector('list', n_grps)
  sizes <- vapply(tbl_grp$.rows, length, integer(1)) * 2 - 1L
  n_row <- sum(sizes)
  row_names <- character(n_row)
  col_names <- colnames(tse)
  k <- 0L
  # box::use(progressr[progressor])
  # seq_i <- seq_len(n_row)

  
  # progressr::with_progress({
  #   p <- progressor(along = seq_i)
  #   out <- BiocParallel::bplapply(
  #     seq_i,
  #     function(i, tse_obj, rows_, prog) {
  #       out <- expand_node(tse_obj[rows_[[i]],])
  #       prog()
  #       out
  #     },
  #     tse_obj = tse,
  #     rows_ = tbl_grp$.rows,
  #     prog = p,
  #     BPPARAM = BPPARAM
  #   )
  # })
  
  
  for (i in cli::cli_progress_along(out, name = "expanding")) {
    k <- k + 1L
    k_ <- k + sizes[i] - 1L
    tse_tmp <- expand_node(tse[tbl_grp$.rows[[i]],])
    row_names[k:k_] <- rownames(tse_tmp)
    out[[i]] <- tse_tmp
    k <- k_
  }
  
  oldAssays <- assays(tse)
  n_assays <- length(oldAssays)
  seq_assay <- seq_len(n_assays)
  NewAssays <- list(
    matrix(0, nrow = n_row,
           ncol = ncol(tse),
           dimnames = list(
             row_names,
             col_names
           ))
  ) |> rep_len(length.out = n_assays)
  names(NewAssays) <- names(oldAssays)
  
  NewLink <- LinkDataFrame(nodeLab = character(n_row),
                           nodeLab_alias = character(n_row),
                           nodeNum = integer(n_row),
                           isLeaf = TRUE,
                           whichTree = "phylo_null")
  NewGR <- GenomicRanges::GRanges(
    vec_rep(as.character(tbl_grp$seqnames), sizes),
    ranges = IRanges::IRanges(start = integer(n_row),
                              end = integer(n_row),
                              name = ""),
    gene_names = vec_rep(as.character(tbl_grp$gene_names), sizes),
    strand = "*")
  NewGR@elementMetadata$descendants <- rep_len(list(""), n_row)
  NewGR@elementMetadata$exons <- rep_len(list(""), n_row)
  names(NewGR) <- row_names
  k <- 0L
  for (i in cli::cli_progress_along(out, name = "collecting")) {
    k <- k + 1L
    k_ <- k + sizes[i] - 1L
    seq_k <- k:k_
    tmp_tse <- out[[i]]
    tmp_assays <- tmp_tse@assays@data
    for (j in seq_assay) {
      NewAssays[[j]][seq_k,] <- tmp_assays[[j]]
    }
    NewLink[seq_k,] <- tmp_tse@rowLinks
    NewGR[seq_k] <- tmp_tse@rowRanges
    k <- k_
  }
  
  trees <- tse@rowTree
  
  tse_out <- SummarizedExperiment(
    assays = NewAssays,
    rowRanges = NewGR
  ) |> as("TreeSummarizedExperiment")
  tse_out@rowLinks <- NewLink
  tse_out@rowTree <- trees
  tse_out
}

print.phylo_null <- function(x) {
  cat("A placeholder Phylogenic object\n")
}


