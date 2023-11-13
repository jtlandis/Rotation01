
box::use(
  utils[read.table],
  stats[hclust, dist],
  stringr[str_extract],
  tidyr[unnest],
  dplyr[select, as_tibble],
  GenomicRanges[GRanges],
  IRanges[IRanges],
  S7[...],
  mods/util[vec_rep]
  
)



#' @param x vector of integer indicies (ordered)
#' @param y vector of integer indicies (ordered)
#' @description x and y should be disjoint and are assumed
#' to be ordered.
#' @keywords internal
#' @return list length of x with indicies of y embeded.
embed_yOnx <- function(x, y) {
  y_max <- max(y) + 1L
  x_out <- vector("list", length(x))
  x <- c(x, y_max)
  y <- c(y, y_max)
  i <- 1L
  lower <- x[2L]
  x_i <- 1L
  for (y_i in seq_along(y)) {
    if (y[y_i] >= lower) {
      x_out[[x_i]] <- y[i:(y_i - 1L)]
      x_i <- x_i + 1L
      i <- y_i
      lower <- x[x_i + 1L]
    }
  }
  x_out
}

#' @export
read_gtf <- function(file, ..., header = F, nrows = Inf) {
  gtf <- readr::read_delim(file = file, delim = "\t", col_names = header, show_col_types = F,
                           progress = T, comment = "#", n_max = nrows, ...)
  
  names(gtf) <- c("seqname", "source", "feature",
                  "start", "end", "score",
                  "strand", "frame", "attribute")
  
  gene_ind <- grep("gene", gtf$feature)
  tran_ind <- grep("transcript", gtf$feature)
  exon_ind <- grep("exon", gtf$feature)
  
  tran_out <- embed_yOnx(tran_ind, exon_ind)
  gene_out <- embed_yOnx(gene_ind, tran_ind)
  
  gene_df <- gtf[gene_ind,]
  
  tran_df <- gtf[tran_ind,]
  tran_df$exons <- lapply(tran_out, FUN = function(i, x) {
    y <- x[i,]
    y$name <- extract_name(y$attribute, type = "exon")
    class(y) <- c("exon_df",class(y))
    y
    }, x = gtf)
  tran_df$n_exons <- vapply(tran_out, length, 1L)
  n <- vapply(gene_out, length, 1L)
  gene_df$trans <- lapply(mapply(`+`,lapply(n,seq_len),c(0L, cumsum(n)[-length(n)])),
                          function(i, x) {
                            y <- x[i,]
                            y$name <- extract_name(y$attribute, type = "transcript")
                            class(y) <- c("transcript_df", class(y))
                            y}, x = tran_df)
  gene_df$n_trans <- n
  gene_df$name <- extract_name(gene_df$attribute, type = "gene")
  class(gene_df) <- c("gene_df", class(gene_df))
  gene_df
  
}




#' @export
read_gtf2GR <- function(file, ngenes = Inf) {
  box::use(GenomicRanges[GRanges, GRangesList, split],
           IRanges[IRanges])
 
  cat("Reading data...\n")
  gtf <- readr::read_delim(
    file = file, delim = "\t", comment = "#",
    progress = T, col_names = c("seqname", "source", "feature",
                                "start", "end", "score",
                                "strand", "frame", "attribute"),
    show_col_types = F, lazy = T
  )
  cat("Indexing Features...")
  gene_ind <- grep("gene", gtf$feature)
  tran_ind <- grep("transcript", gtf$feature)
  exon_ind <- grep("exon", gtf$feature)
  cat("\rFeature Indexing Complete...\n")
  
  if (is.finite(ngenes)) {
    g_limit <- gene_ind[ngenes+1L]
    gene_ind <- gene_ind[1:ngenes]
    tran_ind <- tran_ind[tran_ind<g_limit]
    exon_ind <- exon_ind[exon_ind<g_limit]
    
  }
  
  gtf$name <- character(nrow(gtf))
  cat("Extracting names...")
  gtf$name[gene_ind] <- extract_name(gtf$attribute[gene_ind], "gene")
  cat(".")
  gtf$name[tran_ind] <- extract_name(gtf$attribute[tran_ind], 'transcript')
  cat(".")
  gtf$name[exon_ind] <- extract_name(gtf$attribute[exon_ind], 'exon')
  cat(".")
  cat("\rFeature Names Extracted...\nEmbedding...")
  tran_out <- embed_yOnx(tran_ind, exon_ind)
  gene_out <- embed_yOnx(gene_ind, tran_ind)
  cat("\rIndicies Embeded...\n")
  genes <- gtf$name[gene_ind] 
  GR_n <- length(tran_ind)
  gene_names <- vec_rep(genes, vapply(gene_out, length, integer(1L)))
  trans_names <- gtf$name[tran_ind]
  cat("Forming GenomicRanges...\n")
  
  GR <- GRanges(seqnames = gtf$seqname[tran_ind],
                ranges = IRanges(start = gtf$start[tran_ind],
                                 end = gtf$end[tran_ind],
                                 name = trans_names),
                gene_names = gene_names,
                strand = gtf$strand[tran_ind])
  cat("Installing exons...")
  gtf_sub <- dplyr::select(gtf, seqname, name, start, end, strand)
  exons_sub <- gtf_sub[exon_ind,]
  exons_GR <- as_genomicRange_(exons_sub)
  exons_GR <- GenomicRanges::split(
    exons_GR,
    factor(vec_rep(trans_names,
                   vapply(tran_out, length, integer(1L))),
           trans_names)
    )
  GR$exons <- exons_GR
  cat("... DONE!\n")
  on.exit(gc())
  GR
  
}

#' @export
extract_attribute <- function(str, attribute) {
  
  str_extract(str, sprintf("%s [^;]+", attribute)) |>
    sub(pattern = sprintf("%s ", attribute), replacement = "", x = _)
  
}

#' @export
extract_name <- function(x, type = c("gene", "transcript", "exon")) {
  name <- match.arg(type, c("gene", "transcript", "exon"))
  name <- switch(name,
                 gene = "gene_name",
                 transcript = "transcript_name",
                 exon = "exon_id")
  labels <- str_extract(x, sprintf("%s [^;]+", name)) |>
    sub(sprintf("%s ", name), "", x = _) |>
    gsub('"', '', x = _)
  labels
}

#' @export
as_genomicRange_ <- function(tbl, ...) {
  
  GRanges(seqnames = tbl$seqname,
          ranges = IRanges(start = tbl$start,
                           end = tbl$end,
                           name = tbl$name),
          ...,
          strand = tbl$strand)
  
}

#' @export
as_genomicRange <- function(x) UseMethod('as_genomicRange')

as_genomicRange.gene_df <- function(x) {
  
  as_genomicRange_(
    x,
    score = x$score,
    frame = x$frame,
    trans = do.call(
      GenomicRanges::GRangesList,
      lapply(x$trans, as_genomicRange.transcript_df)
    )
  )
  
}

as_genomicRange.transcript_df <- function(x) {
  
  as_genomicRange_(
    x,
    score = x$score,
    frame = x$frame,
    exons = do.call(
      GenomicRanges::GRangesList,
      lapply(x$exons, as_genomicRange.exon_df)
    )
  )
  
}

as_genomicRange.exon_df <- function(x) {
  
  as_genomicRange_(
    x,
    score = x$score,
    frame = x$frame
  )
  
}

## Classes ----

#' @export
gene_df <- new_S3_class("gene_df")

#' @export
transcript_df <- new_S3_class("transcript_df")

#' @export
exon_df <- new_S3_class("exon_df")
