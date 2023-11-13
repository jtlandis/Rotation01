
#' Useful RNAseq Usage Patterns
#' @name .__module__.
#' @description 
#' helpful functions to perform rnaseq in R. Generally
#' we will use the DESeq2 package on raw counts (unnormalized).
#' 
#' @examples 
#' 
#' #settup
#' set.seed(123)
#' tbl <- tidyr::expand_grid(
#'    data.frame(rows = sprintf('row%04d',1:1000),
#'               lens = sample(100:200, 1000, T)),
#'    data.frame(cols = sprintf('col%04d', 1:100),
#'               trmt = sample(c("cntrl","drug"), 100, T)))
#' 
#' tbl <- dplyr::group_by(tbl, rows, cols) %>%
#'   dplyr::mutate(counts = as.integer(round(runif(dplyr::n(),lens/2, lens) * ifelse(trmt=="drug", runif(1, 1, 2), 1),0)))
#' 
#' #convert to Summarized Experiment from data.frame easily
#' se <- as_se(tbl, genes = rows, cols = cols, vals = counts)
#' 
#' dds <- DESeqDataSet(se, design = ~ trmt)
#' dds <- DESeq(dds)
#' res <- results(dds, name = "trmt_drug_vs_cntrl")
#' res
#' 
#' @export
box::use(
  DESeq2[...],
  SummarizedExperiment[colData, rowData, `colData<-`, `rowData<-`, assay, `assay<-`],
  ./utils[...]
)