
box::use(tidyr[pivot_wider],
         dplyr[select, distinct, left_join],
         methods[as])

#' As SummerizedExperiment
#' @name as_se
#' @param .data an object
#' @export
as_se <- function(.data, ...) UseMethod("as_se")

#' @describeIn as_se
#' Converts a data.frame into a SummarizedExperiment Object.
#' ColumnData and RowData are inferred from data that can be merged
#' into the genes/cols columns uniquely. All other data, unless assigned
#' are ignored.
#' @param genes The column that will become the genes
#' @param cols The column that will become the column
#' @param vals The columns reshaped into matrix data
#' @param ... Unused
as_se.data.frame <- function(.data, genes, cols, vals, ...) {
  box::use(S4Vectors[...])
  obj <- tidyselect::eval_select(rlang::enquo(vals), .data)
  Assays <- vector("list", length(obj))
  for ( i in seq_along(Assays)) {
    Assays[[i]] <- pivot_wider(.data, id_cols = {{genes}}, names_from = {{cols}}, values_from = names(obj)[i])
  }
  Assays <- lapply(Assays, function(x) {r <- x[[1]]; x <- as.matrix(x[,-1L, drop = F]); rownames(x) <- r; x})
  names(Assays) <- names(obj)
  .tmp <- select(.data, -all_of(names(obj)), - {{genes}})
  
  colData <- select(.tmp, {{cols}}) |> distinct() 
  .col <- names(colData)
  .I <- nrow(colData)
  
  for (nm in setdiff(names(.tmp), .col)) {
    .new <- left_join(colData, distinct(.tmp[,c(.col,nm), drop = F]), by = .col)
    if (nrow(.new)==.I) {
      colData <- .new
    }
  }
  colData <- as(as.data.frame(colData), "DataFrame")
  
  .tmp <- select(.data, -all_of(names(obj)), - {{cols}}) 
  rowData <- select(.tmp, {{genes}}) |> distinct()
  .row <- names(rowData)
  .I <- nrow(rowData)
  
  for (nm in setdiff(names(.tmp), .row)) {
    .new <- left_join(rowData, distinct(.tmp[,c(.row,nm), drop = F]), by = .row)
    if (nrow(.new)==.I) {
      rowData <- .new
    }
  }
  rowData <- as(as.data.frame(rowData), "DataFrame")
  SummarizedExperiment::SummarizedExperiment(assays = Assays, rowData = rowData, colData = colData)
  
}

as_se.grouped_df <- function(.data, genes, cols, vals, ...) {
  .data <- dplyr::ungroup(.data)
  as_se.data.frame(.data, {{genes}}, {{cols}}, {{vals}})
}

as_se.default <- function(.data, ...) stop(paste0("No method for object of class <",paste0(class(.data), collapse = "/"),">"), call. = F)