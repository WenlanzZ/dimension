#' IPF Single Cell Altas data
#'
#' The lung data is a subset of patient "001C" profiled in IPF
#' Single Cell Altas. The full data seet can be found on GEO (GSE136831),
#' which was examined  by
#' \href{https://www.biorxiv.org/content/10.1101/759902v1.full.pdf}
#' to build a single cell atlas of Idiopathic Pulmonary Fibrosis (IPF).
#' It contains the 60603 gene expression of 158 cells profiled from
#' "distal lung parenchyma samples obtained from 32 IPF, 18 chronic obstructive
#' pulmonary disease (COPD) and 29 control donor lungs" (Adamset al.2019).
#'
#' @docType data
#'
#' @usage data(lung)
#'
#' @format An object of class Seurat.
#'
#' @keywords datasets
#'
#' @references Adamset  al.(2019),
#' (\href{https://www.biorxiv.org/content/10.1101/759902v1.full.pdf})
#'
#' @source \href{https://www.biorxiv.org/content/10.1101/759902v1.full.pdf}
#'
#' @examples
#' \dontrun{
#' data(lung)
#' }
"lung"
