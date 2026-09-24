#' VDJ class
#'
#' Holds single cell immune receptor data at two levels: every contig, and
#' chains paired per cell.
#'
#' @slot celltype Receptor type: \code{"abTCR"}, \code{"gdTCR"} or \code{"BCR"}.
#' @slot contigs One row per contig, as parsed by \code{\link{readVDJ_10X}},
#'   with AIRR column names. \code{\link{combineChains}} adds a logical
#'   \code{keep} column.
#' @slot metadata One row per cell (or per chain pairing in the expand modes),
#'   built by \code{\link{combineChains}}.
#' @slot species Species name, e.g. \code{"Homo_sapiens"}.
#' @slot clone_data List reserved for clonotype results.
#'
#' @exportClass VDJ
setClass("VDJ", slots=list(celltype='character',
                           contigs="data.frame",
                           metadata="data.frame",
                           species='character',
                           clone_data='list'))



#' Build a VDJ object from a 10x output directory
#'
#' Reads the contigs with \code{\link{readVDJ_10X}} and pairs chains per cell
#' with \code{\link{combineChains}} using \code{handle_multiple = "keep_all"}.
#'
#' @param path Path to a 10x \code{vdj_t}/\code{vdj_b} output directory, or to a
#'   \code{filtered_contig_annotations.csv} file.
#' @param celltype Receptor type: \code{"abTCR"}, \code{"gdTCR"} or \code{"BCR"}.
#' @param species Species name stored on the object.
#' @param input_origin Source of the data. Only \code{"10x"} is parsed; any
#'   other value returns an empty object.
#'
#' @return A \code{\linkS4class{VDJ}} object with \code{contigs} and
#'   \code{metadata} filled in.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom methods new
#' @export
VDJ <- function(path, celltype = 'abTCR',species='Homo_sapiens', input_origin = '10x') {

    object <- new("VDJ",
                  celltype=celltype,
                  species=species,
                  clone_data=list())

    if(input_origin == '10x'){
        object <- readVDJ_10X(object,
                               path=path)
        object <- combineChains(object,handle_multiple = 'keep_all')
    }
    return(object)
}
