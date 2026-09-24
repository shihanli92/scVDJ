# VDJ Object contains two dataframes:
# Pne which contains detailed information of all contigs for cells separately,
# while the metadata contains paired TCRS where each row is a individual cell
setClass("VDJ", slots=list(celltype='character',
                           contigs="data.frame",
                           metadata="data.frame",
                           species='character',
                           clone_data='list'))



#' Builds new VDJ object based on path location
#'
#' @param path string indicating location of datafiles to be parsed
#' @param celltype String indicating the type of immune receptors can be abTCR, gdTCR or BCR
#' @param keep_high_umi Boolean indicating whether to show only most abundant gene for each cell
#'
#' @return Returns a VDJ object
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
