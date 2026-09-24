#' Keep only the top contig on each side
#'
#' Takes the first entry of every pipe-joined metadata column, which is the
#' highest-umi contig, and converts the result to its natural type (e.g. umis
#' become numeric). Row names are kept.
#'
#' @param object A \code{VDJ} object with metadata from \code{combineChains}.
#' @return The object with one value per metadata cell.
#' @noRd
keep_high_umi <- function(object) {

    # Keep first (highest umi) entry of pipe-joined columns; df[] preserves row names
    object@metadata[] <- lapply(object@metadata, function(x){
        if(!is.character(x)) return(x)
        utils::type.convert(stringr::str_split_i(x, pattern = '\\|', i = 1), as.is = TRUE)
    })

    return(object)

}


#' Split one chain side into one row per contig
#'
#' @param metadata Metadata from \code{combineChains}, one row per cell.
#' @param suffix \code{"_vj"} or \code{"_vdj"}.
#' @return A long data frame keyed by \code{barcode} with one row per contig on
#'   that side and the \code{pipe_cols} (suffixed) as columns; umis are numeric.
#'   A cell missing that chain type gets a single \code{NA} row.
#' @noRd
split_side <- function(metadata, suffix) {
    cols <- paste0(pipe_cols, suffix)
    # str_split keeps empty fields (e.g. missing c_call) so all columns stay the same length
    parts <- lapply(metadata[cols], function(x) stringr::str_split(as.character(x), stringr::fixed('|')))
    side <- data.frame(barcode = rep(metadata$barcode, lengths(parts[[1]])))
    for(col in cols){
        side[[col]] <- unlist(parts[[col]])
    }
    side[[paste0('umis', suffix)]] <- as.numeric(side[[paste0('umis', suffix)]])
    return(side)
}


#' Expand cells into every VJ x VDJ contig pairing
#'
#' @param metadata Metadata from \code{combineChains}, one row per cell.
#' @return A data frame with one row per pairing: the cell-level columns
#'   repeated, one value per pipe-joined column, numeric umis and
#'   \code{sum_umis}, the pairing's combined umi count. A cell missing one chain
#'   type gets \code{NA} for that side.
#' @noRd
expand_cell_clones <- function(metadata) {

    vj <- split_side(metadata, '_vj')
    vdj <- split_side(metadata, '_vdj')
    split_cols <- c(paste0(pipe_cols, '_vj'), paste0(pipe_cols, '_vdj'))
    cell_cols <- metadata[, setdiff(names(metadata), split_cols), drop = FALSE]

    expanded <- cell_cols %>%
        dplyr::inner_join(vj, by = 'barcode') %>%
        dplyr::inner_join(vdj, by = 'barcode', relationship = 'many-to-many')
    expanded$sum_umis <- rowSums(expanded[, c('umis_vj', 'umis_vdj')], na.rm = TRUE)

    return(expanded)

}
