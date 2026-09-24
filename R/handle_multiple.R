keep_high_umi <- function(object) {

    # Keep first (highest umi) entry of pipe-joined columns; df[] preserves row names
    object@metadata[] <- lapply(object@metadata, function(x){
        if(!is.character(x)) return(x)
        utils::type.convert(stringr::str_split_i(x, pattern = '\\|', i = 1), as.is = TRUE)
    })

    return(object)

}


# Splits the pipe-joined columns of one chain side into long format: one row per contig,
# keyed by barcode. A cell missing that chain type gets a single NA row.
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


# Expands each metadata row (a single cell) into every pairing of its vj and vdj contigs.
# Pipe-joined columns become one value per row, umis are numeric and sum_umis is the
# combined umi count of the pairing. A cell missing one chain type gets NA for that side.
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
