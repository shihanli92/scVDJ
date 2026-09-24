#' Assign clones to cells
#'
#' Groups cells with the same receptor into clones. A cell's receptor is the set
#' of all its contigs, each identified by its chain and CDR3 sequence (plus its
#' gene calls when \code{match_genes = TRUE}); contig order within a cell does
#' not matter, and contigs with identical keys count once. Cells must match on
#' every contig to share a clone, so a cell with two different TRA chains is not
#' in the same clone as a cell with only one of them.
#'
#' Clones are called separately within each combination of \code{group_cols}
#' (e.g. per donor), so identical receptors in different groups get different
#' clone ids. A cell is identified by \code{group_cols} plus its barcode, so
#' barcodes may repeat across groups. Contigs with a missing or empty CDR3 are
#' ignored; cells with none left get \code{NA}.
#'
#' Input can be contig-level or cell-level. Contig-level data has one row per
#' contig, such as the \code{contigs} slot of a \code{\linkS4class{VDJ}}
#' object. Cell-level data has one row per cell with \code{_vj} and \code{_vdj}
#' columns, such as \code{combineChains} metadata or Seurat \code{meta.data}
#' built from it; multiple chains joined with \code{"|"} are split, each row is
#' one cell, and \code{cell_id_col} is not used. Cell-level data is detected by
#' \code{chain_vj} and \code{chain_vdj} columns without a \code{chain} column.
#'
#' @param contigs Contig-level or cell-level data frame (see Details).
#' @param cdr3 Match CDR3 by amino acid (\code{"aa"}) or nucleotide
#'   (\code{"nt"}) sequence.
#' @param match_genes If \code{TRUE}, contigs must also match on every column in
#'   \code{gene_cols}.
#' @param group_cols Columns to group cells by before calling clones, e.g.
#'   \code{c("donor")}. Each cell must have a single value per column.
#' @param clone_col Name of the clone id column to add.
#' @param prefix Prefix for clone ids, e.g. \code{"clone_"} gives
#'   \code{"clone_1"}.
#' @param count_col Name of the column to add with the number of cells in the
#'   clone.
#' @param cell_id_col,chain_col Cell barcode and chain columns. For cell-level
#'   input these are base names, e.g. \code{"chain"} for \code{chain_vj} and
#'   \code{chain_vdj}; the same applies to the CDR3 and gene columns.
#' @param junction_aa_col,junction_col Amino acid and nucleotide CDR3 columns.
#' @param gene_cols Gene call columns matched when \code{match_genes = TRUE}.
#'   Drop \code{"c_call"} to keep class-switched BCRs in the same clone.
#'   \code{combineChains} metadata has no D calls, so use
#'   \code{c("v_call", "j_call", "c_call")} for cell-level input.
#'
#' @return The input with \code{clone_col} and \code{count_col} added, rows
#'   and row names unchanged; for contig-level input every contig of a cell gets
#'   that cell's clone. Existing columns with those names are replaced.
#'
#' @examples
#' contigs <- data.frame(
#'   barcode = c("A", "A", "B", "B", "C"),
#'   donor = c("d1", "d1", "d1", "d1", "d2"),
#'   chain = c("TRA", "TRB", "TRB", "TRA", "TRB"),
#'   junction_aa = c("CAVF", "CASSF", "CASSF", "CAVF", "CASSF"),
#'   junction = c("tgt1", "tgt2", "tgt2", "tgt1", "tgt2"),
#'   v_call = "V", d_call = "", j_call = "J", c_call = "C")
#' callClones(contigs, cdr3 = "aa", group_cols = "donor")
#'
#' # Cell-level, e.g. seurat_obj@meta.data <- callClones(seurat_obj@meta.data)
#' cells <- data.frame(chain_vj = c("TRA", "TRA|TRA", NA),
#'                     junction_aa_vj = c("CAVF", "CAVF|CAGF", NA),
#'                     chain_vdj = c("TRB", "TRB", NA),
#'                     junction_aa_vdj = c("CASSF", "CASSF", NA),
#'                     row.names = c("cell1", "cell2", "cell3"))
#' callClones(cells)
#'
#' @export
callClones <- function(contigs,
                       cdr3 = c('aa', 'nt'),
                       match_genes = FALSE,
                       group_cols = NULL,
                       clone_col = 'clone_id',
                       prefix = 'clone_',
                       count_col = paste0(clone_col, '_count'),
                       cell_id_col = 'barcode',
                       chain_col = 'chain',
                       junction_aa_col = 'junction_aa',
                       junction_col = 'junction',
                       gene_cols = c('v_call', 'd_call', 'j_call', 'c_call')){

    cdr3 <- match.arg(cdr3)
    cdr3_col <- if(cdr3 == 'aa') junction_aa_col else junction_col
    key_cols <- c(chain_col, cdr3_col, if(match_genes) gene_cols)
    data <- contigs[, setdiff(colnames(contigs), c(clone_col, count_col)), drop = FALSE]

    # Per-cell input (combineChains metadata): <col>_vj / <col>_vdj columns, one row per cell
    cell_level <- !chain_col %in% colnames(data) &&
        all(paste0(chain_col, c('_vj', '_vdj')) %in% colnames(data))

    if(cell_level){
        side_cols <- c(paste0(key_cols, '_vj'), paste0(key_cols, '_vdj'))
        check_cols(c(side_cols, group_cols), data)
        # Each row is a cell; split pipe-joined chains into one contig key each
        cell_ids <- seq_len(nrow(data))
        long <- do.call(rbind, lapply(c('_vj', '_vdj'), function(suffix){
            parts <- lapply(data[paste0(key_cols, suffix)], function(x)
                stringr::str_split(as.character(x), stringr::fixed('|')))
            n <- lengths(parts[[1]])
            data.frame(.cell = rep(cell_ids, n),
                       .cdr3 = unlist(parts[[2]]),
                       .contig_key = do.call(paste, c(lapply(unname(parts), unlist), sep = ':')))
        }))
    } else {
        check_cols(c(key_cols, group_cols, cell_id_col), data)
        # A cell is its group_cols plus barcode, so barcodes may repeat across groups
        cell_ids <- do.call(paste, c(unname(as.list(data[c(group_cols, cell_id_col)])), sep = '\r'))
        long <- data.frame(.cell = cell_ids,
                           .cdr3 = data[[cdr3_col]],
                           .contig_key = do.call(paste, c(unname(as.list(data[key_cols])), sep = ':')))
    }
    long <- long[!is.na(long$.cdr3) & long$.cdr3 != '', ]
    long <- cbind(long, as.data.frame(data[match(long$.cell, cell_ids), group_cols, drop = FALSE]))

    # One receptor key per cell: its unique contig keys sorted (radix = locale independent) and joined
    cells <- long %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(c('.cell', group_cols)))) %>%
        dplyr::summarise(.receptor = paste(sort(unique(.data$.contig_key), method = 'radix'), collapse = '|'),
                         .groups = 'drop') %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(c(group_cols, '.receptor')))) %>%
        dplyr::mutate(.clone = paste0(prefix, dplyr::cur_group_id()),
                      .count = dplyr::n()) %>%
        dplyr::ungroup()

    # Assign by position so row order, row names and data frame class are kept
    idx <- match(cell_ids, cells$.cell)
    data[[clone_col]] <- cells$.clone[idx]
    data[[count_col]] <- cells$.count[idx]
    return(data)
}

# Stops if any of cols is missing from data
check_cols <- function(cols, data) {
    missing_cols <- setdiff(cols, colnames(data))
    if(length(missing_cols) > 0){
        stop('Columns not found: ', paste(missing_cols, collapse = ', '))
    }
}
