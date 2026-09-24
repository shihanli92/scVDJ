#' @include vdj.R
NULL

#' Pair chains per cell
#'
#' Collapses the contigs of a \code{\linkS4class{VDJ}} object into one row per
#' cell, with VJ chains (TRA, TRG, IGK/IGL) and VDJ chains (TRB, TRD, IGH) side
#' by side. Cells with only one chain type are kept, with \code{NA} for the
#' missing side.
#'
#' Within each side, contigs are ordered by decreasing UMIs, with ties broken
#' by \code{tie_break_col}. Metadata columns carry the suffix \code{_vj} or
#' \code{_vdj}. \code{chain}, \code{v_call}, \code{j_call}, \code{c_call},
#' \code{junction}, \code{junction_aa}, \code{umis} and \code{contig_id} hold
#' every contig in that order, joined with \code{"|"}. \code{best_contig} is the
#' top contig, and \code{n_<chain>} columns (e.g. \code{n_TRA}) count contigs
#' per chain.
#'
#' @param object A \code{\linkS4class{VDJ}} object with a filled \code{contigs} slot.
#' @param handle_multiple How to handle cells with more than one contig of a
#'   chain type:
#'   \describe{
#'     \item{\code{"keep_all"}}{One row per cell with every contig pipe-joined.}
#'     \item{\code{"keep_high_umi"}}{One row per cell keeping only the top
#'       contig on each side; UMI columns become numeric.}
#'     \item{\code{"expand_clones"}}{One row per VJ x VDJ contig pairing, with
#'       \code{sum_umis} (combined UMIs) and \code{umi_rank} (rank within the
#'       cell, 1 = highest). Row names are \code{<barcode>_<umi_rank>}.}
#'     \item{\code{"expand_clones_keep_one_clonal"}}{As \code{"expand_clones"},
#'       plus \code{clone_count} (cells sharing the pairing's nucleotide
#'       junctions) and \code{keep}, marking per cell the pairing with the
#'       highest \code{clone_count} (ties by \code{umi_rank}). Only contigs of
#'       kept pairings have \code{keep = TRUE} in the \code{contigs} slot.}
#'   }
#' @param ... Unused.
#' @param contig_id_col,cell_id_col,chain_col Contig id, cell barcode and chain
#'   columns in \code{contigs}.
#' @param variable_gene_col,joining_gene_col,constant_gene_col V, J and C gene
#'   call columns.
#' @param junction_col,junction_aa_col Nucleotide and amino acid junction
#'   (CDR3) columns.
#' @param read_count_col UMI count column used to order contigs.
#' @param tie_break_col Column used to break UMI ties, by decreasing value.
#'   Skipped if \code{NULL} or not in \code{contigs}, in which case ties keep
#'   input order.
#' @param other_cols Extra \code{contigs} columns to carry into
#'   \code{metadata}, taking the top contig's value on each side.
#'
#' @return The \code{VDJ} object with \code{metadata} filled in and a logical
#'   \code{keep} column added to \code{contigs}.
#'
#' @export
setGeneric("combineChains", function(object, handle_multiple = 'keep_all', ...) standardGeneric("combineChains"))

#' Collapse contigs of the given chains to one row per barcode
#'
#' Chains, gene calls, junctions, umis and contig ids are pipe-joined in
#' decreasing umi order, with ties broken by \code{reads} when that column is
#' present (otherwise by input order).
#'
#' @param contigs Contigs with the AIRR column names set by \code{combineChains}.
#' @param chains Chain names to keep, e.g. \code{c("IGK", "IGL")}.
#' @param other_cols Extra columns to carry through, taking the top contig's value.
#' @return A data frame with one row per barcode: \code{barcode},
#'   \code{best_contig} (top contig id), the \code{pipe_cols} and \code{other_cols}.
#' @noRd
collapse_chains <- function(contigs, chains, other_cols) {
    if(!'reads' %in% colnames(contigs)){
        contigs$reads <- 0
    }
    contigs %>%
        dplyr::filter(.data$chain %in% chains) %>%
        dplyr::arrange(dplyr::desc(.data$umis), dplyr::desc(.data$reads)) %>%
        dplyr::group_by(.data$barcode) %>%
        dplyr::summarise(best_contig = dplyr::first(.data$contig_id),
                         dplyr::across(dplyr::all_of(pipe_cols), ~ paste0(.x, collapse = "|")),
                         dplyr::across(dplyr::all_of(other_cols), dplyr::first),
                         .groups = 'drop')
}

# Columns pipe-joined per chain type by collapse_chains
pipe_cols <- c('chain', 'v_call', 'j_call', 'c_call', 'junction', 'junction_aa', 'umis', 'contig_id')

#' Count contigs of each chain type per cell
#'
#' @param object A \code{VDJ} object whose metadata has pipe-joined
#'   \code{chain_vj} and \code{chain_vdj} columns.
#' @param vj,vdj Chain names on each side, e.g. \code{"TRA"} and \code{"TRB"}.
#' @return The object with an integer \code{n_<chain>} metadata column per chain
#'   (e.g. \code{n_TRA}, \code{n_TRB}); 0 when the cell has none.
#' @noRd
count_chains <- function(object, vj, vdj) {
    count_chain <- function(collapsed, chain) {
        vapply(strsplit(collapsed, '|', fixed = TRUE),
               function(x) sum(x == chain, na.rm = TRUE), integer(1))
    }
    for(chain in vj){
        object@metadata[[paste0('n_', chain)]] <- count_chain(object@metadata$chain_vj, chain)
    }
    for(chain in vdj){
        object@metadata[[paste0('n_', chain)]] <- count_chain(object@metadata$chain_vdj, chain)
    }
    return(object)
}

#' @rdname combineChains
setMethod("combineChains", "VDJ", function(object,
                                           handle_multiple = 'keep_all', ...,
                                           contig_id_col = 'contig_id',
                                           cell_id_col ='barcode',
                                           chain_col = 'chain',
                                           variable_gene_col = 'v_call',
                                           joining_gene_col = 'j_call',
                                           constant_gene_col = 'c_call',
                                           junction_col = 'junction',
                                           junction_aa_col = 'junction_aa',
                                           read_count_col = 'umis',
                                           tie_break_col = 'reads',
                                           other_cols = c()
                                           ){

    handle_multiple <- match.arg(handle_multiple, c('keep_all', 'keep_high_umi', 'expand_clones',
                                                    'expand_clones_keep_one_clonal'))

    vj_options <- list(abTCR=c('TRA'),
                       gdTCR=c('TRG'),
                       BCR=c('IGK', 'IGL'))

    vdj_options <- list(abTCR=c('TRB'),
                        gdTCR=c('TRD'),
                        BCR=c('IGH'))

    # Select required columns and rename them to the AIRR names used below
    col_map <- c(contig_id = contig_id_col, barcode = cell_id_col, chain = chain_col,
                 v_call = variable_gene_col, j_call = joining_gene_col, c_call = constant_gene_col,
                 junction = junction_col, junction_aa = junction_aa_col, umis = read_count_col)
    missing_cols <- setdiff(c(col_map, other_cols), colnames(object@contigs))
    if(length(missing_cols) > 0){
        stop('Columns not found in contigs: ', paste(missing_cols, collapse = ', '))
    }
    # Optional column used to break umi ties, skipped if absent
    if(!is.null(tie_break_col) && tie_break_col %in% colnames(object@contigs)){
        col_map <- c(col_map, reads = tie_break_col)
    }
    contigs <- object@contigs %>%
        dplyr::select(dplyr::all_of(col_map), dplyr::all_of(other_cols))

    vj <- collapse_chains(contigs, vj_options[object@celltype][[1]], other_cols)
    vdj <- collapse_chains(contigs, vdj_options[object@celltype][[1]], other_cols)

    if(nrow(vj) == 0 && nrow(vdj) == 0){
        stop('No contigs found, check inputted cell type is correct')
    }

    # Full join keeps cells with only one chain type (missing chain is NA)
    object@metadata <- merge(x=vj, y=vdj, by='barcode', all=TRUE,
                             sort=FALSE, suffixes=c("_vj", "_vdj"))
    object <- count_chains(object,
                           vj=vj_options[object@celltype][[1]],
                           vdj=vdj_options[object@celltype][[1]])
    row.names(object@metadata) <- object@metadata$barcode
    object@contigs$keep <- TRUE

    ##### Ways to handle calls
    # keep_all returns object with all contigs unchanged
    # Removes all genes except highest one
    if(handle_multiple == 'keep_high_umi'){
        object <- keep_high_umi(object)
    }

    # Expand each cell to one row per vj x vdj contig pairing, ranked by combined umis
    if(handle_multiple %in% c('expand_clones', 'expand_clones_keep_one_clonal')){
        object@metadata <- expand_cell_clones(object@metadata) %>%
            dplyr::group_by(.data$barcode) %>%
            dplyr::mutate(umi_rank = rank(-.data$sum_umis, ties.method = 'first')) %>%
            dplyr::ungroup()

        # Count cells sharing each paired ntCDR3 and keep, per cell, the pairing seen in the
        # most cells (ties broken by umis). Only contigs of kept pairings are flagged keep.
        if(handle_multiple == 'expand_clones_keep_one_clonal'){
            object@metadata <- object@metadata %>%
                dplyr::group_by(.data$junction_vj, .data$junction_vdj) %>%
                dplyr::mutate(clone_count = dplyr::n_distinct(.data$barcode)) %>%
                dplyr::group_by(.data$barcode) %>%
                dplyr::mutate(keep = dplyr::row_number() == order(-.data$clone_count, .data$umi_rank)[1]) %>%
                dplyr::ungroup()
            kept <- object@metadata[object@metadata$keep, ]
            object@contigs$keep <- object@contigs[[contig_id_col]] %in%
                c(kept$contig_id_vj, kept$contig_id_vdj)
        }

        object@metadata <- as.data.frame(object@metadata)
        row.names(object@metadata) <- paste0(object@metadata$barcode, '_', object@metadata$umi_rank)
    }

    return(object)
})
