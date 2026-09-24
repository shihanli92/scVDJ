setGeneric("readVDJ_10X", function(object, path) standardGeneric("readVDJ_10X"))

#' Parses VDJ alignments and adds to contigs slot
#'
#' @param object: VDJ object which contains a path string indicating location of 10x VDJ output files.
#' @importFrom utils read.csv
#' @return VDJ object
setMethod("readVDJ_10X", "VDJ",

          function(object, path) {
              # Check that input object is a VDJ object
              if(class(object)[[1]] != 'VDJ'){
                  stop("Input must be a VDJ object")
              }
              # Check if input is a directory, remove trailing '\' and add file names to search
              to_process <- list()
              if(dir.exists(path)){
                  path <- stringr::str_replace(path,'/$', '')

                  if(file.exists(paste0(path, '/filtered_contig_annotations.csv', sep=''))){
                      contigs <- read.csv(paste0(path, '/filtered_contig_annotations.csv', sep=''))
                      to_process <- append(to_process, list(contigs))
                  }
                  else{
                      stop('Could not find filtered_contig_annotations.csv file in directory')
                  }

                  # Check for airr rearrangement tsv file, read and append dataframe to contig dataframe if found
                  if(file.exists(paste0(path, '/airr_rearrangement.tsv', sep=''))){
                      airr <- read.csv(paste0(path, '/airr_rearrangement.tsv', sep=''), sep='\t')
                      airr <- airr[,c("sequence_id", 'sequence', 'sequence_aa', 'v_cigar',
                                      'd_cigar', 'j_cigar', 'c_cigar', 'sequence_alignment',
                                      'germline_alignment', 'junction_length', 'junction_aa_length',
                                      'v_sequence_start', 'v_sequence_end', 'd_sequence_start',
                                      'd_sequence_end', 'j_sequence_start', 'j_sequence_end',
                                      'c_sequence_start', 'c_sequence_end', 'consensus_count',
                                      'duplicate_count')]
                      to_process <- append(to_process, list(airr), 1)
                  }
              }
              # If inputted path is not a directory, check if file name is filtered contig annotations and read if it is
              else{
                  if(grepl('filtered_contig_annotations.csv', path)){
                      contigs <- read.csv(path)
                      to_process <- append(to_process, list(contigs))
                  }
                  else{
                      stop('Could not find filtered_contig_annotations.csv file in directory')
                  }
              }
              # If airr file found merge two files together
              if(length(to_process) > 1){
                  to_return <- merge(x = to_process[[1]], y=to_process[[2]],
                                     by.x='contig_id', by.y='sequence_id', sort=FALSE)
              }
              else{
                  to_return <- to_process[[1]]
              }
              # Rename files to AIRR nomenclature
              to_return <- to_return %>% dplyr::rename(tidyselect::any_of(c(v_call = 'v_gene', j_call='j_gene',
                                                                            d_call='d_gene', c_call='c_gene',
                                                                            junction='cdr3_nt', junction_aa = 'cdr3')))
              # Keep productive
              if('productive' %in% colnames(to_return)){
                  to_return <- to_return[to_return$productive %in% c('TRUE', 'True', 'true', TRUE),]
              }
              to_return <- to_return[to_return$chain != 'Multi',]
              object@contigs <- to_return

              object@contigs <- object@contigs %>%
                  dplyr::group_by(junction) %>%
                  dplyr::mutate(contig_count=dplyr::n()) %>%
                  dplyr::ungroup()

              return(object)})
