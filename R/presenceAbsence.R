#' @title Retrieve dataframe of reference intergenic IDs
#'
#' @description Retrieve a dataframe of reference intergenic IDs. It is generated using the reference intergenic fasta file from Bgee FTP.
#'
#' @param myBgeeMetadata A Class BgeeMetadata object.
#' @param myUserMetadata A Class UserMetadata object.
#'
#' @author Julien Wollbrett
#' @author Julien Roux
#'
#' @return A dataframe containing reference intergenic ids
#' 
#' @noMd
#'
get_ref_intergenic_ids <- function(myBgeeMetadata, myUserMetadata) {
  bgee_intergenic_file <- file.path(get_species_path(myBgeeMetadata, myUserMetadata), myBgeeMetadata@fasta_intergenic_name)
  if (!file.exists(bgee_intergenic_file)) {
    # Download fasta file from Bgee FTP
    download_fasta_intergenic(myBgeeMetadata, myUserMetadata, bgee_intergenic_file)
  }
  bgee_intergenic <- readDNAStringSet(bgee_intergenic_file)
  #keep only intergenic ids from fasta file
  return(as.data.frame(sub("^([^ ]+).*", "\\1", names(bgee_intergenic))))
}

#' @title Generate presence absence
#'
#' @description Generate presence absence calls
#'
#' @param myAbundanceMetadata A descendant object of the Class myAbundanceMetadata.
#' @param myBgeeMetadata A Class BgeeMetadata object.
#' @param myUserMetadata A Class UserMetadata object.
#'
#' @author Julien Wollbrett
#' @author Julien Roux
#'
#' @export
#'
generate_presence_absence <- function(myAbundanceMetadata, myBgeeMetadata, myUserMetadata) {

  # load data
  ref_intergenic <- get_ref_intergenic_ids(myBgeeMetadata, myUserMetadata)
  tool_path <- get_tool_path(myAbundanceMetadata, myBgeeMetadata, myUserMetadata)
  output_path <- get_tool_output_path(myAbundanceMetadata, myBgeeMetadata, myUserMetadata)
  
  # biotype mapping information will depend on summarization at gene level or not 
  biotype_mapping <- ""
  if (myAbundanceMetadata@txOut) {
    biotype_mapping <- load_transcript_to_biotype(myAbundanceMetadata, myBgeeMetadata, myUserMetadata)
  } else {
    biotype_mapping <- load_gene_to_biotype(myAbundanceMetadata, myBgeeMetadata, myUserMetadata)
  }
  
  #run tximport for file with intergenic regions (if myAbundanceMetadata@txOut = FALSE, then tximport will summurarize transcript level estimates at gene level)
  tximportObject <- run_tximport(myAbundanceMetadata, myBgeeMetadata, myUserMetadata)

  # recalculate TPM without intergenic regions and run tximport (if myAbundanceMetadata@txOut = FALSE, then tximport will summurarize transcript level estimates at gene level)
  tximportObject_without_intergenic <- abundance_without_intergenic(myAbundanceMetadata, myBgeeMetadata, myUserMetadata)

  # transform tximport output in order to easily process information
  abundance <- transform_tximport(tximportObject, biotype_mapping)
  abundance_without_intergenic <- transform_tximport(tximportObject_without_intergenic, biotype_mapping)

  # define coding and intergenic abundance subset
  selected_coding <- abundance$biotype %in% "protein_coding"
  selected_intergenic <- (abundance$type == "intergenic" & abundance$id %in% biotype_mapping$id)

  #calculate TPM cutoff
  results <- calculate_abundance_cutoff(abundance, selected_coding, selected_intergenic, myAbundanceMetadata@cutoff)
  abundance_cutoff <- results[1]
  r_cutoff <- results[2]
  
  # generate name of output file
  calls_file_name <- myAbundanceMetadata@gene_calls_file_name
  cutoff_info_file_name <- myAbundanceMetadata@gene_cutoff_file_name
  distribution_file_name <- myAbundanceMetadata@gene_distribution_file_name
  if (myAbundanceMetadata@txOut) {
    calls_file_name <- myAbundanceMetadata@transcript_calls_file_name
    cutoff_info_file_name <- myAbundanceMetadata@transcript_cutoff_file_name
    distribution_file_name <- myAbundanceMetadata@transcript_distribution_file_name
  }

  # generate pdf plot
  pdf(file = file.path(output_path, distribution_file_name), width = 6, height = 5)
  plot_distributions(abundance, selected_coding, selected_intergenic, abundance_cutoff, myUserMetadata)

  # add presence/absence information
  abundance$call <- ifelse(abundance$abundance >= abundance_cutoff, "present", "absent")

  # generate cutoff info file
  cutoff_info_file <- cutoff_info(abundance, "call", abundance_cutoff, r_cutoff, myUserMetadata)
  dev.off()

  # transfert presence/absence annotations to the abundances without intergenic
  abundance_without_intergenic <- merge(abundance_without_intergenic, abundance[, c("id", "call")], by="id")
  
  # Save calls and stats to output folder
  calls_file_path <- file.path(output_path, calls_file_name)
  cutoff_info_file_path <- file.path(output_path, cutoff_info_file_name)
  write.table(abundance_without_intergenic,
              file = calls_file_path,
              quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
  write.table(t(t(cutoff_info_file)),
              file = cutoff_info_file_path,
              quote = FALSE, sep = "\t", col.names = FALSE, row.names = TRUE)
  calls_result <- list()
  calls_result$calls_tsv_path <- calls_file_path
  calls_result$cutoff_info_file_path <- cutoff_info_file_path
  calls_result$abundance_tsv <- file.path(output_path, myAbundanceMetadata@abundance_file)
  calls_result$TPM_distribution_path <- file.path(output_path, distribution_file_name)
  return(calls_result)
  ## t(t(cutoff_info_file)) is a solution to export a vector vertically
}

#' @title Transform txinport object
#'
#' @description transform tximport object in order to easily process information
#'
#' @param tximportObject A tximport Object
#' @param biotype_mapping Mapping between gene or transcript IDs and biotypes
#'
#' @return A tximport Object with biotype and without the countsFromAbundance column
#'
#' @noMd
#' @noRd
#'
transform_tximport <- function (tximportObject, biotype_mapping) {
  tx_df <- as.data.frame(tximportObject)
  tx_df$id <- rownames(tx_df)
  tx_df$countsFromAbundance <- NULL
  abundance <- merge(tx_df, biotype_mapping, by = "id", all = FALSE)
}

