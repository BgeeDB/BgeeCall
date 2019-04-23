#' @title Merge transcriptome file provided by the user with the Bgee intergenic 
#' fasta file.
#'
#' @description This function will create a file corresponding to the concatenation 
#' of the transcriptome fasta file provided by the user and the corresponding 
#' intergenic fasta file created by Bgee.
#'
#' @param myKallistoMetadata A Reference Class KallistoMetadata object.
#' @param myBgeeMetadata A Reference Class BgeeMetadata object.
#' @param myUserMetadata A Reference Class UserMetadata object.
#'
#' @author Julien Wollbrett.
#'
#' @import Biostrings
#'
#' @export
#' 
#' 
#' @examples {
#' bgee <- new('BgeeMetadata', intergenic_release = '0.1')
#' user <- new ('UserMetadata', species_id = '6239')
#' kallisto <- new('KallistoMetadata')
#' user <- setTranscriptomeFromFile(user, system.file("extdata", 
#' "transcriptome.fa", package = "BgeeCall"), 'WBcel235')
#' merge_transcriptome_and_intergenic(kallisto, bgee, user)
#' }
#' 
#' @return save merged file on the hard drive

merge_transcriptome_and_intergenic <- function(myKallistoMetadata, 
    myBgeeMetadata, myUserMetadata) {
    
    # define needed path
    species_dir <- get_species_path(myBgeeMetadata, 
        myUserMetadata)
    bgee_intergenic_file <- file.path(species_dir, 
        myBgeeMetadata@fasta_intergenic_name)
    transcriptome_dir <- get_transcriptome_path(myBgeeMetadata, 
        myUserMetadata)
    transcriptome_with_intergenic_path <- file.path(transcriptome_dir, 
        myKallistoMetadata@full_transcriptome_file)
    
    # test if file containing both transcriptomic and
    # intergenic regions alredy exists
    if (file.exists(transcriptome_with_intergenic_path)) {
        message("File containing both transcriptomic and intergenic 
                regions already exists.\n")
    } else {
        message("Start generation of the file containing both transcriptomic
                and intergenic regions.\n")
        if (!dir.exists(transcriptome_dir)) {
            dir.create(transcriptome_dir, recursive = TRUE)
        }
        
        # read Bgee fasta intergenic file
        if (!file.exists(bgee_intergenic_file)) {
            download_fasta_intergenic(myBgeeMetadata, 
                myUserMetadata, bgee_intergenic_file)
        }
        bgee_intergenic <- readDNAStringSet(bgee_intergenic_file)
        
        # combine transcriptome and intergenic fasta files
        transcriptome_with_intergenic_data <- c(myUserMetadata@transcriptome_object, 
            bgee_intergenic)
        writeXStringSet(x = transcriptome_with_intergenic_data, 
            filepath = transcriptome_with_intergenic_path)
        message("File containing both transcriptomic and intergenic regions has
                been created successfully.\n")
    }
}
