#' @title Load transcript to biotype
#'
#' @description Create a file containing the mapping between transcript IDs, type 
#' (genic or intergenic) and biotypes. Load data from this file.
#' This transcript to biotype mapping is used when expression estimation is not 
#' summurize at gene level but stay at transcript level. Will create the mapping 
#' file if not already existing.
#' 
#' @param myAbundanceMetadata A descendant object of the Class myAbundanceMetadata.
#' @param myBgeeMetadata A Reference Class BgeeMetadata object.
#' @param myUserMetadata A Reference Class UserMetadata object. 
#' This object has to be edited before running kallisto 
#' 
#' @seealso UserMetadata.R BgeeMetadata.R AbundanceMetadata.R
#'
#' @author Julien Wollbrett.
#'
#' @return Mapping between transcript IDs, type (genic or intergenic) and biotypes
#'
#' @import rtracklayer
#'
#' @noMd
#' @noRd
#'
load_transcript_to_biotype <- function(myAbundanceMetadata, 
    myBgeeMetadata, myUserMetadata) {
    column_names <- c("id", "biotype", "type")
    annotation_path <- get_annotation_path(myBgeeMetadata, 
        myUserMetadata)
    transcript_to_biotype_file <- file.path(annotation_path, 
        myAbundanceMetadata@tx2biotype_file)
    
    # check if file already exist
    if (!file.exists(transcript_to_biotype_file)) {
        if (!dir.exists(annotation_path)) {
            dir.create(annotation_path, recursive = TRUE)
        }
        cat(paste0("Generate file ", myAbundanceMetadata@tx2biotype_file, 
            ".\n"))
        # retrieve tx2biotype data frame from annotation
        # file
        gtf = as.data.frame(myUserMetadata@annotation_object)
        gtf_transcript <- gtf[gtf$type == "transcript", 
            ]
        transcript_to_biotype <- as.data.frame(unique(cbind(gtf_transcript$transcript_id, 
            gtf_transcript$transcript_biotype)))
        transcript_to_biotype[, 3] <- "genic"
        names(transcript_to_biotype) <- column_names
        
        # retrieve gene2biotype information from intergenic
        # fasta file
        bgee_intergenic_file <- file.path(get_species_path(myBgeeMetadata, 
            myUserMetadata), myBgeeMetadata@fasta_intergenic_name)
        if (!file.exists(bgee_intergenic_file)) {
            download_fasta_intergenic(myBgeeMetadata, 
                myUserMetadata, bgee_intergenic_file)
        }
        bgee_intergenic <- readDNAStringSet(bgee_intergenic_file)
        # intergenic ID correspond to part of the header
        # before the first space character
        intergenic_to_biotype <- as.data.frame(sub("^([^ ]+).*", 
            "\\1", names(bgee_intergenic)))
        intergenic_to_biotype[, 2] <- NA
        intergenic_to_biotype[, 3] <- "intergenic"
        names(intergenic_to_biotype) <- column_names
        
        # merge both data frame and write file
        transcript_to_biotype <- rbind(transcript_to_biotype, 
            intergenic_to_biotype)
        write.table(transcript_to_biotype, transcript_to_biotype_file, 
            sep = "\t", row.names = FALSE, quote = FALSE)
    }
    return(read.table(transcript_to_biotype_file, header = TRUE))
}
