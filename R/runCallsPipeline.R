# regroup all function used to easily run the RNA-Seq calls presence/absence 
# pipeline

#' @title generate present/absent calls from a UserMetadata object
#'
#' @description Main function allowing to generate present/absent calls for one 
#' library described as a UserMetadata object
#'
#' @param myAbundanceMetadata A Reference Class BgeeMetadata object (optional)
#' @param myBgeeMetadata A Reference Class BgeeMetadata object (optional)
#' @param myUserMetadata A Reference Class UserMetadata object.
#'
#' @author Julien Wollbrett
#' 
#' @export
#' 
#' @return 4 paths to results files
#' 
#' @examples {
#' library(AnnotationHub)
#' ah <- AnnotationHub()
#' ah_resources <- query(ah, c("Ensembl", "Caenorhabditis elegans", "84"))
#' annotation_object <- ah_resources[["AH50789"]]
#' transcriptome_object <- rtracklayer::import.2bit(ah_resources[["AH50453"]])
#' user_BgeeCall <- new("UserMetadata", species_id = "6239")
#' # import annotation and transcriptome in the user_BgeeCall object
#' # it is possible to import them using an S4 object (GRanges, DNAStringSet) 
#' # or a file (gtf, fasta) with methods setAnnotationFromFile() and
#' # setTranscriptomeFromFile()
#' user_BgeeCall <- setAnnotationFromObject(user_BgeeCall, 
#'                                          annotation_object, 
#'                                          "WBcel235_84")
#' user_BgeeCall <- setTranscriptomeFromObject(user_BgeeCall, 
#'                                           transcriptome_object, 
#'                                           "WBcel235")
#' # provide path to the directory of your RNA-Seq library
#' user_BgeeCall <- setRNASeqLibPath(user_BgeeCall, 
#'                  system.file("extdata", "SRX099901_subset", package = "BgeeCall"))
#' calls_output <- run_from_object(
#'              myUserMetadata = user_BgeeCall)                                           
#' }
#'

run_from_object <- function (myAbundanceMetadata = new("KallistoMetadata"), 
                             myBgeeMetadata = new("BgeeMetadata"), 
                             myUserMetadata) {
    if (myAbundanceMetadata@tool_name == "kallisto") {
        run_kallisto(myAbundanceMetadata, myBgeeMetadata, myUserMetadata)
    } else {
        stop(paste0("The myAbundanceMetadata object should be an instance of 
                KallistoMetadata"))
    }
    calls_output <- generate_presence_absence(myAbundanceMetadata, myBgeeMetadata, 
                                              myUserMetadata)
    return(calls_output)
}

#' @title generate present/absent calls from a data frame
#'
#' @description Main function allowing to generate present/absent calls for one 
#' library described in a data frame
#'
#' @param myAbundanceMetadata A Reference Class BgeeMetadata object (optional)
#' @param myBgeeMetadata A Reference Class BgeeMetadata object (optional)
#' @param userMetadataDataFrame A data frame containing all information needed to 
#' generate UserMetadata objects. This data frame must contains 7 columns (species_id, 
#' run_ids, reads_size, rnaseq_lib_path, transcriptome_path, annotation_path, working_path)
#'
#' @author Julien Wollbrett
#' 
#' @export
#' 
#' @return paths to the 4 output files generated per call generation
#' 
#' @seealso run_from_file()
#' 
#' @examples \dontrun{
#' metadata_file <- system.file("userMetadataTemplate.tsv",
#'                                package = "BgeeCall")
#' user_metadata_df <- read.table(metadata_file, header = TRUE, sep = "\t",
#'                               comment.char = "#")
#' run_from_dataframe(userMetadataDataFrame = metadata_file)
#' }
#'
run_from_dataframe <- function (myAbundanceMetadata = new("KallistoMetadata"), 
                                myBgeeMetadata = new("BgeeMetadata"), 
                                userMetadataDataFrame) {
    
    for (row_number in seq_len(userMetadataDataFrame)) {
        
        # init myUserMetadata object
        myUserMetadata <- new("UserMetadata")
        myUserMetadata@species_id <- as.character(
            userMetadataDataFrame[["species_id"]][row_number])
        
        ids <- as.character(userMetadataDataFrame[["run_ids"]][row_number])
        if (length(ids) == 0) {
            myUserMetadata@run_ids <- character(0)
        } else if (length(ids) == 1) {
            if (grepl(", ", ids)) {
                myUserMetadata@run_ids <-strsplit(ids, ", ")
            } else if (grepl(pattern = ",", x = ids)) {
                myUserMetadata@run_ids <-strsplit(ids, ",")
            }
        } else {
            myUserMetadata@run_ids <- ids
        }
        myUserMetadata@reads_size <- as.numeric(
            userMetadataDataFrame[["reads_size"]][row_number])
        myUserMetadata@rnaseq_lib_path <- as.character(
            userMetadataDataFrame[["rnaseq_lib_path"]][row_number])
        myUserMetadata <- setTranscriptomeFromFile(
            userObject = myUserMetadata, 
            transcriptomePath = as.character(
                userMetadataDataFrame[["transcriptome_path"]][row_number]),
            transcriptomeName = "")
        myUserMetadata <- setAnnotationFromFile(
            userObject = myUserMetadata, annotationPath = 
                as.character(userMetadataDataFrame[["annotation_path"]][row_number]),
            annotationName = "")
        myUserMetadata@working_path <- as.character(
            userMetadataDataFrame[["working_path"]][row_number])
        
        # run pipeline
        run_from_object(myAbundanceMetadata, myBgeeMetadata, myUserMetadata)
    }
}

#' @title generate present/absent calls from a file
#'
#' @description Main function allowing to generate present/absent calls for some 
#' libraries described in a file. Each line of the file describes one RNA-Seq 
#' library.
#'
#' @param myAbundanceMetadata A Reference Class BgeeMetadata object (optional)
#' @param myBgeeMetadata A Reference Class BgeeMetadata object (optional)
#' @param userMetadataFile A tsv file describing all user libraries for which 
#' present/absent calls have to be generated. A template of this file named 
#' `userMetadataTemplate.tsv` is available at the root of the `BgeeCall` package.
#' It is a tabular separated value file containing 7 columns :
#' - species_id : species ID
#' - run_ids : (optional) only if you want to generate calls for a subpart of all runs 
#'   of the library 
#' - reads_size (optional) the size of the reads of the library (Default = 50)
#' - rnaseq_lib_path : path to RNA-Seq library directory
#' - transcriptome_path : path to transcriptome file
#' - annotation_path : path to annotation file
#' - working_path : root of the output directory
#'
#' @author Julien Wollbrett
#' 
#' @export
#' 
#' @return paths to the 4 output files generated per call generation
#' 
#' @examples \dontrun{
#' metadata_file <- system.file("userMetadataTemplate.tsv", package = "BgeeCall")
#' run_from_file(userMetadataFile = metadata_file)
#' }
#'
run_from_file <- function (myAbundanceMetadata = new("KallistoMetadata"), 
                           myBgeeMetadata = new("BgeeMetadata"), 
                           userMetadataFile) {
    user_metadata_df <- read.table(userMetadataFile, header = TRUE, sep = "\t",
                                   comment.char = "#")
    output_files <- run_from_dataframe(myAbundanceMetadata, myBgeeMetadata, user_metadata_df)
    return(output_files)
}

