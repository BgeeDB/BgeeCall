# regroup all function used to easily run the RNA-Seq calls presence/absence pipeline

#' @title generate present/absent calls from a UserMetadata object
#'
#' @description Main function allowing to generate present/absent calls for one library described as a UserMetadata object
#'
#' @param myAbundanceMetadata A Reference Class BgeeMetadata object (optional)
#' @param myBgeeMetadata A Reference Class BgeeMetadata object (optional)
#' @param myUserMetadata A Reference Class UserMetadata object.
#'
#' @author Julien Wollbrett
#' 
#' @export
#'
run_from_object <- function (myAbundanceMetadata = new("KallistoMetadata"), myBgeeMetadata = new("BgeeMetadata"), myUserMetadata) {
  #if(!file.exists(file.path(get_tool_output_path(myAbundanceMetadata, myBgeeMetadata, myUserMetadata), myAbundanceMetadata@abundance_file))) {
    if (myAbundanceMetadata@tool_name == "kallisto") {
      run_kallisto(myAbundanceMetadata, myBgeeMetadata, myUserMetadata)
    } else {
      stop(paste0("The myAbundanceMetadata object should be an instance of KallistoMetadata"))
    }
    calls_output <- generate_presence_absence(myAbundanceMetadata, myBgeeMetadata, myUserMetadata)
    return(calls_output)
}

#' @title generate present/absent calls from a data frame
#'
#' @description Main function allowing to generate present/absent calls for one library described in a data frame
#'
#' @param myAbundanceMetadata A Reference Class BgeeMetadata object (optional)
#' @param myBgeeMetadata A Reference Class BgeeMetadata object (optional)
#' @param myUserMetadata A data frame containing all information needed to generate UserMetadata objects
#'
#' @author Julien Wollbrett
#' 
#' @export
#'
run_from_dataframe <- function (myAbundanceMetadata = new("KallistoMetadata"), myBgeeMetadata = new("BgeeMetadata"), userMetadataDataFrame) {
  for (row_number in 1:nrow(userMetadataDataFrame)) {
    
    # init myUserMetadata object
    myUserMetadata <- new("UserMetadata")
    myUserMetadata@species_id <- as.character(userMetadataDataFrame[["species_id"]][row_number])
    ids <- as.character(userMetadataDataFrame[["run_id"]][row_number])
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
    myUserMetadata@reads_size <- as.numeric(userMetadataDataFrame[["reads_size"]][row_number])
    myUserMetadata@rnaseq_lib_path <- as.character(userMetadataDataFrame[["rnaseq_lib_path"]][row_number])
    myUserMetadata <- setTranscriptomeFromFile(userObject = myUserMetadata, 
                                               transcriptomePath = as.character(userMetadataDataFrame[["transcriptome_path"]][row_number]))
    myUserMetadata <- setAnnotationFromFile(userObject = myUserMetadata, 
                                            annotationPath = as.character(userMetadataDataFrame[["annotation_path"]][row_number]))
    myUserMetadata@working_path <- as.character(userMetadataDataFrame[["working_path"]][row_number])

    # run pipeline
    run_from_object(myAbundanceMetadata, myBgeeMetadata, myUserMetadata)
  }
}

#' @title generate present/absent calls from a file
#'
#' @description Main function allowing to generate present/absent calls for some 
#' libraries described in a file. Each line of the file describes one RNA-Seq library.
#'
#' @param myAbundanceMetadata A Reference Class BgeeMetadata object (optional)
#' @param myBgeeMetadata A Reference Class BgeeMetadata object (optional)
#' @param userMetadataFile A tsv file describing all user libraries for which 
#' present/absent calls have to be generated
#'
#' @author Julien Wollbrett
#' 
#' @export
#'
run_from_file <- function (myAbundanceMetadata = new("KallistoMetadata"), 
                           myBgeeMetadata = new("BgeeMetadata"), userMetadataFile) {
  user_metadata_df <- read.table(userMetadataFile, header = TRUE, sep = "\t",
                                 comment.char = "#")
  run_from_dataframe(myAbundanceMetadata, myBgeeMetadata, user_metadata_df)
}

