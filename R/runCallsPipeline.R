# regroup all function used to easily run the
# RNA-Seq calls presence/absence pipeline

#' @title generate present/absent calls
#'
#' @description Main function running the workflow that generates
#' present/absent calls from a file, a data.frame, or objects of the
#' classe UserMetadata (please choose only 1 out of the 3). 
#' This workflow is highly tunable by editing default values of the slots of
#' S4 objects. For more information on how to tune the workflow please have
#' a look at the vignette and the documentation of the classes
#' KallistoMetadata, AbundanceMetadata, UserMetadata and BgeeMetadata
#'
#'
#' @param abundanceMetadata A Class AbundanceMetadata object (optional)
#' allowing to tune your gene quantification abundance analyze
#' @param bgeeMetadata A Class BgeeMetadata object (optional)
#' allowing to choose the version of reference intergenic sequences
#' @param userMetadata A Class UserMetadata object (optional).
#' generate present/absent calls using slots of the UserMetadata class.
#' @param userDataFrame a data.frame comtaining all information to generate
#' present/absent calls. Each line of this data.frame will generate calls for
#' one RNA-Seq library. This data.frame must contains between 4 and 8 columns :
#' \itemize{
#'  \item{species_id : The ensembl species ID}
#'  \item{run_ids : (optional) allows to generate calls for a subpart of all runs
#'   of the library. must be a character or a list of characters}
#'  \item{reads_size (optional) the size of the reads of the library (Default = 51)
#' if the reads size is lower than 51 abundance quantification will be run from an index
#'  generated with a smaller kmer size}
#'  \item{rnaseq_lib_path : path to RNA-Seq library directory}
#'  \item{transcriptome_path : path to transcriptome file}
#'  \item{annotation_path : path to annotation file}
#'  \item{output_dir : (optional)root of the directory where results will be written}
#'  \item{custom_intergenic_path : (optional) To use if the "custom" reference intergenic
#'  release has been selected. Provide the path to the reference intergenic file}
#' }
#'  
#' @param userFile path to a tsv file containing between 4 and 8 columns. these columns are
#' the same than for userDataFrame (see above). a template of this file is
#' available at the root of the package and accessible with the command
#' system.file('userMetadataTemplate.tsv', package = 'BgeeCall')
#' 
#' @param checkTxVersion boolean used to define if BgeeCall check rather transcript version
#' should be removed. Default value is FALSE
#'
#' @author Julien Wollbrett
#'
#' @return paths to the 5 results files (see vignette for more details)
#'
#' @export
#'
#' @seealso AbundanceMetadata, KallistoMetadata, BgeeMetadata, UserMetadata
#'
#' @examples
#' \dontrun{
#' # import gene annotation and transcriptome from AnnotationHub
#' library(AnnotationHub)
#' ah <- AnnotationHub()
#' ah_resources <- query(ah, c('Ensembl', 'Caenorhabditis elegans', '84'))
#' annotation_object <- ah_resources[['AH50789']]
#' transcriptome_object <- rtracklayer::import.2bit(ah_resources[['AH50453']])
#'
#' # instanciate BgeeCall object
#' # add annotation and transcriptome in the user_BgeeCall object
#' # it is possible to import them using an S4 object (GRanges, DNAStringSet)
#' # or a file (gtf, fasta) with methods setAnnotationFromFile() and
#' # setTranscriptomeFromFile()
#' user_BgeeCall <- setAnnotationFromObject(user_BgeeCall,
#'                                          annotation_object,
#'                                          'WBcel235_84')
#' user_BgeeCall <- setTranscriptomeFromObject(user_BgeeCall,
#'                                           transcriptome_object,
#'                                           'WBcel235')
#' # provide path to the directory of your RNA-Seq library
#' user_BgeeCall <- setRNASeqLibPath(user_BgeeCall,
#'                  system.file('extdata', 'SRX099901_subset',
#'                  package = 'BgeeCall'))
#'
#' # run the full BgeeCall workflow
#' calls_output <- generate_calls_workflow(
#'              userMetadata = user_BgeeCall)
#' }
#'
generate_calls_workflow <- function(abundanceMetadata = new("KallistoMetadata"),
                                    bgeeMetadata = new("BgeeMetadata"),
                                    userMetadata = NULL,
                                    userDataFrame = NULL,
                                    userFile = NULL,
                                    checkTxVersion = FALSE) {
    if (is.null(userMetadata) && is.null(userDataFrame) &&
        is.null(userFile)) {
        stop("one of the parameters userMetadata, userDataFrame or userFile
            sould not be NULL")
        
        # run workflow when userMetadata is not null
    } else if (!is.null(userMetadata) && is.null(userDataFrame) &&
        is.null(userFile)) {
        if (isS4(userMetadata)) {
            return(
                run_from_object(
                    myAbundanceMetadata = abundanceMetadata,
                    myBgeeMetadata = bgeeMetadata,
                    myUserMetadata = userMetadata,
                    checkTxVersion = checkTxVersion
                )
            )
        } else if (typeof(userMetadata) == "list" &&
            isS4(userMetadata[[1]])) {
            for (i in seq_along(userMetadata)) {
                results[i] <-
                    run_from_object(
                        myAbundanceMetadata = abundanceMetadata,
                        myBgeeMetadata = bgeeMetadata,
                        myUserMetadata = userMetadata[[i]],
                        checkTxVersion = checkTxVersion
                    )
            }
            return(results)
        } else {
            stop(
                "the parameter userMetadata should only be used to
                provide one or a list of UserMetadata objects"
            )
        }
        # run workflow when userDataFrame is not null
    } else if (is.null(userMetadata) && !is.null(userDataFrame) &&
        is.null(userFile)) {
        return(
            run_from_dataframe(
                myAbundanceMetadata = abundanceMetadata,
                myBgeeMetadata = bgeeMetadata,
                myUserMetadata = myUserMetadata,
                userMetadataDataFrame = userDataFrame,
                checkTxVersion = checkTxVersion
            )
        )
        # run workflow when userDataFrame is not null
    } else if (is.null(userDataFrame) && !is.null(userFile)) {
        if ((!(typeof(userFile) == "character") &&
            file.exists(userFile))) {
            stop(
"Please provide a path to the file that contains all
information allowing to generate UserMetadata objects"
            )
        }
        return(
            run_from_file(
                myAbundanceMetadata = abundanceMetadata,
                myBgeeMetadata = bgeeMetadata,
                myUserMetadata = userMetadata,
                userMetadataFile = userFile,
                checkTxVersion = checkTxVersion
            )
        )
    } else {
        stop(
            "Please use only 1 of the 3 follwowing attributs :
            userMetadata, userDataFrame, userFile"
        )
    }
}

#' @title generate present/absent calls from a UserMetadata object
#'
#' @description function allowing to generate present/absent calls for one
#' library described as a UserMetadata object
#'
#' @param myAbundanceMetadata A Reference Class BgeeMetadata object (optional)
#' @param myBgeeMetadata A Reference Class BgeeMetadata object (optional)
#' @param myUserMetadata A Reference Class UserMetadata object.
#'
#' @author Julien Wollbrett
#'
#' @return 4 paths to results files
#'
#' @seealso generate_calls_workflow
#'
#' @examples {
#' # import annotation and transcriptome from AnnotationHub
#' library(AnnotationHub)
#' ah <- AnnotationHub()
#' ah_resources <- query(ah, c('Ensembl', 'Caenorhabditis elegans', '84'))
#' annotation_object <- ah_resources[['AH50789']]
#' transcriptome_object <- rtracklayer::import.2bit(ah_resources[['AH50453']])
#'
#' # instanciate BgeeCall object
#' # add annotation and transcriptome in the user_BgeeCall object
#' # it is possible to import them using an S4 object (GRanges, DNAStringSet)
#' # or a file (gtf, fasta) with methods setAnnotationFromFile() and
#' # setTranscriptomeFromFile()
#' user_BgeeCall <- setAnnotationFromObject(user_BgeeCall,
#'                                          annotation_object,
#'                                          'WBcel235_84')
#' user_BgeeCall <- setTranscriptomeFromObject(user_BgeeCall,
#'                                           transcriptome_object,
#'                                           'WBcel235')
#' # provide path to the directory of your RNA-Seq library
#' user_BgeeCall <- setRNASeqLibPath(user_BgeeCall,
#'                  system.file('extdata',
#'                  'SRX099901_subset',
#'                  package = 'BgeeCall'))
#'
#' # run the full BgeeCal workflow
#' calls_output <- run_from_object(
#'              myUserMetadata = user_BgeeCall)
#' }
#'
#' @noMd
#' @noRd
#'

run_from_object <-
    function(myAbundanceMetadata = new("KallistoMetadata"),
        myBgeeMetadata = new("BgeeMetadata"),
        myUserMetadata, checkTxVersion) {
      if(checkTxVersion) {
        myAbundanceMetadata@ignoreTxVersion <- should_ignore_tx_version(myUserMetadata)
      }
      if (myAbundanceMetadata@tool_name == "kallisto") {
          run_kallisto(myAbundanceMetadata, myBgeeMetadata,
              myUserMetadata)
      } else {
          stop(
              paste0(
                  "The myAbundanceMetadata object should be an instance of
              KallistoMetadata"
              )
          )
      }
      calls_output <- generate_presence_absence(myAbundanceMetadata,
          myBgeeMetadata, myUserMetadata)
      return(calls_output)
    }

#' @title generate present/absent calls from a data frame
#'
#' @description Function allowing to generate present/absent calls for
#' libraries described in a data frame
#'
#' @param myAbundanceMetadata A Reference Class BgeeMetadata object (optional)
#' @param myBgeeMetadata A Reference Class BgeeMetadata object (optional)
#' @param myUserMetadata A class BgeeMetadata object (optional). Used to modify
#' values slots that are not part of the userMetadataFile
#' @param userMetadataDataFrame A data frame containing all information needed
#' to generate UserMetadata objects. This data frame must contains 7 columns
#' (species_id, run_ids, reads_size, rnaseq_lib_path, transcriptome_path,
#' annotation_path, working_path)
#' @param checkTxVersion boolean used to define if BgeeCall check rather transcript version
#' should be removed. Default value is FALSE
#' 
#' @author Julien Wollbrett
#'
#' @return paths to the 4 output files generated per call generation
#'
#' @seealso run_from_file()
#'
#' @examples {
#' # your data.frame containing 7 columns
#' user_metadata_df
#' run_from_dataframe(userMetadataDataFrame = metadata_file)
#' }
#'
#' @noMd
#' @noRd
#'
run_from_dataframe <-function(myAbundanceMetadata = new("KallistoMetadata"),
        myBgeeMetadata = new("BgeeMetadata"), myUserMetadata = NULL, 
        userMetadataDataFrame, checkTxVersion) {
    all_results <- list()
    for (row_number in seq_len(nrow(userMetadataDataFrame))) {
        ## init myUserMetadata object
        user_metadata <-
            init_userMetadata_from_dataframe(userMetadataDataFrame, myUserMetadata, row_number)
        
        # run pipeline
        result <- run_from_object(myAbundanceMetadata = myAbundanceMetadata, 
            myBgeeMetadata = myBgeeMetadata, myUserMetadata = user_metadata, 
            checkTxVersion = checkTxVersion)
        all_results <- c(all_results, list(result))
    }
    return(all_results)
}

#' @title generate present/absent calls from a file
#'
#' @description function allowing to generate present/absent calls for some
#' libraries described in a file. Each line of the file describes one RNA-Seq
#' library.
#'
#' @param myAbundanceMetadata A Reference Class BgeeMetadata object (optional)
#' @param myBgeeMetadata A Reference Class BgeeMetadata object (optional)
#' @param myUserMetadata A class BgeeMetadata object (optional). Used to modify
#' values slots that are not part of the userMetadataFile
#' @param userMetadataFile A tsv file describing all user libraries for which
#' present/absent calls have to be generated. A template of this file named
#' `userMetadataTemplate.tsv` is available at the root of the `BgeeCall`
#' package. It is a tabular separated value file containing 7 columns :
#' - species_id : species ID
#' - run_ids : (optional) only if you want to generate calls for a subpart
#' of all runs of the library
#' - reads_size (optional) the size of the reads of the library (Default = 50)
#' - rnaseq_lib_path : path to RNA-Seq library directory
#' - transcriptome_path : path to transcriptome file
#' - annotation_path : path to annotation file
#' - working_path : root of the output directory
#' @param checkTxVersion boolean used to define if BgeeCall check rather transcript version
#' should be removed. Default value is FALSE
#' 
#' @author Julien Wollbrett
#'
#' @return paths to the 4 output files generated per call generation
#'
#' @examples {
#' metadata_file <- system.file('path/to/your/file.tsv', package = 'BgeeCall')
#' run_from_file(userMetadataFile = metadata_file)
#' }
#'
#' @noMd
#' @noRd
#'
run_from_file <-
    function(myAbundanceMetadata = new("KallistoMetadata"), myBgeeMetadata = new("BgeeMetadata"), 
             myUserMetadata = NULL, userMetadataFile, checkTxVersion) {
    user_metadata_df <- read.table(userMetadataFile, header = TRUE, sep = "\t", comment.char = "#")
    output_files <- run_from_dataframe(myAbundanceMetadata, myBgeeMetadata, myUserMetadata, 
                                       user_metadata_df, checkTxVersion)
    return(output_files)
}

init_userMetadata_from_dataframe <- function(userMetadataDataFrame, myUserMetadata = NULL, 
                                             row_number) {
    if (is.null(myUserMetadata)) {
        myUserMetadata <- new("UserMetadata")
    }
    myUserMetadata@species_id <-as.character(userMetadataDataFrame[["species_id"]][row_number])
    
    # check if subset of run ids has to be used to generate present/absent
    myUserMetadata@run_ids <- check_run_ids(as.character(userMetadataDataFrame[["run_ids"]][row_number]))
    # check if user provided an output_dir or if the default one will be used
    if ("output_directory" %in% names(userMetadataDataFrame)) {
        output_dir <- as.character(userMetadataDataFrame[["output_directory"]][row_number])
        if (length(output_dir) != 0) {
            if (!dir.exists(output_dir)) {
                dir.create(path = output_dir, recursive = TRUE, showWarnings = TRUE)
            }
            myUserMetadata <- setOutputDir(myUserMetadata, output_dir)
        }
    }
    
    # check if user provided an output_dir or if the default one will be used
    if ("custom_intergenic_path" %in% names(userMetadataDataFrame)) {
        custom_intergenic_path <- as.character(
            userMetadataDataFrame[["custom_intergenic_path"]][row_number])
        if (length(custom_intergenic_path) != 0) {
            myUserMetadata@custom_intergenic_path <-custom_intergenic_path
        }
    }
    
    myUserMetadata@reads_size <-  as.numeric(userMetadataDataFrame[["reads_size"]][row_number])
    myUserMetadata@rnaseq_lib_path <- as.character(
        userMetadataDataFrame[["rnaseq_lib_path"]][row_number])
    myUserMetadata <- setTranscriptomeFromFile(userObject = myUserMetadata, 
                    transcriptomePath = as.character(
                        userMetadataDataFrame[["transcriptome_path"]][row_number]))
    myUserMetadata <- setAnnotationFromFile(userObject = myUserMetadata, 
                    annotationPath = as.character(
                        userMetadataDataFrame[["annotation_path"]][row_number]))
    if (is.na(myUserMetadata@working_path) || myUserMetadata@working_path == '') {
        myUserMetadata@working_path <- as.character(
            userMetadataDataFrame[["working_path"]][row_number])
    }
    return(myUserMetadata)
}
