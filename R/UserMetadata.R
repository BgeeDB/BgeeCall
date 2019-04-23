#' @title UserMetadata S4 class
#' 
#' @description An S4 class containing all metadata that have to be provided by 
#' the user It is mandatory to edit `species_id`, `rnaseq_lib_path`, 
#' `transcriptome_path`, `annotation_name`, `annotation_object` and potentialy 
#' `run_ids` before using the package.
#'
#' @slot species_id The NCBI Taxon Id of the species
#' @slot run_ids A vector of charater. Has to be provided only if a subset of runs
#' present in UserMetadata@rnaseq_lib_path has to be run. If empty, all fastq 
#' files present in the rnaseq_lib_path will be considered as technical 
#' replicates and merged to run one transcript expression estimation analyse.
#' @slot reads_size The size of the reads. If smaller than 
#' `KallistoMetadata@read_size_kmer_threshold`, 
#' an index with a kmer size of 21 nt will be used.
#' @slot rnaseq_lib_path Path to the directory of the RNA-Seq library that 
#' contains fastq files. 
#' @slot transcriptome_name Name of the transcriptome used to generate 
#' arborescence of output repositories.
#' @slot transcriptome_object Object containing transcriptome 
#' @slot annotation_name Name of the annotation used to generate arborescence 
#' of output repositories. 
#' @slot annotation_object Object containing annotations from GTF or GFF file
#' @slot working_path Working directory. By default the working directory is 
#' defined with the `getwd()` function.
#' @slot simple_arborescence logical allowing to create a simple arborescence 
#' of directory. All library ids are on the same directory. Default value is 
#' `FALSE`. Do not use `TRUE` if you plan to generate expression calls for the 
#' same library using different transcriptomes or gene annotations, otherwise 
#' you will overwrite previous results

UserMetadata <- setClass(
    # Set the name for the class
    Class = "UserMetadata",
    
    # Define the slots
    representation = representation(
        species_id = "character",
        run_ids = "character",
        reads_size = "numeric",
        rnaseq_lib_path = "character",
        transcriptome_name = "character",
        transcriptome_object = "DNAStringSet",
        annotation_name = "character",
        annotation_object = "GRanges",
        working_path = "character",
        simple_arborescence = "logical"
    ),
    
    # Set the default values for the slots.
    prototype = prototype(
        working_path = getwd(),
        reads_size = 51,
        simple_arborescence = FALSE
    )
)

#' 
#' @title Set annotation_object of one UserMetadata object
#' 
#' @description Method of the class UserMetadata. Set annotation_object of one 
#' UserMetadata object by using one GRanges object as input.
#' @param userObject The UserMetadata object
#' @param annotationObject object of thr GRanges S4 class
#' @param annotationName (optional) Name of the annotation. Will be used to 
#' create folders. 
#' 
#' @details If no annotationName is provided the name of the file is used to 
#' create folders.
#' 
#' @return An object of the class UserMetadata
#' 
#' @examples {
#' user <- new("UserMetadata")
#' annotation_object <- rtracklayer::import(system.file("extdata", 
#' "annotation.gtf", package = "BgeeCall"))
#' user <- setAnnotationFromObject(user, annotation_object,
#'                                  "annotation_name")
#' }
#' 
#' @export
#' @docType methods
#' @rdname setAnnotationFromObject
#' 
setGeneric("setAnnotationFromObject", 
           function(userObject, annotationObject, 
                    annotationName) {
    standardGeneric("setAnnotationFromObject")
})

#' 
#' @title Set transcriptome_object of one UserMetadata object
#' 
#' @description Method of the class UserMetadata. Set transcriptome_object of 
#' one UserMetadata object 
#' by using one DNAStringSet object as input.
#' 
#' @details Please use a DNAStringSet object as input. This class is defined 
#' in the Biostrings package
#' 
#' @param userObject UserMetadata object
#' @param transcriptomeObject Object of the DNAStringSet S4 class
#' @param transcriptomeName Name of the transcriptome. Will be used to create 
#' transcriptome folders.
#' 
#' @return an object of UserMetadata
#' 
#' @examples {
#' user <- new("UserMetadata")
#' transcriptome_object <- Biostrings::readDNAStringSet(
#'     system.file("extdata", "transcriptome.fa", package = "BgeeCall"))
#' user <- setTranscriptomeFromObject(user, 
#'                  transcriptome_object,
#'                  "transcriptome_name")
#' }
#' 
#' @export
#' @docType methods
#' @rdname setTranscriptomeFromObject
#' 
setGeneric("setTranscriptomeFromObject", function(userObject, 
                                                  transcriptomeObject, 
                                                  transcriptomeName) {
    standardGeneric("setTranscriptomeFromObject")
})

#' 
#' @title Set annotation_object of one UserMetadata object
#' 
#' @description Method of the class UserMetadata. Set annotation_object of 
#' one UserMetadata object  by providing the path to a fasta transcriptome file.
#' 
#' @param userObject The UserMetadata object
#' @param annotationPath Absolute path to the annotation file
#' @param annotationName (optional) Name of the annotation. Will be used to 
#' create folders. 
#' 
#' @details If no annotationName is provided the name of the annotation file 
#' will be used to create folders. 
#' 
#' @return An object of the class UserMetadata
#' 
#' @export
#' @docType methods
#' @rdname setAnnotationFromFile
#' 
#' @examples {
#' # path to gtf annotation file
#' annotation_file <- system.file("extdata", "annotation.gtf", package = "BgeeCall")
#' user <- new("UserMetadata")
#' user <- setAnnotationFromFile(user, annotation_file,
#'                              "annotation_name")
#' }
#'
setGeneric(name="setAnnotationFromFile", 
           def=function(userObject, annotationPath, annotationName) {
               standardGeneric("setAnnotationFromFile")
})

#' 
#' @title Set transcriptome_object of one UserMetadata object
#' 
#' @description Method of the class UserMetadata. Set transcriptome_object of 
#' one UserMetadata object 
#' by providing the path to a fasta transcriptome file.
#' 
#' @param userObject The UserMetadata object
#' @param transcriptomePath Absolute path to the transcriptome file
#' @param transcriptomeName (optional) Name of the trancriptome. Will be used to 
#' create folders. 
#' 
#' @details If no transcriptomeName is provided the name of the transcriptome file 
#' will be used to create folders. 
#' 
#' @return An object of the class UserMetadata
#' 
#' @export
#' @docType methods
#' @rdname setTranscriptomeFromFile
#' 
#' @examples {
#' transcriptome_path <- system.file("extdata", "transcriptome.fa", package = "BgeeCall")
#' user <- new("UserMetadata")
#' user <- setTranscriptomeFromFile(user, transcriptome_path,
#'                                  "transcriptome_name")
#' }
#'
setGeneric(name="setTranscriptomeFromFile", 
           def=function(userObject, transcriptomePath, transcriptomeName) {
    standardGeneric("setTranscriptomeFromFile")
})

#' @title `simple_arborescence` Setter
#' 
#' @description Set value of the `simple_arborescence` slot
#' 
#' @param userObject The UserMetadata object
#' @param simpleArborescence boolean defining if output files will be
#' created a simple arborescence (TRUE) or not (FALSE)
#' 
#' @return An object of the class UserMetadata with new `simple_arborescence`
#'  value
#' 
#' @export
#' @docType methods
#' @rdname setSimpleArborescence
#' 
#' @examples {
#' user <- new("UserMetadata")
#' user <- setSimpleArborescence(user, FALSE)
#' }
#'
setGeneric(name="setSimpleArborescence", 
           def=function(userObject, simpleArborescence) {
               standardGeneric("setSimpleArborescence")
           })

#' @title `simple_arborescence` Getter
#' 
#' @description Get value of the `simple_arborescence` slot
#' 
#' @param userObject The UserMetadata object
#' 
#' @return the value of the `simple_arborescence` slot of the object
#' 
#' @export
#' @docType methods
#' @rdname getSimpleArborescence
#' 
#' @examples {
#' user <- new("UserMetadata")
#' simple_arborescence <- getSimpleArborescence(user)
#' }
#'
setGeneric(name="getSimpleArborescence", def=function(userObject) {
    standardGeneric("getSimpleArborescence")
})

#' @title `run_ids` Setter
#' 
#' @description Method of the class UserMetadata. Set run_ids of 
#' one UserMetadata object by providing the id of all wanted runs
#' 
#' @param userObject The UserMetadata object
#' @param runIds id of all wanted runs
#' 
#' @return An object of the class UserMetadata
#' 
#' @export
#' @docType methods
#' @rdname setRunIds
#' 
#' @examples {
#' user <- new("UserMetadata")
#' user <- setRunIds(user, c("RUN_1", "RUN_2"))
#' }
#'
setGeneric(name="setRunIds", 
           def=function(userObject, runIds) {
               standardGeneric("setRunIds")
           })

#' @title `run_ids` Getter
#' 
#' @description Get value of the `run_ids` slot
#' 
#' @param userObject The UserMetadata object
#' 
#' @return the value of the `run_ids` slot of the object
#' 
#' @export
#' @docType methods
#' @rdname getRunIds
#' 
#' @examples {
#' user <- new("UserMetadata")
#' run_ids <- getRunIds(user)
#' }
#'
setGeneric(name="getRunIds", def=function(userObject) {
    standardGeneric("getRunIds")
})

#' @title `working_path` Setter
#' 
#' @description Set value of the `working_path` slot
#' 
#' @param userObject The UserMetadata object
#' @param workingPath path to the directory wanted as `working_path`
#' 
#' @return An object of the class UserMetadata with new `working_path`
#'  value
#' 
#' @export
#' @docType methods
#' @rdname setWorkingPath
#' 
#' @examples {
#' user <- new("UserMetadata")
#' user <- setWorkingPath(user, getwd())
#' }
#'
setGeneric(name="setWorkingPath", 
           def=function(userObject, workingPath) {
    standardGeneric("setWorkingPath")
})


#' @title `working_path` Getter
#' 
#' @description Get value of the `working_path` slot
#' 
#' @param userObject The UserMetadata object
#' 
#' @return the value of the `working_path` slot of the object
#' 
#' @export
#' @docType methods
#' @rdname getWorkingPath
#' 
#' @examples {
#' user <- new("UserMetadata")
#' working_path <- getWorkingPath(user)
#' }
#'
setGeneric(name="getWorkingPath", def=function(userObject) {
    standardGeneric("getWorkingPath")
})

#' @title `rnaseq_lib_path` Setter
#' 
#' @description Set value of the `rnaseq_lib_path` slot
#' 
#' @param userObject The UserMetadata object
#' @param rnaSeqLibPath path to the directory wanted as `rnaseq_lib_path`
#' 
#' @return An object of the class UserMetadata with new `rnaseq_lib_path`
#'  value
#' 
#' @export
#' @docType methods
#' @rdname setRNASeqLibPath
#' 
#' @examples {
#' user <- new("UserMetadata")
#' user <- setRNASeqLibPath(user, getwd())
#' }
#'
setGeneric(name="setRNASeqLibPath", 
           def=function(userObject, rnaSeqLibPath) {
    standardGeneric("setRNASeqLibPath")
})


#' @rdname setWorkingPath
#' @aliases setWorkingPath,userMetadata,character
setMethod(f="setWorkingPath",
          signature=c(userObject = "UserMetadata", 
                      workingPath = "character"), 
          definition=function(userObject, workingPath) {
              userObject@working_path <- workingPath
              return(userObject)
          })


#' @rdname getWorkingPath
#' @aliases getWorkingPath,userMetadata
setMethod(f="getWorkingPath", 
          signature=c(userObject = "UserMetadata"), 
          definition=function(userObject) {
              return(userObject@working_path)
          })

#' @rdname setRNASeqLibPath
#' @aliases setRNASeqLibPath,userMetadata,character
setMethod(f="setRNASeqLibPath", 
          signature=c(userObject = "UserMetadata", 
                      rnaSeqLibPath = "character"), 
          definition=function(userObject, rnaSeqLibPath) {
              userObject@rnaseq_lib_path <- rnaSeqLibPath
              return(userObject)
          })

#' @rdname setTranscriptomeFromFile
#' @aliases setTranscriptomeFromFile,userMetadata,character,missing
setMethod(f="setTranscriptomeFromFile", 
signature=c(userObject = "UserMetadata", 
            transcriptomePath = "character", 
            transcriptomeName = "missing"), 
definition=function(userObject, transcriptomePath, 
                    transcriptomeName) {
    return(setTranscriptomeFromFile(userObject, transcriptomePath, ""))   
})
    

#' @rdname setTranscriptomeFromFile
#' @aliases setTranscriptomeFromFile,userMetadata,character,character
setMethod(f="setTranscriptomeFromFile", 
          signature=c(userObject = "UserMetadata", 
                      transcriptomePath = "character", 
                      transcriptomeName = "character"), 
          definition=function(userObject, transcriptomePath, 
                              transcriptomeName) {
              if(typeof(transcriptomePath) == "character") {
                  if(file.exists(transcriptomePath)) {
                      userObject@transcriptome_object <- 
                          readDNAStringSet(transcriptomePath)
                  } else {
                      stop(paste0("file ", transcriptomePath, 
                                  " does not exist. Should be the full path to your 
                            transcriptome file"))
                  }
              }
              if (nchar(transcriptomeName) == 0) {
                  userObject@transcriptome_name <- basename(transcriptomePath)
              } else {
                  userObject@transcriptome_name <- transcriptomeName
              }
              return(userObject)
          })

#' @rdname setAnnotationFromFile
#' @aliases setAnnotationFromFile,userMetadata,character,character
setMethod(f="setAnnotationFromFile", 
          signature=c(userObject = "UserMetadata", 
                      annotationPath = "character", annotationName = "missing"),
          definition=function(userObject, annotationPath, annotationName) {
              return(setAnnotationFromFile(userObject, annotationPath, ""))
          })

#' @rdname setAnnotationFromFile
#' @aliases setAnnotationFromFile,userMetadata,character,character
setMethod(f="setAnnotationFromFile", 
          signature=c(userObject = "UserMetadata", 
                      annotationPath = "character", annotationName = "character"),
          definition=function(userObject, annotationPath, annotationName) {
              if(typeof(annotationPath) == "character") {
                  if(file.exists(annotationPath)) {
                      userObject@annotation_object <- rtracklayer::import(annotationPath)
                  } else {
                      stop(paste0("file ", annotationPath, " does not exist. 
                                  Should be the full path to your annotation file"))
                  }
                  }
              if (nchar(annotationName) == 0) {
                  userObject@annotation_name <- basename(annotationPath)
              } else {
                  userObject@annotation_name <- annotationName
              }
              return(userObject)
              })

#' @rdname setTranscriptomeFromObject
#' @aliases setTranscriptomeFromObject,userMetadata,DNAStringSet,character
setMethod(f="setTranscriptomeFromObject", 
          signature=c(userObject = "UserMetadata", 
                      transcriptomeObject = "DNAStringSet", 
                      transcriptomeName = "character"),
          definition=function(userObject, transcriptomeObject, 
                              transcriptomeName) {
              if(isS4(transcriptomeObject)) {
                  userObject@transcriptome_object <- transcriptomeObject
                  userObject@transcriptome_name <- transcriptomeName
              } else {
                  stop("Please provide an object imported using 
                       rtracklayer::import()")
              }
              return(userObject)
              })

#' @rdname setAnnotationFromObject
#' @aliases setAnnotationFromObject,userMetadata,GRanges,character
setMethod(f="setAnnotationFromObject", 
          signature=c(userObject = "UserMetadata", 
                      annotationObject = "GRanges", 
                      annotationName = "character"),
          definition=function(userObject, annotationObject, annotationName = "") {
              if(isS4(annotationObject) && (nchar(annotationName) != 0)) {
                  userObject@annotation_object <- annotationObject
                  userObject@annotation_name <- annotationName
              } else {
                  stop("Please provide an object imported using 
                       rtracklayer::import()")
              }
              return(userObject)
              })

#' @rdname setRunIds
#' @aliases setRunIds,userMetadata,character
setMethod(f="setRunIds",
          signature=c(userObject = "UserMetadata", 
                      runIds = "character"), 
          definition=function(userObject, runIds) {
              userObject@run_ids <- runIds
              return(userObject)
          })

#' @rdname getRunIds
#' @aliases getRunIds,userMetadata
setMethod(f="getRunIds", 
          signature=c(userObject = "UserMetadata"), 
          definition=function(userObject) {
              return(userObject@run_ids)
          })


#' @rdname setSimpleArborescence
#' @aliases setSimpleArborescence,userMetadata,logical
setMethod(f="setSimpleArborescence",
          signature=c(userObject = "UserMetadata", 
                      simpleArborescence = "logical"), 
          definition=function(userObject, simpleArborescence) {
              userObject@simple_arborescence <- simpleArborescence
              return(userObject)
          })

#' @rdname getSimpleArborescence
#' @aliases getSimpleArborescence,userMetadata
setMethod(f="getSimpleArborescence", 
          signature=c(userObject = "UserMetadata"), 
          definition=function(userObject) {
              return(userObject@simple_arborescence)
          })
