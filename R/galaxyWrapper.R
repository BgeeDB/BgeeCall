
#' @title generate calls for galaxy
#' 
#' @description Function generating calls for galaxy. It is a simple adaptor of the
#' function ```generate_calls_workflow``` adapted to galaxy requirements.
#' 
#' @param species_id The NCBI taxon ID of the species
#' @param run_ids The subset of runs from the library used to generate the calls.
#' By default all runs are used.
#' @param rnaseq_lib_path Path to the directory containing the fastq files of the library
#' @param reads_size mean size of the reads. (Default = 51)
#' @param transcriptome_path Path to the transcriptome fastq file
#' @param annotation_path path to the gtf file containing gene annotation
#' @param output_file file containing the path to all result files generated
#' @param custom_intergenic_path (optional) To use if the "custom" reference intergenic
#'  release has been selected. Provides the path to the reference intergenic file
#' 
#' 
#' 

galaxy_calls <- function(species_id = NULL, run_ids = NULL, rnaseq_lib_path = NULL,
                         reads_size = 51, transcriptome_path = NULL, annotation_path = NULL,
                         output_dir = NULL, ignore_tx_version = FALSE,
                         custom_intergenic_path = NULL, userFile = NULL) {
  if (is.null(species_id) || is.null(rnaseq_lib_path) || is.null(transcriptome_path)
      || is.null(annotation_path)) {
    stop("Did not provide the proper input to generate the call. You should provide",
         " either:\n\t - a species ID, the path to your rnaseq library, the path to",
         " the transcriptome, the path to the gene annotation and an output directory",
         " or \n\t - the path to a tsv file containing at least 4 columns as described",
         " in the documentation.")
  }
  myAbundanceMetadata <- new("KallistoMetadata")
  myBgeeMetadata <- new("BgeeMetadata")
  myUserMetadata <- new("UserMetadata", species_id = species_id, run_ids = run_ids,
                       reads_size = reads_size, rnaseq_lib_path = rnaseq_lib_path)
  setAnnotationFromFile(myUserMetadata, annotation_path)
  setTranscriptomeFromFile(myUserMetadata, transcriptome_path)
  check_tx_version <- should_ignore_tx_version(myUserMetadata)
  results <- run_from_object(myAbundanceMetadata = myAbundanceMetadata, myUserMetadata = myUserMetadata,
                  myBgeeMetadata = myBgeeMetadata, checkTxVersion = check_tx_version)
  write.table(results, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
}