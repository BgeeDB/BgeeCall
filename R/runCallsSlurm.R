# regroup all function used to easily run the
# RNA-Seq calls presence/absence pipeline in slurm cluster
#' 
#' @title Generate all indexes for the abundance quantification step 
#' 
#' @description Check all unique lines of the input file to check which
#' indexes have to be generated beore running all abundance quantification.
#' This function is meant to be used with a cluster where the Slurm queuing 
#' system is installed. This step has to be run before the quantification
#' otherwise indexes will be created for each abundance quantification. This
#' will slow down the abundance quantification and can generate errors when 
#' writting the same file at the same time from different nodes.
#' 
#' @param abundanceMetadata A Reference Class BgeeMetadata object (optional)
#' allowing to tune your gene quantification abundance analyze. If no object is
#' provided a new one will be created with default values.
#' @param bgeeMetadata A Reference Class BgeeMetadata object (optional)
#' allowing to choose the version of reference intergenic sequences. If no object
#' is provided a new one will be created with default values.
#' @param userFile Path to the file where each line corresponds to one abundance
#' quantification to be run. The structure of the file is the same than the 
#' `userFile` used as input of  the `generate_calls_workflow` function.
#' A template of this file can be loaded with the command :
#' ```inputFile <- read.table(system.file("userMetadataTemplate.tsv", package = "BgeeCall"), 
#'                 header = TRUE)```
#' It is important to keep the same column names.
#' @param submit_sh_template A template of the bash script used to submit the jobs.
#' By default the submition script provided by rslurm is used. Modify only if module 
#' dependancies have to be added  (like kallisto or R)
#' @param modules A list of modules you want to load in the invironment. Should stay
#' empty except if you need to load R and/or kallisto (e.g module add R)
#' @param slurm_options A named list of options recognized by sbatch. More details in
#' the documentation of the rslurm::slurm_apply function
#' @param rscript_path The location of the Rscript command. If not specified, defaults 
#' to the location of Rscript within the R installation being run.
#' 
#' @return generate index files
#' 
#' @export
#' @import rslurm
#' 
#' @examples 
#' \dontrun{
#' # use function with all default values
#' userFile <- "/path/to/userList.tsv"
#' sjobs <- generate_slurm_indexes(userFile = userFile)
#' }

generate_slurm_indexes <- function(kallistoMetadata = new("KallistoMetadata"), 
                                 bgeeMetadata = new("BgeeMetadata"), 
                                 userMetadata = new("UserMetadata"), userFile, 
                                 submit_sh_template = NULL, slurm_options = NULL,
                                 rscript_path = NULL, modules = NULL) {
  user_df <- read.table(file = userFile, header = TRUE, sep = "\t")
  read_length_threshold <- kallistoMetadata@read_size_kmer_threshold
  user_df[user_df$reads_size < read_length_threshold, "index_type"] <- "short"
  user_df[user_df$reads_size >= read_length_threshold, "index_type"] <- "normal"
  # detect unique conbination of index to generate (transcriptome, index type 
  # and species)
  unique_df <- user_df[!duplicated(user_df[c("species_id", "index_type", 
                                             "transcriptome_path")]),]
  unique_df$index_type <- NULL
  
  #use both rscript_path and modules variables to tune the bash template script
  rscript_path <- generateRScriptPath(rscript_path, modules)
  
  #define function used to generate kallisto index. This function has to be defined inside of the function
  # calling slurm_apply (enclosing environment (http://adv-r.had.co.nz/Environments.html#function-envs))
  index_wrapper <- function(species_id, run_ids, reads_size, rnaseq_lib_path, transcriptome_path, 
                            annotation_path, working_path, output_directory, simple_arborescence) {
    userMetadata <- new("UserMetadata", speciesId = species_id, run_ids = run_ids, 
                        reads_size = reads_size, rnaseq_lib_path = rnaseq_lib_path, 
                        output_directory = output_directory)
    if(!is.null(working_path)) {
      userMetadata@working_path <- working_path
    }
    if(!is.null(simple_arborescence)) {
      userMetadata@simple_arborescence <- simple_arborescence
    }
    userMetadata <- setTranscriptomeFromFile(userMetadata, transcriptomePath = transcriptome_path)
    merge_transcriptome_and_intergenic(myKallistoMetadata = kallistoMetadata, myBgeeMetadata = bgeeMetadata, 
                                       myUserMetadata = userMetadata)
    create_kallisto_index(myKallistoMetadata = kallistoMetadata, myBgeeMetadata = bgeeMetadata, 
                          myUserMetadata = userMetadata)
    message("run index generation for species ",userMetadata@species_id)
  }
  
  sjobs <- rslurm::slurm_apply(f = index_wrapper, params = unique_df, jobname = "generate_index", 
                              nodes = 1, cpus_per_node = 1, submit = TRUE, 
                              add_objects = c("kallistoMetadata", "bgeeMetadata", "userMetadata"), 
                              sh_template = submit_sh_template, rscript_path = rscript_path, slurm_options = slurm_options)
  return(sjobs)
}




abundance_wrapper <- function(species_id, run_ids, reads_size, rnaseq_lib_path, transcriptome_path, 
                          annotation_path, working_path, output_directory, simple_arborescence) {
  userMetadata <- new("UserMetadata", speciesId = species_id, run_ids = run_ids, 
                      reads_size = reads_size, rnaseq_lib_path = rnaseq_lib_path, 
                      working_path = working_path, output_directory = output_directory,
                      simple_arborescence = simple_arborescence)
  userMetadata <- setTranscriptomeFromFile(userMetadata, transcriptomePath = transcriptome_path)
  run_from_object(myKallistoMetadata = kallistoMetadata, myBgeeMetadata = bgeeMetadata, 
                        myUserMetadata = userMetadata)
  message("generate present/absent calls")
}

# internal function hacking the rscript_path variable of rslurm to automatically add
# module dependancies in the bash script generated by rslurm
generateRScriptPath <- function(rscript_path, modules) {
  if(is.null(rscript_path)) {
    rscript_path <- file.path(R.home("bin"), "Rscript")
  }
  if (!is.null(modules)) {
    for (i in seq(modules)) {
      rscript_path <- paste0(modules[i],"\n",rscript_path)
    }
  }
  return(rscript_path)
}

