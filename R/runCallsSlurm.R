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
#' This function also generate tx2gene and gene2biotype mapping files.
#' 
#' @param kallistoMetadata A Reference Class KallistoMetadata object (optional)
#' allowing to tune your gene quantification abundance analyze. If no object is
#' provided a new one will be created with default values.
#' @param bgeeMetadata A Reference Class BgeeMetadata object (optional)
#' allowing to choose the version of reference intergenic sequences. If no object
#' is provided a new one will be created with default values.
#' @param userMetadata A Class UserMetadata object (optional).
#' If no object is provided a new one will be created with default values.
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
#' @param submit Whether or not to submit the job to the cluster with sbatch. Default 
#' value is TRUE
#' @param nodes The (maximum) number of cluster nodes to spread the calculation over. 
#' slurm_apply automatically divides params in chunks of approximately equal size to 
#' send to each node. Less nodes are allocated if the parameter set is too small to 
#' use all CPUs on the requested nodes. By default this number is 10.


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
                                 rscript_path = NULL, modules = NULL, submit = TRUE, 
                                 nodes = 10) {
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
                            annotation_path, output_directory=NULL, custom_intergenic_path=NULL) {
    userMetadata@species_id <- as.character(species_id) 
    userMetadata@run_ids <- as.character(run_ids)
    userMetadata@reads_size <- as.integer(reads_size)
    userMetadata@rnaseq_lib_path <- as.character(rnaseq_lib_path)
    outputDir <- as.character(output_directory)
    if (length(outputDir) != 0) {
      userMetadata <- setOutputDir(userMetadata, outputDir)
    }
    customIntergenicPath <- as.character(custom_intergenic_path)
    if (length(customIntergenicPath) != 0) {
      userMetadata@custom_intergenic_path <- customIntergenicPath
    }
    userMetadata <- setTranscriptomeFromFile(userMetadata, transcriptomePath = as.character(transcriptome_path))
    userMetadata <- setAnnotationFromFile(userMetadata, annotationPath =  as.character(annotation_path))
    merge_transcriptome_and_intergenic(myKallistoMetadata = kallistoMetadata, myBgeeMetadata = bgeeMetadata, 
                                       myUserMetadata = userMetadata)
    create_kallisto_index(myKallistoMetadata = kallistoMetadata, myBgeeMetadata = bgeeMetadata, 
                          myUserMetadata = userMetadata)
    #create tx2gene file
    tx2gene <- create_tx2gene(myAbundanceMetadata = kallistoMetadata, myBgeeMetadata = bgeeMetadata, 
                              myUserMetadata = userMetadata)
    #create gene2biotype file
    create_gene_to_biotype(myAbundanceMetadata = kallistoMetadata, myBgeeMetadata = bgeeMetadata, 
                           myUserMetadata = userMetadata)
  }
  
  sjobs <- rslurm::slurm_apply(f = index_wrapper, params = unique_df, jobname = "generate_index", 
                              nodes = nodes, cpus_per_node = 1, submit = submit, 
                              add_objects = c("kallistoMetadata", "bgeeMetadata", "userMetadata"), 
                              sh_template = submit_sh_template, rscript_path = rscript_path, slurm_options = slurm_options)
  return(sjobs)
}

#' 
#' @title Generate present/absent calls on slurm queuing system 
#' 
#' @description This function is meant to be used with a cluster where the Slurm queuing 
#' system is installed. It processes all steps to generate present/absent calls at RNA-Seq 
#' library level. This function does not generate the kallisto indexes. If they are not 
#' already generated please run function ```generate_slurm_indexes``` first.
#' Steps of present/absent gene expression calls generation are :
#' \itemize{
#'   \item Quantifying abundances of transcripts from RNA-Seq libraries
#'   \item Summarizing abundance at gene level 
#'   \item generate present/absent expression calls
#' } 
#' 
#' 
#' @param kallistoMetadata A Reference Class KallistoMetadata object (optional)
#' allowing to tune your gene quantification abundance analyze. If no object is
#' provided a new one will be created with default values.
#' @param bgeeMetadata A Reference Class BgeeMetadata object (optional)
#' allowing to choose the version of reference intergenic sequences. If no object
#' is provided a new one will be created with default values.
#' @param userMetadata A Class UserMetadata object (optional).
#' If no object is provided a new one will be created with default values.
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
#' @param submit Whether or not to submit the job to the cluster with sbatch. Default 
#' value is TRUE
#' @param nodes The (maximum) number of cluster nodes to spread the calculation over. 
#' slurm_apply automatically divides params in chunks of approximately equal size to 
#' send to each node. Less nodes are allocated if the parameter set is too small to 
#' use all CPUs on the requested nodes. By default this number is 10.


#' 
#' @return generate calls
#' 
#' @export
#' @import rslurm
#' 
#' @examples 
#' \dontrun{
#' # use function with all default values
#' userFile <- "/path/to/userList.tsv"
#' sjobs <- generate_slurm_calls(userFile = userFile)
#' }

generate_slurm_calls <- function(kallistoMetadata = new("KallistoMetadata"), 
                                   bgeeMetadata = new("BgeeMetadata"), 
                                   userMetadata = new("UserMetadata"), userFile, 
                                   submit_sh_template = NULL, slurm_options = NULL,
                                   rscript_path = NULL, modules = NULL, submit = TRUE, 
                                   nodes = 10, checkTxVersion = FALSE) {
  user_df <- read.table(file = userFile, header = TRUE, sep = "\t")

  #use both rscript_path and modules variables to tune the bash template script
  rscript_path <- generateRScriptPath(rscript_path, modules)
  
  #define function used to generate present/absent calls. This function has to be defined inside of the function
  # calling slurm_apply (enclosing environment (http://adv-r.had.co.nz/Environments.html#function-envs))
  calls_wrapper <- function(species_id, run_ids, reads_size, rnaseq_lib_path, transcriptome_path, 
                            annotation_path, output_directory=NULL, custom_intergenic_path=NULL) {
    userMetadata@species_id <- as.character(species_id) 
    userMetadata@run_ids <- check_run_ids(as.character(run_ids))
    userMetadata@reads_size <- as.integer(reads_size)
    userMetadata@rnaseq_lib_path <- as.character(rnaseq_lib_path)
    outputDir <- as.character(output_directory)
    if (length(outputDir) != 0) {
      userMetadata <- setOutputDir(userMetadata, outputDir)
    }
    customIntergenicPath <- as.character(custom_intergenic_path)
    if (length(customIntergenicPath) != 0) {
      userMetadata@custom_intergenic_path <- customIntergenicPath
    }
    userMetadata <- setTranscriptomeFromFile(userMetadata, transcriptomePath = as.character(transcriptome_path))
    userMetadata <- setAnnotationFromFile(userMetadata, annotationPath =  as.character(annotation_path))
    if(checkTxVersion) {
      kallistoMetadata@ignoreTxVersion <- should_ignore_tx_version(userMetadata)
    }
    run_kallisto(myKallistoMetadata = kallistoMetadata, myBgeeMetadata = bgeeMetadata, 
                 myUserMetadata = userMetadata)
    generate_presence_absence(myAbundanceMetadata = kallistoMetadata, myBgeeMetadata = bgeeMetadata, 
                          myUserMetadata = userMetadata)
  }
  
  sjobs <- rslurm::slurm_apply(f = calls_wrapper, params = user_df, jobname = "generate_calls", 
                               nodes = nodes, cpus_per_node = 1, submit = submit, 
                               add_objects = c("kallistoMetadata", "bgeeMetadata", "userMetadata", 
                               "checkTxVersion"), sh_template = submit_sh_template, 
                               rscript_path = rscript_path, slurm_options = slurm_options)
  return(sjobs)
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

