## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval = FALSE------------------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("BgeeCall")

## ---- message = FALSE, warning = FALSE----------------------------------------
library(BgeeCall)

## ---- eval=FALSE--------------------------------------------------------------
#  library("ShortRead")
#  # keep 48.000 reads
#  sampler <- FastqSampler(file.path("absolute_path","/SRX099901/SRR350955.fastq.gz"), 48000)
#  set.seed(1); SRR350955 <- yield(sampler)
#  writeFastq(object = SRR350955,
#             file =file.path( "absolute_path","SRX099901_subset", "SRR350955_subset.fastq.gz"),
#             mode = "w", full = FALSE, compress = TRUE)

## ---- message = FALSE, warning = FALSE----------------------------------------
ah <- AnnotationHub::AnnotationHub()
ah_resources <- AnnotationHub::query(ah, c("Ensembl", "Caenorhabditis elegans", "84"))
annotation_object <- ah_resources[["AH50789"]]
# remove MtDNA not tag as C. elegans genome
annotation_object <- GenomeInfoDb::dropSeqlevels(annotation_object, "MtDNA", "coarse")
transcriptome_object <- rtracklayer::import.2bit(ah_resources[["AH50453"]])

## ---- message = FALSE, warning = FALSE----------------------------------------
# create an object of class UserMetadata and specify the species ID
user_BgeeCall <- new("UserMetadata", species_id = "6239")
# import annotation and transcriptome in the user_BgeeCall object
# it is possible to import them using an S4 object (GRanges, DNAStringSet) or a file (gtf, fasta)
user_BgeeCall <- setAnnotationFromObject(user_BgeeCall, annotation_object, "WBcel235_84")
user_BgeeCall <- setTranscriptomeFromObject(user_BgeeCall, transcriptome_object, "WBcel235")
# provide path to the directory of your RNA-Seq library
user_BgeeCall <- setRNASeqLibPath(user_BgeeCall, 
                                  system.file("extdata", 
                                              "SRX099901_subset", 
                                              package = "BgeeCall"))

## ---- eval = FALSE------------------------------------------------------------
#  calls_output <- generate_calls_workflow(userMetadata = user_BgeeCall)

## ---- echo=FALSE--------------------------------------------------------------
user_BgeeCall <- setWorkingPath(user_BgeeCall, system.file("extdata", package = "BgeeCall"))
user_BgeeCall<- setSimpleArborescence(user_BgeeCall, TRUE)
calls_output <- generate_presence_absence(myUserMetadata = user_BgeeCall)

## ---- message = FALSE, warning = FALSE----------------------------------------
head(read.table(calls_output$calls_tsv_path, header = TRUE), n = 5)

## ---- message = FALSE, warning = FALSE----------------------------------------
read.table(calls_output$cutoff_info_file_path)

## ---- message = FALSE, warning = FALSE----------------------------------------
head(read.table(calls_output$abundance_tsv, header = TRUE), n = 5)

## ---- eval = FALSE------------------------------------------------------------
#  openPDF(calls_output$TPM_distribution_path)

## ---- message = FALSE, warning = FALSE----------------------------------------
read.table(calls_output$S4_slots_summary, header = TRUE, sep = "\t")

## ---- eval=FALSE--------------------------------------------------------------
#  calls_output <- generate_calls_workflow(userFile = "path_to_your_file.tsv")

## ---- eval = FALSE------------------------------------------------------------
#  # generate kallisto indexes
#  generate_slurm_indexes(userFile = "path_to_your_file.tsv")
#  # generate expression calls
#  generate_slurm_calls(userFile = "path_to_your_file.tsv")

## -----------------------------------------------------------------------------
list_intergenic_release()

## -----------------------------------------------------------------------------
# create BgeeMetadata object and define one reference intergenic release
bgee <- new("BgeeMetadata", intergenic_release = "0.1")
# change the reference intergenic release of your BgeeMetadata object
bgee <- setIntergenicRelease(bgee, "0.2")

## -----------------------------------------------------------------------------
list_bgee_ref_intergenic_species(myBgeeMetadata = bgee)

## -----------------------------------------------------------------------------
list_community_ref_intergenic_species()

## ---- eval=FALSE--------------------------------------------------------------
#  # create a BgeeMetadata object using the community release
#  bgee <- new("BgeeMetadata", intergenic_release = "community")
#  calls_output <- generate_calls_workflow(bgeeMetadata = bgee, userMetadata = user_BgeeCall)

## ---- eval=FALSE--------------------------------------------------------------
#  bgee <- new("BgeeMetadata", intergenic_release = "custom")
#  user_BgeeCall@custom_intergenic_path = "path/to/custom/ref_intergenic.fa.gz"
#  calls_output <- generate_calls_workflow(bgeeMetadata = bgee, userMetadata = user_BgeeCall)

## ---- eval=FALSE--------------------------------------------------------------
#  kallisto <- new("KallistoMetadata", txOut = TRUE)
#  calls_output <- generate_calls_workflow(myAbundanceMetadata = kallisto, userMetadata = user_BgeeCall)

## ---- eval=FALSE--------------------------------------------------------------
#  kallisto <- new("KallistoMetadata", download_kallisto = TRUE)
#  calls_output <- generate_calls_workflow(myAbundanceMetadata = kallisto, userMetadata = user_BgeeCall)

## ---- eval=FALSE--------------------------------------------------------------
#  kallisto <- new("KallistoMetadata", single_end_parameters = "-t 3 --single -l 150 -s 30", pair_end_parameters = "-t 2 -b --seed 36")
#  calls_output <- generate_calls_workflow(myAbundanceMetadata = kallisto, userMetadata = user_BgeeCall)

## ---- eval=FALSE--------------------------------------------------------------
#  # libraries with reads smaller than 70bp will use the index with kmer size = 15
#  kallisto <- new("KallistoMetadata", read_size_kmer_threshold = 70)
#  calls_output <- generate_calls_workflow(myAbundanceMetadata = kallisto, userMetadata = user_BgeeCall)

## ---- eval=FALSE--------------------------------------------------------------
#  # RNA-Seq run SRR350955_subsetof from the RNA-Seq library will be used to generate the calls
#  user_BgeeCall <- setRunIds(user_BgeeCall, c("SRR350955_subset"))
#  calls_output <- run_from_object(myUserMetadata = user_BgeeCall)

## -----------------------------------------------------------------------------
kallisto <- new("KallistoMetadata", cutoff = 0.1)

## -----------------------------------------------------------------------------
# use intergenic approach with default cutoff ratio 0.05
kallisto <- new("KallistoMetadata", cutoff_type = "intergenic")
# use intergenic approach with cutoff ratio 0.01
kallisto <- new("KallistoMetadata", cutoff_type = "intergenic", cutoff = 0.01)

## -----------------------------------------------------------------------------
user_BgeeCall@verbose <- FALSE

## ---- eval=FALSE--------------------------------------------------------------
#  user_BgeeCall <- setSimpleArborescence(userObject = user_BgeeCall, simpleArborescence = FALSE)
#  calls_output <- run_from_object(myUserMetadata = user_BgeeCall)

## -----------------------------------------------------------------------------
user_BgeeCall <- setOutputDir(user_BgeeCall, "path/to/calls/for/this/library/")

## ---- eval=FALSE--------------------------------------------------------------
#  # run 50 jobs in parallel
#  generate_slurm_indexes(userFile = "path/to/file.tsv", nodes = 50)

## ---- eval=FALSE--------------------------------------------------------------
#  # create temporary files but do not submit the jobs
#  generate_slurm_indexes(userFile = "path/to/file.tsv", submit = FALSE)

## ---- eval=FALSE--------------------------------------------------------------
#  # add slurm options to the sbatch script
#  slurm_options_index <- list(account = "account", time = "2:00:00", partition = "partition", mem = "30G")
#  generate_slurm_indexes(userFile = "path/to/file.tsv", slurm_options = slurm_options_index)

## ---- eval=FALSE--------------------------------------------------------------
#  # load R 3.6.1 and kallisto in a cluster environment where software has to loaded manually
#  modules <- c("module add R/3.6.1;", "module add kallisto;")
#  generate_slurm_indexes(userFile = "path/to/file.tsv", modules = modules)

## ---- eval=FALSE--------------------------------------------------------------
#  # create BgeeCall objects and use them to generate indexes
#  kallistoMetadata <- new("KallistoMetadata", download_kallisto=TRUE)
#  userMetadata <- new("UserMetadata", working_path = "/path/to/working/dir")
#  bgeeMetadata <- new("BgeeMetadata", intergenic_release = "0.1")
#  generate_slurm_indexes(userFile = "path/to/file.tsv", kallistoMetadata = kallistoMetadata, bgeeMetadata = bgeeMetadata, userMetadata = userMetadata)

## ---- eval=FALSE--------------------------------------------------------------
#  mergingLibraries <- merging_libraries(userFile = "path/to/userFile.tsv", approach = "BH", condition = "species_id", cutoff = 0.05, outDir = "path/to/output_directory/")

## ---- eval=FALSE--------------------------------------------------------------
#  # merging libraries
#  mergingLibraries <- merging_libraries(userFile = "path/to/userFile.tsv", approach = "fdr_inverse", condition = "species_id", cutoff = 0.05, outDir = "path/to/output_directory/")

## ---- eval=FALSE--------------------------------------------------------------
#  
#  # merging libraries where p-values were calculated
#  mergingLibraries_pValue <- merging_libraries(userFile = "path/to/userFile.tsv", approach = "BH", condition = c("species_id", "anatEntity", "strain"), cutoff = 0.005, outDir = "path/to/output_directory/")
#  
#  # merging libraries where q-values were calculated
#  mergingLibraries_qValue <- merging_libraries(userFile = "path/to/userFile.tsv", approach = "fdr_inverse", condition = c("species_id", "devStage" ,"anatEntity", "strain"), cutoff = 0.001, outDir = "path/to/output_directory/")
#  

## ----sessioninfo--------------------------------------------------------------
sessionInfo()

## ----cleanup_after, echo=FALSE, message=FALSE, warning=FALSE------------------
unlink(BgeeCall:::get_kallisto_dir_path(kallisto, user_BgeeCall), recursive = TRUE)
unlink(file.path(getWorkingPath(user_BgeeCall), paste0(getIntergenicPrefix(bgee), "*")), recursive = TRUE)

