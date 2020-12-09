#' @title generate gene expression calls with BgeeCall
#' 
#' @keywords RNA-Seq, abundance quantification, present/absent calls
#'
#' @description BgeeCall allows to generate present/absent gene expression calls without 
#' using an arbitrary cutoff like TPM<1. Calls are generated based on reference intergenic 
#' sequences. These sequences are generated based on expression of all RNA-Seq libraries 
#' of each species integrated in Bgee (https://bgee.org).
#' 
#' @details Thes most important functions are :
#' \itemize{
#'   \item generate_calls_workflow : generate present/absent calls on a computer
#'   \item generate_slurm_indexes : generate kallisto indexes for a list of libraries 
#'   on a cluster with slurm queuing system.
#'   \item generate_slurm_calls : generate present/absent calls for a list of libraries 
#'   on a cluster with slurm queuing system. Indexes have to be generated first with the function
#'   `generate_slurm_indexes`
#'   \item merging_libraries : merge calls from different libraries corresponding to the 
#'   same condition. Extremely useful if different libraries correspond to same condition 
#'   (e.g. same anatomical entity from same species)
#' }
#' For more details please have a look at the vignette with the command \strong{vignette("BgeeCall")}
#' 
#' @seealso https://github.com/BgeeDB/BgeeCall
#' 
#' @author Julien Wollbrett
#'
#' @docType package
#' @name BgeeCall
NULL
