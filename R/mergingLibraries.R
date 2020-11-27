#' @title Checking user file 
#'
#' @description Check if the user provide a file with at least one
#' column (species_id) specified/fill-up to perform the merging of libraries
#'
#' @param userFile A file provided by the user
#'
#' @author Sara Fonseca Costa
#'
#' @return User data frame
#' @return Vector name of columns detected
#' 
#' @import data.table
#' @import sjmisc
#' 
#' @noMd
#' @noRd
#' 
checkUserFile <- function(userFile){
  
  ## read and check columns of input file
  userFile <- fread(userFile)
  idContition <- c("species_id", "anatEntity", "devStage", "sex","strain")
  userFileCheck <- idContition %in% names(userFile)
  
  ## no column is detected in the input file
  if (any(userFileCheck) == FALSE){
    stop("At least one condition (species_id) should exist and be fill-up in order to do the merging.", "\n",
         "Please verify your input file.", "\n")
  } 
  ## exist columns specified but missing species_id column
  else if (any(userFileCheck) == TRUE & userFileCheck[1] == FALSE){
    stop("The condition/s where detected to perform the merging, but the species_id is missing", "\n",
         "Please verify your input file.", "\n")
  }
  ## exist column species_id but missing information (is NA)
  else if (userFileCheck[1] == TRUE & is_empty(userFile$species_id) == TRUE) {
    stop("The row species_id is specified but empty", "\n")
    
  } else {
    conditionDataTable <- idContition[which(userFileCheck=="TRUE")]
  }
  return(list(userFile, conditionDataTable))
}


#' @title Approaches used to merge/combine libraries
#'
#' @description For a set of libraries calculate the adjusted p-values using BH method
#' for each gene id OR calculate the inverse fdr correction to be applied to each gene ID 
#' in a set of libraries that use q-values.
#' The call of genes as present is done based on the assumption that at least one library 
#' have an adjusted p-value < cutoff or have a q-value < fdr_inverse.
#'
#' @param AllFiles List with all libraries to be treated
#' @param approach Approach selected to merge libraries, BH or fdr_inverse
#' @param cutoff Threshold value to be applied to call expressed genes
#'
#' @author Sara Fonseca Costa
#'
#' @return data frame with minimum pValue or qValue detected (dependent of the approach selected) and the calls
#' 
#' @import data.table
#' @import dplyr
#' 
#' @noMd
#' @noRd
#' 
approachesMerging <- function(AllFiles, approach, cutoff){
  
  librariesData <- try(do.call("cbind", AllFiles), silent = TRUE) 
  ## select gene_id
  select_info <- librariesData[,1]
  
  if (approach == "BH"){
    
    ## select all pValues column from all libraries
    select_pValue = librariesData[, grepl("^pValue", names(librariesData))]
    if(is_empty(select_pValue) == TRUE){
      stop("Select the appropriated method for your quantitative metric: BH for p-values OR fdr_inverse for q-values", "\n")
    }
    ## calculate the p.adjusted values using a vector of original pValues for each gene_id
    collect_padjValues <- apply(select_pValue, 1, function (x) p.adjust(x[1:length(select_pValue)], method = "BH"))
    collect_padjValues <- as.data.frame(t(collect_padjValues))
    collect_padjValues$minimum_pValue <- do.call(pmin, c(collect_padjValues, list(na.rm = TRUE))) 
    ## data frame with all information (gene_id + all adjusted pvalues of all libraries)
    allInfo <- data.frame(select_info, collect_padjValues)
    ## provide info just about id and minimum_pValue detected
    allInfo <- allInfo %>% select("select_info", "minimum_pValue")
    colnames(allInfo)[1] <- c("id")
    
    ## perform the calls
    allInfo$call <- ifelse(allInfo$minimum_pValue <= as.numeric(cutoff), "present", "absent")
    ## replace NA (that are values = 0 in the abundance) per absent calls
    allInfo$call[is.na(allInfo$call)] <- "absent"
    
    return(allInfo)
    
  } else if (approach == "fdr_inverse"){
    
    ## select all qValues column from all libraries
    select_qValue = librariesData[, grepl("^qValue", names(librariesData))]
    if(is_empty(select_qValue) == TRUE){
      stop("Select the appropriated method for your quantitative metric: BH for p-values OR fdr_inverse for q-values", "\n")
    }
    select_qValue$minimum_qValue <- do.call(pmin, c(select_qValue, list(na.rm = TRUE))) 
    ## data frame with all information (gene_id + minimum q-Value detected of all libraries)
    allInfo <- data.frame(select_info, select_qValue)
    ## provide info just about id and minimum_qValue detected
    allInfo <- allInfo %>% select("select_info", "minimum_qValue")
    colnames(allInfo)[1] <- c("id")
    
    ## Calculate FDR inverse
    fdrInverse <- (1-((1-as.numeric(cutoff))^(1/(length(select_qValue)-1))))
   
    ## calls using the fdr_inverse
    allInfo$call <- ifelse(allInfo$minimum_qValue <= fdrInverse, "present", "absent")
    ## replace NA (that are values = 0 in the abundance) per absent calls
    allInfo$call[is.na(allInfo$call)] <- "absent"
   
    return(allInfo)
   
  } else {
    stop("The approach provided is not recognized.", "\n", "Please choose: BH (for p-values) or fdr_inverse (for q-values).", "\n")
  }
}


#' @title Calls of expression in combined libraries
#'
#' @description Merging/combine libraries based in a condition specified by the user.
#' The merging can be done using the p-values of the libraries, by applying the BH method,
#' or using the q-values of the libraries using the fdr_inverse method.
#'
#' @param userFile A file provided by the user with correspondent conditions
#' @param approach Approach used to do the merging of libraries
#' @param condition Condition/s where the merging should be done
#' @param cutoff Cutoff that should be applied to call P/A genes
#' @param outDir Directory where the output files should be saved
#'
#' @author Sara Fonseca Costa
#' 
#' @export
#' 
#' @examples 
#' \dontrun{
#' callsMerging_species <- merging_libraries(userFile = 'PATH_USER_FILE', approach = 'BH', condition = 'species_id', cutoff = 0.05, outDir = 'PATH_OUTPUT')
#' callsMerging_species_sex <- merging_libraries(userFile = 'PATH_USER_FILE', approach = 'fdr_inverse', condition = c(species_id, sex), cutoff = 0.01, outDir = 'PATH_OUTPUT')
#' callsMerging_all <- merging_libraries(userFile = 'PATH_USER_FILE', approach = 'fdr_inverse', condition = c(species_id, anatEntity, devStage, sex, strain), cutoff = 0.05, outDir = 'PATH_OUTPUT')
#' }
#' 
#' @return A dataframe containing the minimum quantitative value (p-value or q-value) and 
#' the calls to each gene id for the referent condition.
#' 
#' 
merging_libraries <- function(userFile = NULL, approach = "BH", condition = "species_id", cutoff = 0.05, outDir = NULL) {

  ## check user input
  fileUser <- checkUserFile(userFile = userFile)
  
  ## check if condition/s provided in the user input file match the conditions specified by the user in the condition argument
  userFileconditionDetected <- fileUser[[2]]
  userConditionsToMerge <- condition %in% userFileconditionDetected

  if ('FALSE' %in% userConditionsToMerge){
    stop("One of the parameters provided in the condition argument is not detected in the user input file.", "\n",
         "Check the input file OR the parameters name provide in the condition that should be: ", "\n", "species_id","\n", "anatEntity", "\n", "devStage", "\n", "sex", "\n", "strain", "\n")
  } else {
    
    if ((length(condition)==1) == TRUE & condition[1] == "species_id"){
      
      message("Merging will be done for condition c(species)", "\n")
      
      for (species in unique(fileUser[[1]]$species_id)) {
      ## collect libraries from each species
      librariesSpecies <- fileUser[[1]]$output_directory[fileUser[[1]]$species_id == species]
      ## collect files of individual call
      AllFiles <- list.files(librariesSpecies, pattern="gene_level_abundance\\+calls.tsv", full.names=T, recursive = TRUE)
      message("Using ", length(AllFiles), " libraries for the condition: species_id = ", species ,"\n")
      AllFiles <- lapply(AllFiles, read.delim)
      ## perform the calls and export information of the merging
      callsFile <- approachesMerging(AllFiles = AllFiles, approach = approach, cutoff = cutoff)
      write.table(callsFile, file = paste0(outDir, "/Calls_merging_",approach, "_", "cutoff=", cutoff,"_species=", species, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
      }
      } else if ((length(condition)==2) == TRUE & "species_id" %in% condition[1:length(condition)] & "anatEntity" %in% condition[1:length(condition)]) {
        
        message("Merging will be done for condition c(species, anatEntity)", "\n")
        if(is_empty(fileUser[[1]]$anatEntity) == TRUE){
          stop("The column anatEntity is empty")
        }
        
        for (species in unique(fileUser[[1]]$species_id)) {
          for (anatEntity in unique(fileUser[[1]]$anatEntity)) {
            ## collect libraries from each species and from each anatomical entity
            librariesSpecies <- fileUser[[1]]$output_directory[fileUser[[1]]$species_id == species & fileUser[[1]]$anatEntity == anatEntity]
            if (is_empty(librariesSpecies) == TRUE){
            } else {
            ## collect files of individual call
            AllFiles <- list.files(librariesSpecies, pattern="gene_level_abundance\\+calls.tsv", full.names=T, recursive = TRUE)
            message("Using ", length(AllFiles), " libraries for the condition: c(", species ,",", anatEntity, ")" ,"\n")
            AllFiles <- lapply(AllFiles, read.delim)
            ## perform the calls and export information of the merging
            callsFile <- approachesMerging(AllFiles = AllFiles, approach = approach, cutoff = cutoff)
            write.table(callsFile, file = paste0(outDir, "/Calls_merging_",approach, "_", "cutoff=", cutoff,"_species=", species, "_anatEntity=",anatEntity,".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
            }
          }
        }
        } else if ((length(condition)==2) == TRUE & "species_id" %in% condition[1:length(condition)] & "devStage" %in% condition[1:length(condition)]) {
          message("Merging will be done for condition c(species, devStage)", "\n")
          if(is_empty(fileUser[[1]]$devStage) == TRUE){
            stop("The column devStage is empty")
          }
          
          for (species in unique(fileUser[[1]]$species_id)) {
            for (devStage in unique(fileUser[[1]]$devStage)) {
              ## collect libraries from each species and from each developmental stage
              librariesSpecies <- fileUser[[1]]$output_directory[fileUser[[1]]$species_id == species & fileUser[[1]]$devStage == devStage]
              if (is_empty(librariesSpecies) == TRUE){
              } else {
              ## collect files of individual call
              AllFiles <- list.files(librariesSpecies, pattern="gene_level_abundance\\+calls.tsv", full.names=T, recursive = TRUE)
              message("Using ", length(AllFiles), " libraries for the condition: c(", species ,",", devStage, ")" ,"\n")
              AllFiles <- lapply(AllFiles, read.delim)
              ## perform the calls and export information of the merging
              callsFile <- approachesMerging(AllFiles = AllFiles, approach = approach, cutoff = cutoff)
              write.table(callsFile, file = paste0(outDir, "/Calls_merging_",approach, "_", "cutoff=", cutoff,"_species=", species, "_devStage=", devStage, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
              }
            }
          }
          } else if ((length(condition)==2) == TRUE & "species_id" %in% condition[1:length(condition)] & "sex" %in% condition[1:length(condition)]){
            message("Merging will be done for condition c(species, sex)", "\n")
            if(is_empty(fileUser[[1]]$sex) == TRUE){
              stop("The column sex is empty")
            }
              
              for (species in unique(fileUser[[1]]$species_id)) {
                for (sex in unique(fileUser[[1]]$sex)) {
                  ## collect libraries from each species and from each sex
                  librariesSpecies <- fileUser[[1]]$output_directory[fileUser[[1]]$species_id == species & fileUser[[1]]$sex == sex]
                  if (is_empty(librariesSpecies) == TRUE){
                  } else {
                  ## collect files of individual call
                  AllFiles <- list.files(librariesSpecies, pattern="gene_level_abundance\\+calls.tsv", full.names=T, recursive = TRUE)
                  message("Using ", length(AllFiles), " libraries for the condition: c(", species ,",", sex, ")" ,"\n")
                  AllFiles <- lapply(AllFiles, read.delim)
                  ## perform the calls and export information of the merging
                  callsFile <- approachesMerging(AllFiles = AllFiles, approach = approach, cutoff = cutoff)
                  write.table(callsFile, file = paste0(outDir, "/Calls_merging_",approach, "_", "cutoff=", cutoff,"_species=", species, "_sex=", sex, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
                  }
                }
              }
              } else if ((length(condition)==2) == TRUE & "species_id" %in% condition[1:length(condition)] & "strain" %in% condition[1:length(condition)]){
                message("Merging will be done for condition c(species, strain)", "\n")
                if(is_empty(fileUser[[1]]$strain) == TRUE){
                  stop("The column strain is empty")
                }
                  
                  for (species in unique(fileUser[[1]]$species_id)) {
                    for (strain in unique(fileUser[[1]]$strain)) {
                      ## collect libraries from each species and from each strain
                      librariesSpecies <- fileUser[[1]]$output_directory[fileUser[[1]]$species_id == species & fileUser[[1]]$strain == strain]
                      if (is_empty(librariesSpecies) == TRUE){
                      } else {
                      ## collect files of individual call
                      AllFiles <- list.files(librariesSpecies, pattern="gene_level_abundance\\+calls.tsv", full.names=T, recursive = TRUE)
                      message("Using ", length(AllFiles), " libraries for the condition: c(", species ,",", strain, ")" ,"\n")
                      AllFiles <- lapply(AllFiles, read.delim)
                      ## perform the calls and export information of the merging
                      callsFile <- approachesMerging(AllFiles = AllFiles, approach = approach, cutoff = cutoff)
                      write.table(callsFile, file = paste0(outDir, "/Calls_merging_",approach, "_", "cutoff=", cutoff,"_species=", species, "_strain=", strain, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
                      }
                    }
                  }
                  } else if ((length(condition)==3) == TRUE & "species_id" %in% condition[1:length(condition)] & "anatEntity" %in% condition[1:length(condition)]& "devStage" %in% condition[1:length(condition)]){
                    message("Merging will be done for condition c(species, anatEntity, devStage)", "\n")
                    if(is_empty(fileUser[[1]]$anatEntity) == TRUE | is_empty(fileUser[[1]]$devStage) == TRUE){
                      stop("The column anatEntity or devStage are empty")
                    }
                    
                    for (species in unique(fileUser[[1]]$species_id)) {
                      for (anatEntity in unique(fileUser[[1]]$anatEntity)) {
                        for (devStage in unique(fileUser[[1]]$devStage)) {
                          
                          ## collect libraries from each species and from each anatomical entity and from each developmental stage
                          librariesSpecies <- fileUser[[1]]$output_directory[fileUser[[1]]$species_id == species & fileUser[[1]]$anatEntity == anatEntity & fileUser[[1]]$devStage == devStage]
                          if (is_empty(librariesSpecies) == TRUE){
                          } else {
                          ## collect files of individual call
                          AllFiles <- list.files(librariesSpecies, pattern="gene_level_abundance\\+calls.tsv", full.names=T, recursive = TRUE)
                          message("Using ", length(AllFiles), " libraries for the condition: c(", species ,",", anatEntity,",", devStage, ")" ,"\n")
                          AllFiles <- lapply(AllFiles, read.delim)
                          ## perform the calls and export information of the merging
                          callsFile <- approachesMerging(AllFiles = AllFiles, approach = approach, cutoff = cutoff)
                          write.table(callsFile, file = paste0(outDir, "/Calls_merging_",approach, "_", "cutoff=", cutoff,"_species=", species, "_anatEntity=", anatEntity, "_devStage=", devStage, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
                          }
                        }
                      }
                    }
                    } else if ((length(condition)==3) == TRUE & "species_id" %in% condition[1:length(condition)] & "anatEntity" %in% condition[1:length(condition)]& "sex" %in% condition[1:length(condition)]){
                      message("Merging will be done for condition c(species, anatEntity, sex)", "\n")
                      if(is_empty(fileUser[[1]]$anatEntity) == TRUE | is_empty(fileUser[[1]]$sex) == TRUE){
                        stop("The column anatEntity or sex are empty")
                      }
                        
                        for (species in unique(fileUser[[1]]$species_id)) {
                          for (anatEntity in unique(fileUser[[1]]$anatEntity)) {
                            for (sex in unique(fileUser[[1]]$sex)) {
                              
                              ## collect libraries from each species, from each anatomical entity and from sex
                              librariesSpecies <- fileUser[[1]]$output_directory[fileUser[[1]]$species_id == species & fileUser[[1]]$anatEntity == anatEntity & fileUser[[1]]$sex == sex]
                              if (is_empty(librariesSpecies) == TRUE){
                              } else {
                              ## collect files of individual call
                              AllFiles <- list.files(librariesSpecies, pattern="gene_level_abundance\\+calls.tsv", full.names=T, recursive = TRUE)
                              message("Using ", length(AllFiles), " libraries for the condition: c(", species ,",", anatEntity,",", sex, ")" ,"\n")
                              AllFiles <- lapply(AllFiles, read.delim)
                              ## perform the calls and export information of the merging
                              callsFile <- approachesMerging(AllFiles = AllFiles, approach = approach, cutoff = cutoff)
                              write.table(callsFile, file = paste0(outDir, "/Calls_merging_",approach, "_", "cutoff=", cutoff,"_species=", species, "_anatEntity=", anatEntity, "_sex=", sex, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
                              }
                            }
                          }
                        }
                        } else if ((length(condition)==3) == TRUE & "species_id" %in% condition[1:length(condition)] & "anatEntity" %in% condition[1:length(condition)]& "strain" %in% condition[1:length(condition)]){
                          message("Merging will be done for condition c(species, anatEntity, strain)", "\n")
                          if(is_empty(fileUser[[1]]$anatEntity) == TRUE | is_empty(fileUser[[1]]$strain) == TRUE){
                            stop("The column anatEntity or strain are empty")
                          }
                          
                          for (species in unique(fileUser[[1]]$species_id)) {
                            for (anatEntity in unique(fileUser[[1]]$anatEntity)) {
                              for (strain in unique(fileUser[[1]]$strain)) {
                                ## collect libraries from each species, from each anatomical entity and from strain
                                librariesSpecies <- fileUser[[1]]$output_directory[fileUser[[1]]$species_id == species & fileUser[[1]]$anatEntity == anatEntity & fileUser[[1]]$strain == strain]
                                if (is_empty(librariesSpecies) == TRUE){
                                } else {
                                ## collect files of individual call
                                AllFiles <- list.files(librariesSpecies, pattern="gene_level_abundance\\+calls.tsv", full.names=T, recursive = TRUE)
                                message("Using ", length(AllFiles), " libraries for the condition: c(", species ,",", anatEntity,",", strain, ")" ,"\n")
                                AllFiles <- lapply(AllFiles, read.delim)
                                ## perform the calls and export information of the merging
                                callsFile <- approachesMerging(AllFiles = AllFiles, approach = approach, cutoff = cutoff)
                                write.table(callsFile, file = paste0(outDir, "/Calls_merging_",approach, "_", "cutoff=", cutoff,"_species=", species, "_anatEntity=", anatEntity, "_strain", strain, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
                                }
                              }
                            }
                          }
                          } else if ((length(condition)==3) == TRUE & "species_id" %in% condition[1:length(condition)] & "devStage" %in% condition[1:length(condition)]& "sex" %in% condition[1:length(condition)]){
                            message("Merging will be done for condition c(species, devStage, sex)", "\n")
                            if(is_empty(fileUser[[1]]$devStage) == TRUE | is_empty(fileUser[[1]]$sex) == TRUE){
                              stop("The column devStage or sex are empty")
                            }
                            
                            for (species in unique(fileUser[[1]]$species_id)) {
                              for (devStage in unique(fileUser[[1]]$devStage)) {
                                for (sex in unique(fileUser[[1]]$sex)) {
                                  ## collect libraries from each species, from each developmental stage and from sex
                                  librariesSpecies <- fileUser[[1]]$output_directory[fileUser[[1]]$species_id == species & fileUser[[1]]$devStage == devStage & fileUser[[1]]$sex == sex]
                                  if (is_empty(librariesSpecies) == TRUE){
                                  } else {
                                  ## collect files of individual call
                                  AllFiles <- list.files(librariesSpecies, pattern="gene_level_abundance\\+calls.tsv", full.names=T, recursive = TRUE)
                                  message("Using ", length(AllFiles), " libraries for the condition: c(", species ,",", devStage,",", sex, ")" ,"\n")
                                  AllFiles <- lapply(AllFiles, read.delim)
                                  ## perform the calls and export information of the merging
                                  callsFile <- approachesMerging(AllFiles = AllFiles, approach = approach, cutoff = cutoff)
                                  write.table(callsFile, file = paste0(outDir, "/Calls_merging_",approach, "_", "cutoff=", cutoff,"_species=", species, "_devStage=", devStage, "_sex=", sex, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
                                  }
                                }
                              }
                            }
                            } else if ((length(condition)==3) == TRUE & "species_id" %in% condition[1:length(condition)] & "devStage" %in% condition[1:length(condition)]& "strain" %in% condition[1:length(condition)]){
                              message("Merging will be done for condition c(species, devStage, strain)", "\n")
                              if(is_empty(fileUser[[1]]$devStage) == TRUE | is_empty(fileUser[[1]]$strain) == TRUE){
                                stop("The column devStage or strain are empty")
                              }
                              
                              for (species in unique(fileUser[[1]]$species_id)) {
                                for (devStage in unique(fileUser[[1]]$devStage)) {
                                  for (strain in unique(fileUser[[1]]$strain)) {
                                    ## collect libraries from each species, from each developmental stage and from strain
                                    librariesSpecies <- fileUser[[1]]$output_directory[fileUser[[1]]$species_id == species & fileUser[[1]]$devStage == devStage & fileUser[[1]]$strain == strain]
                                    if (is_empty(librariesSpecies) == TRUE){
                                    } else {
                                    ## collect files of individual call
                                    AllFiles <- list.files(librariesSpecies, pattern="gene_level_abundance\\+calls.tsv", full.names=T, recursive = TRUE)
                                    message("Using ", length(AllFiles), " libraries for the condition: c(species_id, devStage, strain) = c(", species ,",", devStage,",", strain, ")" ,"\n")
                                    AllFiles <- lapply(AllFiles, read.delim)
                                    ## perform the calls and export information of the merging
                                    callsFile <- approachesMerging(AllFiles = AllFiles, approach = approach, cutoff = cutoff)
                                    write.table(callsFile, file = paste0(outDir, "/Calls_merging_",approach, "_", "cutoff=", cutoff,"_species=", species, "_devStage=", devStage, "_strain=", strain, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
                                    }
                                  }
                                }
                              }
                              } else if ((length(condition)==3) == TRUE & "species_id" %in% condition[1:length(condition)] & "sex" %in% condition[1:length(condition)]& "strain" %in% condition[1:length(condition)]){
                                message("Merging will be done for condition c(species, sex, strain)", "\n")
                                if(is_empty(fileUser[[1]]$sex) == TRUE | is_empty(fileUser[[1]]$strain) == TRUE){
                                  stop("The column sex or strain are empty")
                                }
                                
                                for (species in unique(fileUser[[1]]$species_id)) {
                                  for (sex in unique(fileUser[[1]]$sex)) {
                                    for (strain in unique(fileUser[[1]]$strain)) {
                                      ## collect libraries from each species, from each sex and from strain
                                      librariesSpecies <- fileUser[[1]]$output_directory[fileUser[[1]]$species_id == species & fileUser[[1]]$sex == sex & fileUser[[1]]$strain == strain]
                                      if (is_empty(librariesSpecies) == TRUE){
                                      } else {
                                      ## collect files of individual call
                                      AllFiles <- list.files(librariesSpecies, pattern="gene_level_abundance\\+calls.tsv", full.names=T, recursive = TRUE)
                                      message("Using ", length(AllFiles), " libraries for the condition: c(", species ,",", sex,",", strain, ")" ,"\n")
                                      AllFiles <- lapply(AllFiles, read.delim)
                                      ## perform the calls and export information of the merging
                                      callsFile <- approachesMerging(AllFiles = AllFiles, approach = approach, cutoff = cutoff)
                                      write.table(callsFile, file = paste0(outDir, "/Calls_merging_",approach, "_", "cutoff=", cutoff,"_species=", species, "_sex=", sex, "_strain=", strain, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
                                      }
                                    }
                                  }
                                }
                                } else if ((length(condition)==4) == TRUE & "species_id" %in% condition[1:length(condition)] & "anatEntity" %in% condition[1:length(condition)] & "devStage" %in% condition[1:length(condition)] & "sex" %in% condition[1:length(condition)]){
                                  message("Merging will be done for condition c(species, anatEntity, devStage, sex)", "\n")
                                  if(is_empty(fileUser[[1]]$anatEntity) == TRUE | is_empty(fileUser[[1]]$devStage) == TRUE | is_empty(fileUser[[1]]$sex) == TRUE){
                                    stop("The column anatEntity or devStage or sex are empty")
                                  }
                                  
                                  for (species in unique(fileUser[[1]]$species_id)) {
                                    for (anatEntity in unique(fileUser[[1]]$anatEntity)) {
                                      for (devStage in unique(fileUser[[1]]$devStage)) {
                                        for (sex in unique(fileUser[[1]]$sex)){
                                          ## collect libraries from each species, from each anatEntity, from each devStage and from sex
                                          librariesSpecies <- fileUser[[1]]$output_directory[fileUser[[1]]$species_id == species & fileUser[[1]]$anatEntity == anatEntity & fileUser[[1]]$devStage == devStage & fileUser[[1]]$sex == sex]
                                          if (is_empty(librariesSpecies) == TRUE){
                                          } else {
                                          ## collect files of individual call
                                          AllFiles <- list.files(librariesSpecies, pattern="gene_level_abundance\\+calls.tsv", full.names=T, recursive = TRUE)
                                          message("Using ", length(AllFiles), " libraries for the condition: c(", species ,",", anatEntity,",", devStage,",", sex, ")" ,"\n")
                                          AllFiles <- lapply(AllFiles, read.delim)
                                          ## perform the calls and export information of the merging
                                          callsFile <- approachesMerging(AllFiles = AllFiles, approach = approach, cutoff = cutoff)
                                          write.table(callsFile, file = paste0(outDir, "/Calls_merging_",approach, "_", "cutoff=", cutoff,"_species=", species, "_anatEntity=", anatEntity, "_devStage=", devStage, "_sex=", sex,".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
                                          }
                                        }
                                      }
                                    }
                                  }
                                } else if ((length(condition)==4) == TRUE & "species_id" %in% condition[1:length(condition)] & "anatEntity" %in% condition[1:length(condition)] & "devStage" %in% condition[1:length(condition)] & "strain" %in% condition[1:length(condition)]){
                                  message("Merging will be done for condition c(species, anatEntity, devStage, strain)", "\n")
                                  if(is_empty(fileUser[[1]]$anatEntity) == TRUE | is_empty(fileUser[[1]]$devStage) == TRUE | is_empty(fileUser[[1]]$strain) == TRUE){
                                    stop("The column anatEntity or devStage or strain are empty")
                                  }
                                  
                                  for (species in unique(fileUser[[1]]$species_id)) {
                                    for (anatEntity in unique(fileUser[[1]]$anatEntity)) {
                                      for (devStage in unique(fileUser[[1]]$devStage)) {
                                        for (strain in unique(fileUser[[1]]$strain)){
                                          ## collect libraries from each species, from each anatEntity, from each devStage and from each strain
                                          librariesSpecies <- fileUser[[1]]$output_directory[fileUser[[1]]$species_id == species & fileUser[[1]]$anatEntity == anatEntity & fileUser[[1]]$devStage == devStage & fileUser[[1]]$strain == strain]
                                          if (is_empty(librariesSpecies) == TRUE){
                                          } else {
                                            ## collect files of individual call
                                            AllFiles <- list.files(librariesSpecies, pattern="gene_level_abundance\\+calls.tsv", full.names=T, recursive = TRUE)
                                            message("Using ", length(AllFiles), " libraries for the condition: c(", species ,",", anatEntity,",", devStage,",", strain, ")" ,"\n")
                                            AllFiles <- lapply(AllFiles, read.delim)
                                            ## perform the calls and export information of the merging
                                            callsFile <- approachesMerging(AllFiles = AllFiles, approach = approach, cutoff = cutoff)
                                            write.table(callsFile, file = paste0(outDir, "/Calls_merging_",approach, "_", "cutoff=", cutoff,"_species=", species, "_anatEntity=", anatEntity, "_devStage=", devStage, "_strain=", strain,".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
                                          }
                                        }
                                      }
                                    }
                                  }
                                } else if ((length(condition)==4) == TRUE & "species_id" %in% condition[1:length(condition)] & "anatEntity" %in% condition[1:length(condition)] & "sex" %in% condition[1:length(condition)] & "strain" %in% condition[1:length(condition)]){
                                  message("Merging will be done for condition c(species, anatEntity, sex, strain)", "\n")
                                  if(is_empty(fileUser[[1]]$anatEntity) == TRUE | is_empty(fileUser[[1]]$sex) == TRUE | is_empty(fileUser[[1]]$strain) == TRUE){
                                    stop("The column anatEntity or sex or strain are empty")
                                  }
                                  
                                  for (species in unique(fileUser[[1]]$species_id)) {
                                    for (anatEntity in unique(fileUser[[1]]$anatEntity)) {
                                      for (sex in unique(fileUser[[1]]$sex)) {
                                        for (strain in unique(fileUser[[1]]$strain)){
                                          ## collect libraries from each species, from each anatEntity, from sex and from each strain
                                          librariesSpecies <- fileUser[[1]]$output_directory[fileUser[[1]]$species_id == species & fileUser[[1]]$anatEntity == anatEntity & fileUser[[1]]$sex == sex & fileUser[[1]]$strain == strain]
                                          if (is_empty(librariesSpecies) == TRUE){
                                          } else {
                                            ## collect files of individual call
                                            AllFiles <- list.files(librariesSpecies, pattern="gene_level_abundance\\+calls.tsv", full.names=T, recursive = TRUE)
                                            message("Using ", length(AllFiles), " libraries for the condition: c(", species ,",", anatEntity,",", sex,",", strain, ")" ,"\n")
                                            AllFiles <- lapply(AllFiles, read.delim)
                                            ## perform the calls and export information of the merging
                                            callsFile <- approachesMerging(AllFiles = AllFiles, approach = approach, cutoff = cutoff)
                                            write.table(callsFile, file = paste0(outDir, "/Calls_merging_",approach, "_", "cutoff=", cutoff,"_species=", species, "_anatEntity=", anatEntity, "_sex=", sex, "_strain=", strain,".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
                                          }
                                        }
                                      }
                                    }
                                  }
                                } else if ((length(condition)==4) == TRUE & "species_id" %in% condition[1:length(condition)] & "devStage" %in% condition[1:length(condition)] & "sex" %in% condition[1:length(condition)] & "strain" %in% condition[1:length(condition)]){
                                  message("Merging will be done for condition c(species, devStage, sex, strain)", "\n")
                                  if(is_empty(fileUser[[1]]$devStage) == TRUE | is_empty(fileUser[[1]]$sex) == TRUE | is_empty(fileUser[[1]]$strain) == TRUE){
                                    stop("The column devStage or sex or strain are empty")
                                  }
                                  
                                  for (species in unique(fileUser[[1]]$species_id)) {
                                    for (devStage in unique(fileUser[[1]]$devStage)) {
                                      for (sex in unique(fileUser[[1]]$sex)) {
                                        for (strain in unique(fileUser[[1]]$strain)){
                                          ## collect libraries from each species, from each devStage, from sex and from each strain
                                          librariesSpecies <- fileUser[[1]]$output_directory[fileUser[[1]]$species_id == species & fileUser[[1]]$devStage == devStage & fileUser[[1]]$sex == sex & fileUser[[1]]$strain == strain]
                                          if (is_empty(librariesSpecies) == TRUE){
                                          } else {
                                            ## collect files of individual call
                                            AllFiles <- list.files(librariesSpecies, pattern="gene_level_abundance\\+calls.tsv", full.names=T, recursive = TRUE)
                                            message("Using ", length(AllFiles), " libraries for the condition: c(", species ,",", devStage,",", sex,",", strain, ")" ,"\n")
                                            AllFiles <- lapply(AllFiles, read.delim)
                                            ## perform the calls and export information of the merging
                                            callsFile <- approachesMerging(AllFiles = AllFiles, approach = approach, cutoff = cutoff)
                                            write.table(callsFile, file = paste0(outDir, "/Calls_merging_",approach, "_", "cutoff=", cutoff,"_species=", species, "_devStage=", devStage, "_sex=", sex, "_strain=", strain,".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
                                          }
                                        }
                                      }
                                    }
                                  }
                                } else if ((length(condition)==5) == TRUE & "species_id" %in% condition[1:length(condition)] & "anatEntity" %in% condition[1:length(condition)] & "devStage" %in% condition[1:length(condition)] & "sex" %in% condition[1:length(condition)] & "strain" %in% condition[1:length(condition)]){
                                    message("Merging will be done for condition c(species, anatEntity, devStage, sex, strain)", "\n")
                                    if(is_empty(fileUser[[1]]$anatEntity) == TRUE | is_empty(fileUser[[1]]$devStage) == TRUE | is_empty(fileUser[[1]]$sex) == TRUE | is_empty(fileUser[[1]]$strain) == TRUE){
                                      stop("The column anatEntity or devStage or sex or strain are empty")
                                    }
                                    
                                    for (species in unique(fileUser[[1]]$species_id)) {
                                      for (anatEntity in unique(fileUser[[1]]$anatEntity)) {
                                        for (devStage in unique(fileUser[[1]]$devStage)) {
                                          for (sex in unique(fileUser[[1]]$sex)){
                                            for (strain in unique(fileUser[[1]]$strain)) {
                                              ## collect libraries from each species, from each anatEntity, from each devStage, from sex and from each strain
                                              librariesSpecies <- fileUser[[1]]$output_directory[fileUser[[1]]$species_id == species & fileUser[[1]]$anatEntity == anatEntity & fileUser[[1]]$devStage == devStage & fileUser[[1]]$sex == sex & fileUser[[1]]$strain == strain]
                                              if (is_empty(librariesSpecies) == TRUE){
                                              } else {
                                              ## collect files of individual call
                                              AllFiles <- list.files(librariesSpecies, pattern="gene_level_abundance\\+calls.tsv", full.names=T, recursive = TRUE)
                                              message("Using ", length(AllFiles), " libraries for the condition: c(", species ,",", anatEntity,",", devStage,",", sex, ",", strain, ")" ,"\n")
                                              AllFiles <- lapply(AllFiles, read.delim)
                                              ## perform the calls and export information of the merging
                                              callsFile <- approachesMerging(AllFiles = AllFiles, approach = approach, cutoff = cutoff)
                                              write.table(callsFile, file = paste0(outDir, "/Calls_merging_",approach, "_", "cutoff=", cutoff,"_species=", species, "_anatEntity=", anatEntity, "_devStage=", devStage, "_sex=", sex, "_strain=", strain,".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
                                              }
                                            }
                                          }
                                        }
                                      }
                                    }
                                  } else if ((length(condition)<=4) == TRUE & !"species_id" %in% condition[1:length(condition)]){
                                    stop("In the condition argument is missing species_id to perform the merging across libraries.", "\n")
                                  }
  }
}
