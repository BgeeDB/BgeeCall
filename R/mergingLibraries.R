#' @title Checking user file 
#'
#' @description Check if the user provide a file with at least one
#' column (species_id) specified/fill-up to perform the merging of libraries
#'
#' @param userFile A file provided by the user
#' @param condition condition to considere for merging
#'
#' @author Sara Fonseca Costa
#'
#' @return User data frame
#' @return Vector name of columns detected
#' 
#' @import data.table
#' @import sjmisc
#' @import dplyr
#' 
#' @noMd
#' @noRd
#' 
checkUserFile <- function(userFile, condition){
  
  ## read and check columns of input file
  userFileCheck <- condition %in% names(userFile)
  
  if(!("species_id" %in% colnames(userFile)) | !("species_id" %in% condition)) {
    stop("the condition \"species_id\" should exist and be fill-up in order to do the merging.\n",
         "Please verify that \"species_id\" is present in your input file and that ",
         "the condition is choosed to merge.")
  }
  
  if(is_empty(userFile[,"species_id"]) ==TRUE) {
    stop("Some rows have empty value for the \"species_id\" column.\n",
         "\"species_id\" should always be provided.")
  }
  ## no column is detected in the input file
  if (FALSE %in% userFileCheck){
    stop("All provided merging condition should correspond to column name.\n",
         "condition to merge: ", paste(condition, ""), "\n",
         "column names :", paste(colnames(userFile), ""))
  }
}


#' @title Approaches used to merge/combine libraries
#'
#' @description For a set of libraries calculate the adjusted p-values using BH method
#' for each gene id OR calculate the inverse fdr correction to be applied to each gene ID 
#' in a set of libraries that use q-values.
#' The call of genes as present is done based on the assumption that at least one library 
#' have an adjusted p-value < cutoff or have a q-value < fdr_inverse.
#'
#' @param allFiles List with all libraries to be treated
#' @param approach Approach selected to merge libraries, BH or fdr_inverse
#' @param cutoff Threshold value to be applied to call expressed genes
#'
#' @author Sara Fonseca Costa
#' @author Alessandro Brandulas Cammarata
#'
#' @return data frame with minimum pValue or qValue detected (dependent of the approach selected) and the calls
#' 
#' @import data.table
#' @import dplyr
#' 
#' @noMd
#' @noRd
#' 
approachesMerging <- function(allFiles, approach, cutoff, weights=FALSE, weightValues=c()){
  
  librariesData <- try(do.call("cbind", allFiles), silent = TRUE) 
  ## select gene_id
  select_info <- librariesData[,1]
  
  if (approach != "fdr_inverse"){
    
    ## select all pValues column from all libraries
    select_pValue = librariesData[, grepl("^pValue", names(librariesData))]
    if(is_empty(select_pValue) == TRUE){
      stop("Select the appropriated method for your quantitative metric: BH, Mean or Median for p-values OR fdr_inverse for q-values", "\n")
    }
    ## calculate the p.adjusted values using a vector of original pValues for each gene_id
    if (approach == "BH"){
      collect_padjValues <- apply(select_pValue, 1, function (x) p.adjust(x[1:length(select_pValue)], method = "BH"))
      collect_padjValues <- as.data.frame(t(collect_padjValues))
      collect_padjValues$minimum_pValue <- do.call(pmin, c(collect_padjValues, list(na.rm = TRUE)))
    } else if (approach == "Mean"){
      collect_padjValues <- apply(select_pValue, 1, function (x) x[1:length(select_pValue)])
      collect_padjValues <- as.data.frame(t(collect_padjValues))
      if (weights == TRUE){
        collect_padjValues$minimum_pValue <- Pvalue_averaging(collect_padjValues, method = "mean", w_values = weightValues)
      } else {
        collect_padjValues$minimum_pValue <- Pvalue_averaging(collect_padjValues, method = "mean")
      }
    } else if (approach == "Median"){
      collect_padjValues <- apply(select_pValue, 1, function (x) x[1:length(select_pValue)])
      collect_padjValues <- as.data.frame(t(collect_padjValues))
      if (weights == TRUE){
        collect_padjValues$minimum_pValue <- Pvalue_averaging(collect_padjValues, method = "median", w_values = weightValues)
      } else {
      collect_padjValues$minimum_pValue <- Pvalue_averaging(collect_padjValues, method = "median")
      }
    }
    ## data frame with all information (gene_id + all adjusted pvalues of all libraries)
    allInfo <- data.frame(select_info, collect_padjValues)
    ## provide info just about id and minimum_pValue detected
    allInfo <- allInfo %>% dplyr::select("select_info", "minimum_pValue")
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
    allInfo <- allInfo %>% dplyr::select("select_info", "minimum_qValue")
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
#' @param cutoff Cutoff that should be applied to call Present/Absent genes
#' @param outDir Directory where the output files should be saved
#'
#' @author Sara Fonseca Costa
#' @author Alessandro Brandulas Cammarata
#' 
#' @export
#' 
#' @examples 
#' \dontrun{
#' callsMerging_species <- merging_libraries(userFile = 'PATH_USER_FILE', approach = 'BH', 
#' condition = 'species_id', cutoff = 0.05, outDir = 'PATH_OUTPUT')
#' callsMerging_species_sex <- merging_libraries(userFile = 'PATH_USER_FILE', approach = 'fdr_inverse', 
#' condition = c(species_id, sex), cutoff = 0.01, outDir = 'PATH_OUTPUT')
#' callsMerging_all <- merging_libraries(userFile = 'PATH_USER_FILE', approach = 'fdr_inverse', 
#' condition = c(species_id, anatEntity, devStage, sex, strain), cutoff = 0.05, outDir = 'PATH_OUTPUT')
#' }
#' 
#' @return A dataframe containing the minimum quantitative value (p-value or q-value) and 
#' the calls to each gene id for the referent condition.
#' 
#' 
merging_libraries <- function(userFile = NULL, approach = "BH", condition = "species_id", cutoff = 0.05, outDir = NULL, weights=FALSE) {

  ## check user input
  userFile <- read.table(file = userFile, header = TRUE, sep = "\t")
  checkUserFile(userFile = userFile, condition = condition)

  # retrieve unique condition combination to merge
  uniqueCondition <- unique(userFile[condition])
  
  # init variable used to keep only libraries corresponding to one condition
  filtered_libraries <- userFile
  
  for(currentRowIndex in seq(nrow(uniqueCondition))) {
    currentRow <- uniqueCondition[currentRowIndex,]
    # init name of the file where merged libraries will be stored depending on condition
    conditionFileName = ""
    
    # for each unique condition
    for(currentColumnIndex in seq(length(currentRow))) {
      # update name of file for each condition
      conditionFileName = paste0(conditionFileName, "_", condition[currentColumnIndex], "=", 
                                 currentRow[[currentColumnIndex]])
      filtered_libraries <- filtered_libraries[filtered_libraries[[condition[currentColumnIndex]]] == 
                                                 currentRow[[currentColumnIndex]],]
    }
    
    #retrieve all abundance files corresponding to condition
    if (approach == "BH" || approach == "Median" || approach == "Mean"){
      allFiles <- list.files(path = file.path(filtered_libraries$output_directory), pattern="gene_level_abundance\\+calls.tsv", 
                             full.names=T, recursive = TRUE)  
    } else {
      allFiles <- list.files(path = file.path(filtered_libraries$output_directory), pattern="gene_level_abundance\\+calls_qValue.tsv", 
                             full.names=T, recursive = TRUE)
    }
    message("Using ", length(allFiles), " libraries for condition: ", paste(condition,""), 
      " with values: ",paste(currentRow, ""))
    allFiles <- lapply(allFiles, read.delim)
    
    #merge calls based on condition
    if(weights == FALSE){
      callsFile <- approachesMerging(allFiles = allFiles, approach = approach, cutoff = cutoff, weights=FALSE)
    } else {
      weightValues = userFile["weights"]
      callsFile <- approachesMerging(allFiles = allFiles, approach = approach, cutoff = cutoff, weights=TRUE, weightValues = weightValues)
    }
    
    #write file with merged results
    write.table(callsFile, file = paste0(outDir, "/Calls_merging_",approach, "_", "cutoff=", cutoff, 
                  conditionFileName, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    
    #take all rows of the input file into consideration before filtering on an other set of condition
    filtered_libraries <- userFile
  }
}

#' @title Pvalue averaging
#'
#' @description Combined p-values for each gene Id using the mean or median averaging method
#'
#' @param pval_collect A data frame containing the p-values for each gene id of each library
#' @param method Method used to calculate the mean p-value
#' @param w_values A vector of weights to be used in the mean p-value calculation
#'
#' @author Alessandro Brandulas Cammarata
#' 
#' @export
#' 
#' @examples 
#' \dontrun{
  #' Pvalue_averaging(pval_collect = pval_collect, w_values = c(0.5, 0.5), method = "mean")
#' }
#' 
#' @return A dataframe containing the consensus p-value for each gene id
#' 
Pvalue_averaging <- function(pval_collect, w_values=c(), method="mean"){
  #Calculate the mean p-value for each gene_id
  corrected_mean_pval = c()
  #If no weights are provided, the mean p-value is calculated without weights
  if(!length(w_values)){
    for(row in 1:nrow(pval_collect)){
      if (method == "mean"){
      mean_pval = mean(as.double(pval_collect[row,]), na.rm = TRUE)
      } else if (method == "median"){
      mean_pval = median(as.double(pval_collect[row,]), na.rm = TRUE)
      }
      corrected_mean_pval = append(corrected_mean_pval, 2*mean_pval)
      corrected_mean_pval[corrected_mean_pval > 1] = 1
    }
  } else {
  #If weights are provided, the mean p-value is calculated with weights
    for(row in 1:nrow(pval_collect)){
      pvals = as.double(pval_collect[row,])
      if (method == "mean"){
        total_weight = sum(w_values)
        w_pvals = pvals * (w_values/total_weight)
        mean_pval = sum(w_pvals)
      } else if (method == "median"){
        mean_pval = weighted.median(pvals, w_values, na.rm = TRUE)
      }
      corrected_mean_pval = append(corrected_mean_pval, 2*mean_pval)
      corrected_mean_pval[corrected_mean_pval > 1] = 1
  }
  }
  #Return the corrected mean p-values for each gene_id
  return(corrected_mean_pval)
}
