#' @title Gather statistical information
#'
#' @description Collect the statistics provided by the gene_cutoff_info_file from each individual library,
#' in order to generate a global summary file.
#' @param userFile A data frame containing all information of each library
#' @param outDir Output directory where the generated file should be saved
#' @return A tsv file
#'
#' @author Sara Fonseca Costa
#' 
#' @export
#' 
#' @noMd
#' @noRd

get_summary_stats <- function(userFile, outDir){
  
  userFile <- read.table(userFile, header=TRUE, sep="\t")
  
  ## create the output summary file
  summaryStats <- file.path(outDir, "summary_Stats_All_Libraries.tsv")
  if (!file.exists(summaryStats)){
    file.create(summaryStats)
  } 
  
  collectInfo <- c()
  for (pathLib in unique(userFile$output_directory)) {
    ## Note: multiple info files can exist per library (we have at the moment 3 different approaches to call expressed genes)
    files <- list.files(path = pathLib, pattern = "^gene_cutoff_info_file")
    
    for (library in unique(files)) {
      
      readFiles <- read.table(file.path(pathLib, library), header = FALSE)
      selectStats <- t(readFiles[1:11, ])[2,]
      approach <- readFiles[12,1]
      cutoff <- readFiles[12,2]
      
      if(approach == "pValueCutoff"){
        mean <- readFiles[13,2]
        standard_deviation <- readFiles[14,2]
      } else {
        mean <- NA
        standard_deviation <- NA
      }
      selectStats <- as.data.frame(t(selectStats))
      selectStats$approach <- approach 
      selectStats$cutoff <- cutoff
      selectStats$mean <- mean
      selectStats$standard_deviation <- standard_deviation
      collectInfo <- rbind(collectInfo, selectStats)
    }
  }
  ## write all information of all libraries to the summary file
  colnames(collectInfo) <- c("libraryId", "cutoffTPM", "proportionGenicPresent", "numberGenicPresent", "numberGenic", 
                             "proportionCodingPresent", "numberPresentCoding", "numberCoding", "proportionIntergenicPresent", 
                             "numberIntergenicPresent", "numberIntergenic", "approach", "cutoff", "mean", "standard_deviation")
  write.table(collectInfo, file = summaryStats, sep="\t", row.names = FALSE, quote = FALSE)
}
