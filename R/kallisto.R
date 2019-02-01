#' @title Create kallisto indexes.
#'
#' @description This function creates kallisto indexes. 2 indexes are created. One with default kmer value and one with kmer size of 21nt.
#'
#' @param myKallistoMetadata A Reference Class KallistoMetadata object.
#' @param myBgeeMetadata A Reference Class BgeeMetadata object.
#' @param myUserMetadata A Reference Class UserMetadata object.
#'
#' @author Julien Wollbrett.
#'
#' @export
#'
create_kallisto_index <- function (myKallistoMetadata, myBgeeMetadata, myUserMetadata) {

  # define needed path
  species_path <- get_species_path(myBgeeMetadata, myUserMetadata)
  kallisto_path <- get_tool_path(myKallistoMetadata, myBgeeMetadata, myUserMetadata)
  index_path <- get_tool_transcriptome_path(myKallistoMetadata, myBgeeMetadata, myUserMetadata)
  transcriptome_index_path <- file.path(index_path, myKallistoMetadata@index_file)
  transcriptome_k21_index_path <- file.path(index_path, myKallistoMetadata@k21_index_file)
  kallisto_exec <- getKallistoPath(myKallistoMetadata)
  transcriptome_path <- file.path(get_transcriptome_path(myBgeeMetadata, myUserMetadata), myKallistoMetadata@full_transcriptome_file)

  # test existence of files/directories
  if (!dir.exists(index_path)) {
    dir.create(index_path, recursive = T)
  }
  if ((myUserMetadata@reads_size >= myKallistoMetadata@read_size_kmer_threshold && file.exists(transcriptome_index_path)) || (myUserMetadata@reads_size < myKallistoMetadata@read_size_kmer_threshold && file.exists(transcriptome_k21_index_path))) {
    cat("Index file already exist. No need to create a new one.\n")
  } else {

    cat("Start generation of kallisto index files.\n")

    #create kallisto index with default kmer size
    if (myUserMetadata@reads_size >= myKallistoMetadata@read_size_kmer_threshold && !file.exists(transcriptome_index_path)) {
      kallisto_command <- paste0(kallisto_exec, " index -i ", transcriptome_index_path, " ", transcriptome_path)
      system(kallisto_command)
    }

    #create kallisto index with kmer size equal to 21
    if (myUserMetadata@reads_size < myKallistoMetadata@read_size_kmer_threshold && !file.exists(transcriptome_k21_index_path)) {
      kallisto_k21_command <- paste0(kallisto_exec, " index -k 21 -i ", transcriptome_k21_index_path, " ", transcriptome_path)
      system(kallisto_k21_command)
    }
    cat(paste0("kallisto index files have been succesfully created for species ", myUserMetadata@species_id, ".\n"))
  }
}

#' @title Run one kallisto abundance analyse
#'
#' @description Run kallisto and all preliminary steps if needed like :
#' - creation of transcriptome with intergenic (if needed)
#' - installation of kallisto (if needed)
#' - index creation (if needed)
#' - run kallisto quantification
#'
#' @param myKallistoMetadata A Reference Class KallistoMetadata object.
#' @param myBgeeMetadata A Reference Class BgeeMetadata object.
#' @param myUserMetadata A Reference Class UserMetadata object. This object has to be edited before running kallisto @seealso UserMetadata.R
#'
#' @author Julien Wollbrett.
#'
#' @export
#'
run_kallisto <- function (myKallistoMetadata, myBgeeMetadata, myUserMetadata) {

  # define path needed in this function
  kallisto_exec_path <- getKallistoPath(myKallistoMetadata)
  kallisto_index_dir <- get_tool_transcriptome_path(myKallistoMetadata, myBgeeMetadata, myUserMetadata)
  kallisto_index_path <- file.path(file.path(kallisto_index_dir, myKallistoMetadata@index_file))
  kallisto_output_path <- get_tool_output_path(myKallistoMetadata, myBgeeMetadata, myUserMetadata)
  output_abundance_file <- file.path(kallisto_output_path, myKallistoMetadata@abundance_file)
  
  # test if kallisto has to be run
  if (file.exists(output_abundance_file)) {
    cat(paste0("kallisto abundance file already exists for library ", basename(myUserMetadata@rnaseq_lib_path), ". No need to generate it again.\n"))
  } else {

    #create transcriptome containing bgee intergenic regions
    bgee_transcriptome <- merge_transcriptome_and_intergenic(myKallistoMetadata, myBgeeMetadata, myUserMetadata)

    # test if kallisto has to be installed
    if(myKallistoMetadata@download_kallisto) {
      if (is_kallisto_installed(myKallistoMetadata) == 1) {
        cat("It is the first time you try to use Kallisto downloaded from this package. Kallisto has to be downloaded.
        This version of Kallisto will only be used inside of this package.
        It will have no impact on your potential already installed version of Kallisto.\n")
        download_kallisto(myKallistoMetadata)
      }
    }

    #create index
    create_kallisto_index(myKallistoMetadata, myBgeeMetadata, myUserMetadata)

    # create output directory if not already existing
    if (!dir.exists(kallisto_output_path)) {
      dir.create(kallisto_output_path, recursive = T)
    }

    # if read size < 50nt use transcriptome index with small kmer size
    if (myUserMetadata@reads_size < 50) {
      kallisto_index_path <- file.path(file.path(kallisto_index_dir, myKallistoMetadata@k21_index_file))
    }

    # check library folder and test if _1 and _2 files are present
    fastq_files <- get_merged_fastq_file_names(myUserMetadata)
  
  
    kallisto_parameters <- myKallistoMetadata@single_end_parameters
    # if paired-end analyses
    if (grepl("_1.", strsplit(fastq_files, " ")[1])) {
      kallisto_parameters <- myKallistoMetadata@pair_end_parameters
    }
    kallisto_command <- paste(kallisto_exec_path, "quant -i", kallisto_index_path, "-o", kallisto_output_path, kallisto_parameters, fastq_files, sep = " ")
    cat(paste0("Will run kallisto using this command line : ", kallisto_command, "\n"))
    system(kallisto_command)
    if(myKallistoMetadata@ignoreTxVersion) {
      cat(paste0("remove transcript version info in ", myKallistoMetadata@abundance_file," ", myKallistoMetadata@tool_name, " abundance file.\n"))
      removeTxVersionFromAbundance(myKallistoMetadata, myBgeeMetadata, myUserMetadata)
    }
  }
}



#' @title Check if kallisto is already installed
#'
#' @description Check if kallisto is already installed.
#'
#' @param myKallistoMetadata A Reference Class KallistoMetadata object.
#'
#' @author Julien Wollbrett.
#'
#' @return a boolean
#'
#' @examples{
#' myKallistoMetadata <- new("KallistoMetadata")
#' logical <- is_kallisto_installed(myKallistoMetadata)
#' }
#' 
#' @noMd
#' @noRd
#' 

is_kallisto_installed <- function(myKallistoMetadata) {
  if ( file.exists(getKallistoPath(myKallistoMetadata ))) {
    return(0)
  } else {
    return (1)
  }
}

#' @title Download binary version of kallisto.
#'
#' @description Check your OS and download correct binary version of kallisto.
#'
#' @param myKallistoMetadata A Reference Class KallistoMetadata object.
#'
#' @author Julien Wollbrett.
#'
#' @export
#'
download_kallisto <- function(myKallistoMetadata) {
  if (dir.exists(myKallistoMetadata@kallisto_dir)) {
    message("kallisto directory already present. Kallisto do not need to be downloaded and installed again.")
  } else {
    dir.create(myKallistoMetadata@kallisto_dir, recursive = TRUE)
    os_version <- get_os()
    message(paste0("\nDownloading kallisto for ", os_version, "..."))

    # download .gz archive depending on the OS
    if(os_version == 'linux') {
      temp_path <- file.path(myKallistoMetadata@kallisto_dir, "temp.gz")
      success <- download.file(url = myKallistoMetadata@kallisto_linux_url, destfile=temp_path, mode='wb')
      untar(temp_path, exdir= getwd())
    } else if(os_version == 'osx') {
      temp_path <- file.path(myKallistoMetadata@kallisto_dir, "temp.gz")
      success <- download.file(url = myKallistoMetadata@kallisto_osx_url, destfile=temp_path, mode='wb')
      untar(temp_path, exdir=getwd())
    } else if (os_version == 'windows') {
      temp_path <- file.path(myKallistoMetadata@kallisto_dir, "temp.zip")
      success <- download.file(url = myKallistoMetadata@kallisto_windows_url, destfile=temp_path, mode='wb')
      unzip(temp_path, exdir=getwd())
    } else {
      stop("kallisto can not be downloaded on your computer. Only linux, OSX, and windows OS are compatible with this functionality.\n
         If you want to use this package please install your own version of Kallisto and run this command :\n
         \tmyKallistoMetadata@download_kallisto <- FALSE\n")
    }
    # delete downloaded archive of kallisto
    file.remove(temp_path)
    message("Kallisto has been succesfully installed.")
  }
}


