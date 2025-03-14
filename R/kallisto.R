#' @title Create kallisto indexes.
#'
#' @description This function creates kallisto indexes. 
#' Two indexes can be created depending on the reads size (see 
#' `myKallistoMetadata@read_size_kmer_threshold` and `UserMetadata@reads_size` 
#' for more information). One with default kmer value (31 nt) and one with 
#' kmer size of 15 nt. In order to generate.
#'
#' @param myKallistoMetadata A Reference Class KallistoMetadata object.
#' @param myBgeeMetadata A Reference Class BgeeMetadata object.
#' @param myUserMetadata A Reference Class UserMetadata object.
#' @param transcriptome_path path to the transcriptome fasta file. 
#' If no path is provided the default path created using BgeeCall
#' will be used. IMPORTANT : in BgeeCall the transcriptome used
#' to generate present/absent calls contains both intergenic sequences 
#' downloaded from Bgee and the reference transcriptome. If this 
#' function is run to generate present/absent then 'transcriptome_path' 
#' has to be empty
#'
#' @author Julien Wollbrett.
#'
#' @export
#' 
#' @examples 
#' \dontrun{
#' # first a transcriptome is needed. Here it is downloaded from AnnotationHub
#' library(AnnotationHub)
#' ah <- AnnotationHub()
#' ah_resources <- query(ah, c('Ensembl', 'Caenorhabditis elegans', '84'))
#' 
#' # kallisto can not deal with S4 objects. A Path to a transcriptome file is 
#' # required
#' transcriptome_object <- rtracklayer::import.2bit(ah_resources[['AH50453']])
#' transcriptome_path <- file.path(getwd(),'transcriptome.fa')
#' Biostrings::writeXStringSet(transcriptome_object, transcriptome_path)
#' 
#' 
#' 
#' # initialize objects needed to create destination folder
#' bgee <- new('BgeeMetadata')
#' user <- new('UserMetadata', species_id = '6239')
#' kallisto <- new('KallistoMetadata')
#' 
#' # generate transcriptome index
#' create_kallisto_index(kallisto, bgee, user, transcriptome_path)
#' }
#'
#' @return create kallisto index and save it on the hard drive
#'
create_kallisto_index <- function(myKallistoMetadata, 
    myBgeeMetadata, myUserMetadata, transcriptome_path = "") {
    
    # define needed path
    species_path <- get_species_path(myBgeeMetadata, 
        myUserMetadata)
    kallisto_path <- get_tool_path(myKallistoMetadata, 
        myBgeeMetadata, myUserMetadata)
    index_path <- get_tool_transcriptome_path(myKallistoMetadata, 
        myBgeeMetadata, myUserMetadata)
    transcriptome_index_path <- file.path(index_path, 
        myKallistoMetadata@index_file)
    transcriptome_k15_index_path <- file.path(index_path, 
        myKallistoMetadata@k15_index_file)
    kallisto_exec <- get_kallisto_program_path(myKallistoMetadata, 
        myUserMetadata)
    
    # If no transcriptome path is provided the default
    # path is used for index creation
    if (transcriptome_path == "") {
        transcriptome_path <- file.path(get_transcriptome_path(myBgeeMetadata, 
            myUserMetadata), myKallistoMetadata@full_transcriptome_file)
    }
    
    # test if kallisto has to be installed
    if (myKallistoMetadata@download_kallisto) {
        if (is_kallisto_installed(myKallistoMetadata, 
            myUserMetadata) == 1) {
            message("It is the first time you try to use Kallisto downloaded 
from this package. Kallisto has to be downloaded. This version of Kallisto 
will only be used inside of this package. It will have no impact on your 
potential already installed version of Kallisto.\n")
            download_kallisto(myKallistoMetadata, myUserMetadata)
        }
    } else {
        tryCatch({
            system2(command = "kallisto")
        }, error = function(w) {
  stop("kallisto is not installed. You should either 
    - automatically install a version of kallisto used only by this package (see vignette for more details)
    - install kallisto on your system following official website instructions (https://pachterlab.github.io/kallisto/download)")
        })
    }
    
    # test existence of files/directories
    if (!dir.exists(index_path)) {
        dir.create(index_path, recursive = TRUE)
    }
    need_to_generate_index <- TRUE
    if(!isTRUE(myKallistoMetadata@overwrite_index)) {
        if ((myUserMetadata@reads_size >= myKallistoMetadata@read_size_kmer_threshold && 
            file.exists(transcriptome_index_path)) || 
            (is.na(myUserMetadata@reads_size) == "TRUE" && 
             file.exists(transcriptome_index_path)) ||
            (myUserMetadata@reads_size < 
            myKallistoMetadata@read_size_kmer_threshold && 
            file.exists(transcriptome_k15_index_path))) {
            if(isTRUE(myUserMetadata@verbose)) {
                message("Index file already exist. No need to create a new one.\n")
            }
            need_to_generate_index <- FALSE
        } 
    } 
    if (need_to_generate_index) {
        if(isTRUE(myUserMetadata@verbose)) {
            message("Start generation of kallisto index files.\n")
        }
        
        # create kallisto index with default kmer size also when user not provide the read_size information
        if (myUserMetadata@reads_size >= myKallistoMetadata@read_size_kmer_threshold || is.na(myUserMetadata@reads_size) == "TRUE" && 
            !file.exists(transcriptome_index_path)) {
            kallisto_args <- paste0(" index -i ", transcriptome_index_path,
                " ", transcriptome_path)
            system2(command = kallisto_exec, args = kallisto_args)
        }
        
        # create kallisto index with kmer size equal to 15
        if (myUserMetadata@reads_size < myKallistoMetadata@read_size_kmer_threshold && 
            is.na(myUserMetadata@reads_size) == "FALSE" &&
            !file.exists(transcriptome_k15_index_path)) {
            kallisto_k15_args <- paste0(" index -k 15 -i ", 
                transcriptome_k15_index_path, " ", transcriptome_path)
            system2(command = kallisto_exec, args = kallisto_k15_args)
        }
        if(isTRUE(myUserMetadata@verbose)) {
            message("kallisto index files have been succesfully created for species ", 
                myUserMetadata@species_id, ".\n")
        }
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
#' @param myUserMetadata A Reference Class UserMetadata object. 
#' This object has to be edited before running kallisto @seealso 
#' UserMetadata.R
#' @param transcriptome_path path to the transcriptome fasta file. 
#' If no path is provided the default path created using BgeeCall
#' will be used. IMPORTANT : in BgeeCall the transcriptome used
#' to generate present/absent calls contains both intergenic sequences 
#' downloaded from Bgee and the reference transcriptome. 
#'
#' @author Julien Wollbrett.
#' 
#' @return NULL
#'
#' @export
#' 
#' @examples 
#' \dontrun{
#' # first a transcriptome is needed. Here it is downloaded from AnnotationHub
#' library(AnnotationHub)
#' ah <- AnnotationHub()
#' ah_resources <- query(ah, c('Ensembl', 'Caenorhabditis elegans', '84'))
# ,
#' 
#' # kallisto can not deal with S4 objects. Path to transcriptome file is 
#' # required
#' transcriptome_object <- rtracklayer::import.2bit(ah_resources[['AH50453']])
#' transcriptome_path <- file.path(getwd(),'transcriptome.fa')
#' Biostrings::writeXStringSet(transcriptome_object, transcriptome_path)
#' 
#' # initialize objects needed to create destination folder
#' bgee <- new('BgeeMetadata')
#' user <- new('UserMetadata', species_id = '6239')
#' user <- setRNASeqLibPath(user, system.file( 
#'                      'extdata', 'SRX099901_subset', 
#'                      package = 'BgeeCall'))
#' kallisto <- new('KallistoMetadata')
#' 
#' # generate transcriptome index
#' run_kallisto(kallisto, bgee, user, transcriptome_path)
#' }
#' 
#' @return create kallisto output files and save them on the hard drive
#'
run_kallisto <- function(myKallistoMetadata, myBgeeMetadata, 
    myUserMetadata, transcriptome_path = "") {
    
    # define path needed in this function
    kallisto_exec_path <- get_kallisto_program_path(myKallistoMetadata, 
        myUserMetadata)
    kallisto_index_dir <- get_tool_transcriptome_path(myKallistoMetadata, 
        myBgeeMetadata, myUserMetadata)
    kallisto_index_path <- file.path(file.path(kallisto_index_dir, 
        myKallistoMetadata@index_file))
    
    # use the standard output dir or the one defined by the user
    kallisto_output_path <- get_tool_output_path(myKallistoMetadata, 
                                                myBgeeMetadata, myUserMetadata)

    output_abundance_file <- file.path(kallisto_output_path, 
        myKallistoMetadata@abundance_file)
    
    # test if abundance file already exists
    need_to_generate_abundance = TRUE
    if (!isTRUE(myKallistoMetadata@overwrite_quant)) {
        if (file.exists(output_abundance_file)) {
            if(isTRUE(myUserMetadata@verbose)) {
                message("kallisto abundance file already exists for library ", 
                    basename(myUserMetadata@rnaseq_lib_path), 
                    ". No need to generate a new one.")
            }
            need_to_generate_abundance = FALSE
        }
    }
    if (need_to_generate_abundance) {
        # create transcriptome containing bgee intergenic regions
        if (transcriptome_path == "") {
            merge_transcriptome_and_intergenic(myKallistoMetadata, 
                myBgeeMetadata, myUserMetadata)
        }
    
        # test if kallisto has to be installed
        if (myKallistoMetadata@download_kallisto) {
            if (is_kallisto_installed(myKallistoMetadata, 
                myUserMetadata) == 1) {
                if(isTRUE(myUserMetadata@verbose)) {
                    message("It is the first time you try to use Kallisto downloaded 
from this package. Kallisto has to be downloaded. This version of Kallisto 
will only be used inside of this package. It will have no impact on your 
potential already installed version of Kallisto.")
                }
                download_kallisto(myKallistoMetadata, myUserMetadata)
            }
        }
    
        # create index
        create_kallisto_index(myKallistoMetadata, myBgeeMetadata, 
            myUserMetadata, transcriptome_path)
    
        # create output directory if not already existing
        if (!dir.exists(kallisto_output_path)) {
            dir.create(kallisto_output_path, recursive = TRUE)
        }
    
        # if read size < 50nt use transcriptome index with
        # small kmer size
        if (myUserMetadata@reads_size < myKallistoMetadata@read_size_kmer_threshold && 
            is.na(myUserMetadata@reads_size) == "FALSE") {
            kallisto_index_path <- file.path(file.path(kallisto_index_dir, 
                myKallistoMetadata@k15_index_file))
        }
 
        # check library folder and test if _1 and _2 files
        # are present
        fastq_files <- get_merged_fastq_file_names(myUserMetadata)
        kallisto_parameters <- myKallistoMetadata@single_end_parameters
        # if paired-end analyses
        if (grepl("_1.", lapply(strsplit(fastq_files, " "), 
            basename))[1]) {
            kallisto_parameters <- myKallistoMetadata@pair_end_parameters
        }
        kallisto_args <- paste("quant -i",  kallisto_index_path, 
            "-o", kallisto_output_path, kallisto_parameters, fastq_files, sep = " ")
      
        message("Will run kallisto using this command line : ", 
            paste(kallisto_exec_path, kallisto_args))
        
        # can not use system2 function for encrypted libraries as it needs to well manage piped 
        # commands
        if(is_encrypted_library(myUserMetadata)) {
          system(paste("echo \"", paste(kallisto_exec_path, kallisto_args), "\" | bash"))
        } else {
          system2(command = kallisto_exec_path, args = kallisto_args)
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
#' myKallistoMetadata <- new('KallistoMetadata')
#' logical <- is_kallisto_installed(myKallistoMetadata)
#' }
#' 
#' @noMd
#' @noRd
#' 

is_kallisto_installed <- function(myKallistoMetadata, 
    myUserMetadata) {
    if (file.exists(get_kallisto_program_path(myKallistoMetadata, 
        myUserMetadata))) {
        return(0)
    } else {
        return(1)
    }
}

#' @title Download binary version of kallisto.
#'
#' @description Check your OS and download correct binary 
#' version of kallisto.
#'
#' @param myKallistoMetadata A Reference Class KallistoMetadata object.
#' @param myUserMetadata A Reference Class UserMetadata object.
#'
#' @author Julien Wollbrett.
#' 
#' @examples {
#'   kallisto <- new('KallistoMetadata')
#'   user <- new('UserMetadata')
#'   download_kallisto(kallisto, user)
#' }
#' 
#' @export
#' 
#' @return save uncompressed executable of kallisto on the hard drive
#'
download_kallisto <- function(myKallistoMetadata, myUserMetadata) {
    kallisto_dir <- get_kallisto_dir_path(myKallistoMetadata, 
        myUserMetadata)
    if (dir.exists(kallisto_dir)) {
        if(isTRUE(myUserMetadata@verbose)) {
            message("kallisto directory already present. 
                Kallisto do not need to be downloaded and installed again.")
        }
    } else {
        dir.create(kallisto_dir, recursive = TRUE)
        os_version <- get_os()
        if(isTRUE(myUserMetadata@verbose)) {
            message(paste0("\nDownloading kallisto for ", os_version, "..."))
        }
        
        # download .gz archive depending on the OS
        if (os_version == "linux") {
            temp_path <- file.path(kallisto_dir, "temp.gz")
            success <- download.file(url = myKallistoMetadata@kallisto_linux_url, 
                destfile = temp_path, mode = "wb")
            untar(temp_path, exdir = myUserMetadata@working_path)
            
        } else if (os_version == "osx") {
            temp_path <- file.path(kallisto_dir, "temp.gz")
            success <- download.file(url = myKallistoMetadata@kallisto_osx_url, 
                destfile = temp_path, mode = "wb")
            untar(temp_path, exdir = myUserMetadata@working_path)
            
        } else if (os_version == "windows") {
            temp_path <- file.path(kallisto_dir, "temp.zip")
            success <- download.file(url = myKallistoMetadata@kallisto_windows_url, 
                destfile = temp_path, mode = "wb")
            unzip(temp_path, exdir = myUserMetadata@working_path)
            
        } else {
            stop("kallisto can not be downloaded on your computer. 
Only linux, OSX, and windows OS are compatible with this functionality.\n
If you want to use this package please install your own version of Kallisto ")
        }
        
        # delete downloaded archive of kallisto
        file.remove(temp_path)
        message("Kallisto has been succesfully installed.")
    }
}


