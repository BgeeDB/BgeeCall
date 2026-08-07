#' @title Validate a droplet library run identifier
#'
#' @description Checks the validity of the run_id argument in dropletMetadata.
#' The run id is used as the output prefix for most output files,
#' so an empty or missing value silently drops the argument from the
#' command line.
#'
#' @param run_id Value of the `run_id` slot of a DropletMetadata object.
#'
#' @return `run_id` if valid, otherwise raises an error.
#'
#' @noMd
#' @noRd
#'
check_run_id <- function(run_id) {
    if (!is.character(run_id) || length(run_id) != 1L || is.na(run_id) ||
        !nzchar(run_id)) {
        stop("DropletMetadata@run_id must be a single non-empty character ",
            "string since it is used as the prefix of every bustools output",
             "file.")
    }
    if (grepl("[/\\\\]", run_id)) {
        stop("DropletMetadata@run_id must not contain a path separator. ",
            "Provided value : ", run_id)
    }
    invisible(run_id)
}

#' @title Sets the paths of the count matrix files written by bustools
#'
#' @description Generates the paths for the outputs of the bustools count step.
#'
#' @param dge_dir Directory containing the bustools outputs.
#' @param run_id Run Identifier of the droplet library
#'
#' @return A named list with the `prefix`, `mtx`, `genes` and `barcodes` paths.
#'
#' @noMd
#' @noRd
#'
bus_output_paths <- function(dge_dir, run_id) {
    check_run_id(run_id)
    prefix <- file.path(dge_dir, run_id)
    return(list(
        prefix = prefix,
        mtx = paste0(prefix, ".mtx"),
        genes = paste0(prefix, ".genes.txt"),
        barcodes = paste0(prefix, ".barcodes.txt")))
}

#' @title Ordering paired FASTQ files for kallisto bus
#'
#' @description Order the paired FASTQ files correctly for kallisto bus and
#' checks that there is the same number of r1 and r2 and that they exist.
#'
#' @param fastq_r1 Character vector of read 1 (barcode and UMI) files.
#' @param fastq_r2 Character vector of read 2 (biological sequence) files.
#' @param check_exists Logical. Check that every file exists (default TRUE).
#'
#' @return A character vector of interleaved file paths.
#'
#' @noMd
#' @noRd
#'
Order_fastq_r1_r2 <- function(fastq_r1, fastq_r2, check_exists = TRUE) {
    if (length(fastq_r1) == 0 && length(fastq_r2) == 0) {
        stop("both DropletMetadata@fastq_r1_path and ",
            "DropletMetadata@fastq_r2_path must contain at least one file.")
    }
    if (length(fastq_r1) != length(fastq_r2)) {
        stop("DropletMetadata@fastq_r1_path and DropletMetadata@fastq_r2_path ",
            "must have the same length (", length(fastq_r1), " and ",
            length(fastq_r2), " provided). kallisto bus expects one read 2 ",
            "file for each read 1 file.")
    }
    if (isTRUE(check_exists)) {
        all_files <- c(fastq_r1, fastq_r2)
        absent <- all_files[!file.exists(all_files)]
        if (length(absent) > 0) {
            stop("fastq file(s) not found : ", paste(absent, collapse = ", "))
        }
    }
    return(as.vector(rbind(fastq_r1, fastq_r2)))
}


#' @title Retrieve execution path or download Bustools
#' 
#' @description Checks if bustools binary already exists in the system or if absent and download flag
#' is set to TRUE, it downloads bustools from github.
#'
#' @param bustools_metadata An object of class BustoolsMetadata containing the required metadata to download
#' and run bustools.
#' @param myUserMetadata An object of class UserMetadata containing user-specific metadata such as working directory.
#' @return The absolute path to the bustools executable.
#' @export
get_bustools_path <- function(bustools_metadata, myUserMetadata) {

    # Get OS
    os_version <- get_os()
    is_windows <- os_version == "windows"

    # Find name of executable based on OS  
    bustools_exec_name <- ifelse(is_windows,
                                bustools_metadata@windows_bustools_name,
                                bustools_metadata@unix_bustools_name)
    
    system_bustools_path <- Sys.which(bustools_exec_name)
    if (system_bustools_path != "") {
        if(isTRUE(myUserMetadata@verbose)) {
            message("Found existing bustools executable in the system at: ", system_bustools_path,".
                    Using this version for the analysis.")
        }
        return(as.character(system_bustools_path))
    } else if (isFALSE(bustools_metadata@download_bustools)) {
        stop("Bustools executable not found in the system. Please set download_bustools to TRUE in the BustoolsMetadata
        object to allow the package to download it for you.")
    }
    bustools_local_dir <- file.path(myUserMetadata@working_path, bustools_metadata@bustools_dir)
    bustools_excutable_path <- file.path(bustools_local_dir, bustools_exec_name)
    if (file.exists(bustools_excutable_path)) {
        if(isTRUE(myUserMetadata@verbose)) {
            message("Found existing bustools executable in the local directory at: ", bustools_excutable_path,".
                    Using this version for the analysis.")
        }
        return(bustools_excutable_path)
    }

    # Create local directory if it does not exist
    if (!dir.exists(bustools_local_dir)) {
        dir.create(bustools_local_dir, recursive = TRUE)
    }

    # Determine the correct URL based on the OS
    bustools_url <- switch(os_version,
                            "windows" = bustools_metadata@bustools_windows_url, 
                            "linux" = bustools_metadata@bustools_linux_url,
                            "osx" = bustools_metadata@bustools_osx_url,
                            stop("Unsupported operating system: ", os_version, ". Supported OS are: windows, linux,
                                osx. If you want to use this package please install your own version of Kallisto"))
    
    # Reserving temporary file path for the downloaded archive
    temp_file <- tempfile()

    # Download the archive
    tryCatch({
        download.file(url = bustools_url, destfile = temp_file, mode = "wb", quiet = !isTRUE(myUserMetadata@verbose))},
            error = function(e) {
                stop("Failed to download bustools from ", bustools_url, ": ", e$message)
            }
        )
    # Uncompress bustools archive
    if (is_windows) {
        unzip(temp_file, exdir = myUserMetadata@working_path)
    } else {
        untar(temp_file, exdir = myUserMetadata@working_path)
    }
    unlink(temp_file)
    if (!file.exists(bustools_excutable_path)) {
        stop("Bustools executable not found after extraction. Please check the contents of the downloaded archive and ensure it contains the expected executable.")
    }

    # In linux, change permissions to make the file executable
    if (!is_windows) {
        Sys.chmod(bustools_excutable_path, mode = "0755")
    }
    if(isTRUE(myUserMetadata@verbose)) {
        message("Bustools successfully downloaded and installed at: ", bustools_excutable_path)
    }
    return(bustools_excutable_path)
}
