
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
