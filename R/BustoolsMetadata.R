#' @title BustoolsMetadata S4 class
#'
#' @description An S4 class that allows to get the required metadata to download, install and 
#' run bustools to generate a count matrix for droplet-based single-cell RNA-seq data.
#'
#' @slot download_bustools A logical allowing to use an already installed version
#' @slot bustools_windows_url URL to the binary of bustools for windows
#' @slot bustools_linux_url URL to the binary of bustools for linux
#' @slot bustools_osx_url URL to the binary of bustools for MacOS
#' @slot bustools_dir Name of the directory where bustools will be installed
#' @slot unix_bustools_name Name of the bustools executable in linux and macOS
#' @slot windows_bustools_name Name of the bustools executable in windows
#'
#' @exportClass BustoolsMetadata
BustoolsMetadata <- setClass(
    Class = "BustoolsMetadata",

    representation = representation(
        download_bustools = "logical",
        bustools_windows_url = "character",
        bustools_linux_url = "character",
        bustools_osx_url = "character",
        bustools_dir = "character",
        unix_bustools_name = "character",
        windows_bustools_name = "character",
        threads = "character"
    ),
    prototype = prototype(
        download_bustools = TRUE,
        bustools_windows_url = "https://github.com/BUStools/bustools/releases/download/v0.45.1/bustools_windows-master.zip",
        bustools_linux_url = "https://github.com/BUStools/bustools/releases/download/v0.45.1/bustools_linux-master.tar.gz",
        bustools_osx_url = "https://github.com/BUStools/bustools/releases/download/v0.45.1/bustools_mac-master.tar.gz",
        bustools_dir = "bustools",
        unix_bustools_name = "bustools",
        windows_bustools_name = "bustools.exe",
        threads = "4"
    ),
    validity = function(object) {
        if (!is.logical(object@download_bustools)) {
            return("download_bustools must be a logical value (TRUE or FALSE).")
        }
        return(TRUE)
    }
)
