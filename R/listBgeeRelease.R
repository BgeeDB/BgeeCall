#' @title List Bgee releases available to use with bgeeCall package
#'
#' @description Returns information on available Bgee releases, the access URL for FTP, and the date of release
#'
#' @param release A character specifying a targeted release number (e.g., "14"). If not specified, all available releases are shown.
#'
#' @return A data frame with information on Bgee releases available to use with bgeeCall package.
#'
#' @examples{
#'  listBgeeRelease()
#' }
#'
#' @author Julien Roux, Julien Wollbrett
#' @export

# Function displaying the user a data frame describing all releases available for bgeeCall package
listBgeeRelease <- function(release=NULL){
  cat("Downloading release information from Bgee...\n")
  allReleases <- .getBgeeRelease()
  if (length(release) == 1){
    if (sum(allReleases$release == 1)) {
      cat(paste0("Only displaying information from targeted release ",
                 release, "\n"))
      allReleases <- allReleases[allReleases$release == release, ]
    } else {
      stop("ERROR: The specified release number is invalid or is not available for BgeeDB.")
    }
  }
  ## Only return the columns of interest to the user
  return(allReleases[, c("release","releaseDate", "FTPURL", "referenceIntergenicFastaURL", "minimumVersionBgeeCall", "messageToUsers")])
}

# Function returning a data frame describing all releases available for BgeeCall package
.getBgeeRelease <- function(removeFile=TRUE){
  ## query FTP to get file describing all releases
  releaseUrl <- 'ftp://ftp.bgee.org/bgeeCall_release.tsv'
  success <- try(download.file(releaseUrl, quiet = TRUE,
                               destfile=file.path(getwd(), 'release.tsv.tmp')))
  if (success == 0 & file.exists(file.path(getwd(), 'release.tsv.tmp'))) {
    file.rename(from=file.path(getwd(), 'release.tsv.tmp'),
                to=file.path(getwd(), 'release.tsv'))
    allReleases <- read.table("release.tsv", header=TRUE, sep="\t")
    if (removeFile == TRUE){
      file.remove(file.path(getwd(), 'release.tsv'))
    }
  } else {
    file.remove(file.path(getwd(), 'release.tsv.tmp'))
    stop("ERROR: File describing releases could not be downloaded from FTP.")
  }
  ## Keep release available with the current version of the package
  allAvailableReleases <- allReleases[  sapply(as.character(allReleases$minimumVersionBgeeCall), compareVersion, as.character(packageVersion("BgeeCall"))) <= 0,]
  return(allAvailableReleases)
}
