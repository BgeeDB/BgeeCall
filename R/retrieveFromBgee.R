#' @title List species having Bgee reference intergenic sequences
#'
#' @description Return information related to species having Bgee reference intergenic 
#' sequences available for the selected Bgee intergenic release:
#' - speciesId : the NCBI species ID of the species
#' - specieName : scientific species name
#' - numberOfLibraries : number of libraries used to generate these reference intergenic
#' sequences
#' If a BgeeMetadata object is provided this function retrieve the list of species using 
#' BgeeMetadata@intergenic_release.
#' If only a `release` is provided it will use it to retrieve the list of species.
#' If none of them are provided the default Bgee reference intergenic release will be used.
#'
#' @param myBgeeMetadata A Reference Class BgeeMetadata object
#' @param release A Bgee reference intergenic release name
#' @export
#' 
#' @return list all species having reference intergenic sequences available in the selected
#' release
#' 
#' @examples{
#' bgee <- new("BgeeMetadata")
#' list_bgee_ref_intergenic_species(myBgeeMetadata = bgee)
#' list_bgee_ref_intergenic_species(release = '0.2')
#' list_bgee_ref_intergenic_species()
#' }
#'
#' @author Julien Wollbrett
#' @export
#'

list_bgee_ref_intergenic_species <- function(myBgeeMetadata = NULL, release = NULL, speciesId = NULL) {
    non_bgee_releases <- c("community", "custom")
    species_info_file <- "species_info.tsv"
    if (!is.null(myBgeeMetadata)) {
        intergenic_release <- myBgeeMetadata@intergenic_release
        all_releases <- myBgeeMetadata@all_releases
    } else if (!(is.null(release))) {
        all_releases <- list_intergenic_release()
        intergenic_release <- release
    } else {
        stop("You should provide one BgeeMetadata object or a Bgee reference intergenic release name.")
    }
    if (intergenic_release %in% non_bgee_releases) {
        stop("`community` and `custom` are not Bgee reference intergenic releases.")
    } else if (!(intergenic_release %in% all_releases$release)) {
        stop("The selected reference intergenic release does not exist.")
    }
    url <- paste0(all_releases$FTPURL[all_releases$release == intergenic_release],
                  species_info_file)
    species <- read.table(file = url, header = TRUE, sep = "\t")
    return(species)
}

