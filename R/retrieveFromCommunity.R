#' @title List species having reference intergenic sequences created by the BgeeCall community
#'
#' @description Return information related to species having reference intergenic 
#' sequences created by the BgeeCall community
#' - speciesId : the NCBI species ID of the species
#' - url : url to the reference intergenic fasta file
#' - numberOfLibraries : number of libraries used to generate these reference intergenic
#' sequences
#'
#' @import jsonlite
#' 
#' @return list all species having reference intergenic sequences created by the community
#' 
#' @examples{
#' list_community_ref_intergenic_species()
#' }
#'
#' @author Julien Wollbrett
#' @export
#'
#'
list_community_ref_intergenic_species <- function() {
    community <- fromJSON(txt = "https://sandbox.zenodo.org/api/records/?communities=test_community2")
    datasets <- NULL
    datasets_index <- 1
    for (records_index in 1:nrow(community)) {
        record <- community[records_index,]
       keywords <- unlist(record$metadata$keywords)
       speciesId <- unlist(strsplit( x= as.character(
           keywords[grep("speciesId", keywords)]), split = ":"))[2]
       type_intergenic <- unlist(strsplit(x= as.character(
           keywords[grep("type", keywords)]), split = ":"))[2]
       number_libraries <- unlist(strsplit(x= as.character(
           keywords[grep("numberOfLibraries", keywords)]), split = ":"))[2]
       urls <- as.data.frame(record$files)
       # test presence of all mandatory metadata
       if (!(is.na(urls) || is.na(speciesId) || is.na(type_intergenic) || 
             type_intergenic != "ref_intergenic")) {
            # mandatory to have at least one file when uploading new dataset
            for (url_index in 1:nrow(urls)) {
                if(grepl(basename(urls[url_index,]$links$download), 
                         paste0(speciesId,"_intergenic.fa.gz"))) {
                    datasets$speciesId[datasets_index] <- as.character(speciesId)
                    datasets$url[datasets_index] <- as.character(urls[url_index,]$links$download)
                    if(is.null(number_libraries)) {
                        datasets$numberOfLibraries[datasets_index] <- NA
                    } else {
                        datasets$numberOfLibraries[datasets_index] <- 
                            as.character(number_libraries)
                   }
                datasets_index <- datasets_index + 1
                }
            }
        }
    }
    return(as.data.frame(datasets))
}

retrieve_community_ref_intergenic_url <- function(speciesId, speciesDataSet = NULL) {
    species_dataset <- speciesDataSet
    if(is.null(species_dataset)) {
        species_dataset <- retrieve_community_ref_intergenic_species()
    }
    if(! (as.character(speciesId) %in% as.character(species_dataset$speciesId))) {
        stop("No reference intergenic sequences available for speciesId ", speciesId, ".")
    }
    file <- species_dataset$url[species_dataset$speciesId == as.character(speciesId)]
    return(as.character(file))
}
