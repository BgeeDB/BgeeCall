#' @title List species having reference intergenic sequences created by the BgeeCall community
#'
#' @description Return information related to species having reference intergenic
#' sequences created by the BgeeCall community
#' - speciesId : the NCBI species ID of the species
#' - url : url to the reference intergenic fasta file
#' - numberOfLibraries : number of libraries used to generate these reference intergenic
#' sequences
#'
#' @author Julien Wollbrett
#'
#' @import jsonlite
#'
#' @export
#'
#' @return list all species having reference intergenic sequences created by the community
#'
#' @examples{
#' list_community_ref_intergenic_species()
#' }
#'
#'
list_community_ref_intergenic_species <- function() {
    community <-
        jsonlite::fromJSON(txt = "https://zenodo.org/api/records/?communities=bgee_intergenic")
    datasets <- NULL
    datasets_index <- 1
    for (records_index in seq_len(nrow(community))) {
        record <- community[records_index, ]
        keywords <- unlist(record$metadata$keywords)
        species_id <- unlist(strsplit(x = as.character(keywords[grep("speciesId", keywords)]), split = ":"))[2]
        annotation_version <- unlist(strsplit(x = as.character(keywords[grep("annotationVersion", keywords)]), split = ":"))[2]
        genome_version <- unlist(strsplit(x = as.character(keywords[grep("genomeVersion", keywords)]), split = ":"))[2]
        number_libraries <- unlist(strsplit(x = as.character(keywords[grep("numberOfLibraries", keywords)]), split = ":"))[2]
        kallisto_version <- unlist(strsplit(x = as.character(keywords[grep("kallistoVersion", keywords)]), split = ":"))[2]
        urls <- as.data.frame(record$files)
        # test presence of all mandatory metadata
        if (!(is.null(urls$links$download) ||
            is.null(species_id))) {
            # mandatory to have at least one file when uploading new dataset
            for (url_index in seq_len(nrow(urls))) {
                if (grepl(basename(urls[url_index, ]$links$download),
                    paste0("ref_intergenic.fa.gz"))) {
                    datasets$speciesId[datasets_index] <- as.character(species_id)
                    if (is.null(number_libraries)) {
                        datasets$numberOfLibraries[datasets_index] <- NA
                    } else {
                        datasets$numberOfLibraries[datasets_index] <-
                            as.character(number_libraries)
                    }
                    if (is.null(annotation_version)) {
                        datasets$annotationVersion[datasets_index] <- NA
                    } else {
                        datasets$annotationVersion[datasets_index] <-
                            as.character(annotation_version)
                    }
                    if (is.null(genome_version)) {
                        datasets$genomeVersion[datasets_index] <- NA
                    } else {
                        datasets$genomeVersion[datasets_index] <-
                            as.character(genome_version)
                    }
                    if (is.null(kallisto_version)) {
                        datasets$kallistoVersion[datasets_index] <- NA
                    } else {
                        datasets$kallistoVersion[datasets_index] <-
                            as.character(kallisto_version)
                    }
                    datasets$url[datasets_index] <-
                        as.character(urls[url_index, ]$links$download)
                    datasets_index <- datasets_index + 1
                }
            }
        }
    }
    return(as.data.frame(datasets))
}

retrieve_community_ref_intergenic_url <-
    function(speciesId, speciesDataSet = NULL) {
        species_dataset <- speciesDataSet
        if (is.null(species_dataset)) {
            species_dataset <- list_community_ref_intergenic_species()
        }
        if (!(as.character(speciesId) %in% as.character(species_dataset$speciesId))) {
            stop(
                "No reference intergenic sequences available for speciesId ",
                speciesId,
                " in the community release"
            )
        }
        file <-
            species_dataset$url[species_dataset$speciesId == as.character(speciesId)]
        return(as.character(file))
    }
