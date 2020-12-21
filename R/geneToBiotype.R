#' @title Load gene to biotype
#'
#' @description Load the file containing the mapping between gene IDs and biotypes. 
#' This function also create the mapping file if it does not already exist
#' @param myAbundanceMetadata A descendant object of the Class myAbundanceMetadata.
#' @param myBgeeMetadata A Reference Class BgeeMetadata object.
#' @param myUserMetadata A Reference Class UserMetadata object. This object has 
#' to be edited before running kallisto @seealso UserMetadata.R
#'
#' @author Julien Wollbrett.
#'
#' @return The mapping between gene IDs and biotypes
#' 
#' @noMd
#' @noRd
#'
load_gene_to_biotype <- function(myAbundanceMetadata = new("AbundanceMetadata"), 
    myBgeeMetadata = new("BgeeMetadata"), myUserMetadata) {
    annotation_path <- get_annotation_path(myBgeeMetadata, 
                                           myUserMetadata)
    gene_to_biotype_file <- file.path(annotation_path, 
                                      myAbundanceMetadata@gene2biotype_file)
    create_gene_to_biotype(myAbundanceMetadata, myBgeeMetadata, myUserMetadata)
    return(read.table(gene_to_biotype_file, header = TRUE))
}

#' @title Create gene to biotype mapping
#'
#' @description Create a file containing the mapping between gene IDs and biotypes. 
#' @param myAbundanceMetadata A descendant object of the Class myAbundanceMetadata.
#' @param myBgeeMetadata A Reference Class BgeeMetadata object.
#' @param myUserMetadata A Reference Class UserMetadata object. This object has 
#' to be edited before running kallisto @seealso UserMetadata.R
#'
#' @author Julien Wollbrett.
#'
#' @return The mapping between gene IDs and biotypes
#'
#' @import rtracklayer
#' 
#' @return create the mapping between genes and biotypes
#' 
#' @noMd
#' @noRd
#'
create_gene_to_biotype <- function(myAbundanceMetadata = new("AbundanceMetadata"), 
                                 myBgeeMetadata = new("BgeeMetadata"), myUserMetadata) {   
    column_names <- c("id", "biotype", "type")
    annotation_path <- get_annotation_path(myBgeeMetadata, 
        myUserMetadata)
    gene_to_biotype_file <- file.path(annotation_path, 
        myAbundanceMetadata@gene2biotype_file)
    
    # check if file already exist
    if (!file.exists(gene_to_biotype_file)) {
        if (!dir.exists(annotation_path)) {
            dir.create(annotation_path, recursive = TRUE)
        }
        if(isTRUE(myUserMetadata@verbose)) {
            message("Generate file ", myAbundanceMetadata@gene2biotype_file, 
            ".\n")
        }
        # retrieve gene2biotype data frame from annotation
        # file
        gtf <- as.data.frame(myUserMetadata@annotation_object)
        gtf_gene <- gtf[gtf$source != "intergenic", 
            ]
        
        if (myUserMetadata@gtf_source == "ensembl"){
            gene_to_biotype <- as.data.frame(unique(cbind(gtf_gene$gene_id, 
                                                          gtf_gene$gene_biotype)))
        } else if (myUserMetadata@gtf_source == "gencode"){
            gene_to_biotype <- as.data.frame(unique(cbind(gtf_gene$gene_id, 
                                                          gtf_gene$gene_type)))
        } else {
            warning("The annotation file should be provided from ensembl or gencode.")
        }
        gene_to_biotype[, 3] <- "genic"
        names(gene_to_biotype) <- column_names
        
        # retrieve gene2biotype information from intergenic
        # fasta file
        bgee_intergenic_file <- retrieve_intergenic_path(myBgeeMetadata, 
            myUserMetadata)
        bgee_intergenic <- readDNAStringSet(bgee_intergenic_file)
        # intergenic ID correspond to part of the header
        # before the first space character
        intergenic_to_biotype <- as.data.frame(sub("^([^ ]+).*", 
            "\\1", names(bgee_intergenic)))
        intergenic_to_biotype[, 2] <- NA
        intergenic_to_biotype[, 3] <- "intergenic"
        names(intergenic_to_biotype) <- column_names
        
        # merge both data frame and write file
        gene_to_biotype <- rbind(gene_to_biotype, intergenic_to_biotype)
        write.table(gene_to_biotype, gene_to_biotype_file, 
            sep = "\t", row.names = FALSE, quote = FALSE)
    }
}
