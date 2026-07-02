#' @title DropletMetadata s4 class
#'
#' @description An S4 class that contains all metadata needed to run the present/absent calls
#' for pseudobulked droplet-based single-cell RNA-seq data.
#'
#' @slot sequencing_technology Character string indicating the single-cell target technology 
#' (e.g., "10xV2", "10xV3", "DropSeq"). Essential for bustools processing (check kallisto kb --list for 
#' supported technologies).
#' @slot celltype_annotation A data.frame mapping unique cell barcodes to their 
#' corresponding biological groupings, such as cell types or clusters.
#' @slot count_matrix An optional highly sparse matrix (dgCMatrix) provided by the user 
#' if the computationally expensive fastq processing phase is bypassed.
#' @slot fastq_r1_path Character vector of paths to Read 1 FASTQ files (barcodes/UMIs).
#' @slot fastq_r2_path Character vector of paths to Read 2 FASTQ files (biological sequences).
#' @slot h5ad_path Optional character path to an AnnData (.h5ad) file containing the count matrix and metadata.
#' 
#' @exportClass DropletMetadata

DropletMetadata <- setClass(
    Class = "DropletMetadata",

    representation = representation(
        species_id = "character",
        run_id = "character",
        sequencing_technology = "character",
        celltype_annotation = "data.frame",
        count_matrix = "ANY",
        fastq_r1_path = "character",
        fastq_r2_path = "character",
        h5ad_path = "character",
        whitelist_path = "character"
    ),

    prototype = prototype(
        sequencing_technology = "10xv3",
        celltype_annotation = data.frame(),
        count_matrix = NULL,
        fastq_r1_path = character(0),
        fastq_r2_path = character(0),
        h5ad_path = character(0),
        whitelist_path = character(0)
    ),
    validity = function(object) {
        if (nrow(object@celltype_annotation) > 0) {
            if (!all(c("barcode", "celltype") %in% colnames(object@celltype_annotation))) {
                return("the celltype_annotation data.frame must contain a 'barcode' and 'celltype' column.")
            }   
        }
        if (is.null(object@count_matrix) && (length(object@h5ad_path) == 0) && 
            (length(object@fastq_r1_path) == 0 || length(object@fastq_r2_path) == 0)) {
                return("Please provide either a pre-computed count matrix, an h5ad file, or the paths to both R1 and R2 FASTQ files.")
            }
        return(TRUE)
    }
)