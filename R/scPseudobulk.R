#' @title Aggregate a single-cell count matrix into pseudobulk counts per
#' cell type
#'
#' @description Sums the UMI counts of all the cells of each cell type,
#' turning the genes by cells matrix into a genes by cell types matrix.
#'
#' @param count_matrix Sparse matrix of UMI counts, genes in rows and cells
#' in columns, as returned by `read_bustools_matrix`.
#' @param celltype_annotation A data.frame holding a `barcode` and a
#' `celltype` column as returned by `resolve_celltype_annotation`.
#'
#' @param verbose Logical. Report the number of cells summed per cell type.
#'
#' @return A numeric matrix of summed UMI counts, genes in rows and cell
#' types in columns, cell types sorted alphabetically.
#'
#' @noMd
#' @noRd
#'
sc_pseudobulk_counts <- function(count_matrix, celltype_annotation,
    verbose = FALSE) {
    if (is.null(rownames(count_matrix)) || is.null(colnames(count_matrix))) {
        stop("the count matrix must have gene identifiers as row names and ",
            "cell barcodes as column names.")
    }
    absent <- setdiff(celltype_annotation$barcode, colnames(count_matrix))
    if (length(absent) > 0) {
        stop(length(absent), " annotated barcode(s) are missing from the ",
            "count matrix, for instance ", absent[1], ". The annotation ",
            "must be keyed on the matrix barcodes, see ",
            "resolve_celltype_annotation.")
    }

    cell_types <- sort(unique(as.character(celltype_annotation$celltype)))
    pseudobulk <- matrix(0, nrow = nrow(count_matrix),
        ncol = length(cell_types),
        dimnames = list(rownames(count_matrix), cell_types))
    for (cell_type in cell_types) {
        barcodes <- celltype_annotation$barcode[
            celltype_annotation$celltype == cell_type]

        pseudobulk[, cell_type] <- Matrix::rowSums(
            count_matrix[, barcodes, drop = FALSE])
        if (isTRUE(verbose)) {
            message("Summed ", length(barcodes), " cell(s) into the '",
                cell_type, "' pseudobulk sample.")
        }
    }
    return(pseudobulk)
}

#' @title Build the abundance tables of one pseudobulked cell type
#'
#' @description Turns the summed UMI counts of one cell type into the two
#' abundance data frames the calls generation works on, mirroring the two
#' tximport objects of the bulk pipeline :
#'
#' * `abundance` holds every feature. Its CPM are computed over all the
#' features, reference intergenic regions included, exactly as the bulk TPM
#' are computed over the merged transcriptome plus intergenic index. This
#' common scale is what makes genes comparable to the intergenic background,
#' so this is the table the statistics run on.
#' 
#' * `abundance_without_intergenic` holds the genic features only, with CPM
#' renormalised over them through `countToTpm`.
#' This is the table written to the calls output file.
#'
#' @param pseudobulk_counts Named numeric vector of summed UMI counts of one
#' cell type, one entry per feature (genes and intergenic regions).
#' @param biotype_mapping A data.frame holding the `id`, `biotype` and `type`
#' columns, as returned by `load_gene_to_biotype`.
#' @param verbose Logical. Report features missing from the mapping.
#'
#' @return A list of two data frames, `abundance` (columns `id`, `counts`,
#' `abundance`, `biotype`, `type`) and `abundance_without_intergenic` (same
#' columns, genic features only, CPM renormalised).
#'
#' @noMd
#' @noRd
#'
sc_abundance_tables <- function(pseudobulk_counts, biotype_mapping,
    verbose = FALSE) {
    if (is.null(names(pseudobulk_counts))) {
        stop("pseudobulk_counts must be named by feature identifier.")
    }
    total_counts <- sum(pseudobulk_counts)
    if (total_counts == 0) {
        stop("every pseudobulk count is zero, no abundance can be computed.")
    }

    abundance <- data.frame(id = names(pseudobulk_counts),
        counts = as.numeric(pseudobulk_counts),
        stringsAsFactors = FALSE)

    # Computes CPM over intergenic regions plus genic regions.
    abundance$abundance <- abundance$counts / total_counts * 1e+06

    # Join the biotype mapping table with the abundance table, dropping
    # features that are absent from the mapping.
    unmapped <- sum(!abundance$id %in% biotype_mapping$id)
    if (unmapped > 0 && isTRUE(verbose)) {
        message(unmapped, " feature(s) of the count matrix are absent from ",
            "the gene to biotype mapping and were dropped.")
    }
    abundance <- merge(abundance, biotype_mapping, by = "id", all = FALSE)

    if (!any(abundance$type == "intergenic")) {
        stop("no reference intergenic feature found in the pseudobulk ",
            "counts. The kallisto index must be built from the merged ",
            "transcriptome plus intergenic fasta and the transcript to ",
            "gene file must map the intergenic regions to themselves, ",
            "otherwise there is no background to call presence against.")
    }

    abundance_without_intergenic <-
        abundance[abundance$type == "genic", , drop = FALSE]
    
    # Renormalise CPM by excluding intergenic regions this time.
    abundance_without_intergenic$abundance <- countToTpm(
        abundance_without_intergenic$counts,
        rep_len(1, nrow(abundance_without_intergenic)))

    return(list(abundance = abundance,
        abundance_without_intergenic = abundance_without_intergenic))
}
