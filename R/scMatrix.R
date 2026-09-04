#' @title Read a bustools count matrix, transpose it (genesXbarcodes),
#' checks dimensions, checks for duplicates and returns a sparse matrix
#' of UMI counts
#'
#' @description Reads the three files written by `bustools count` and returns
#' a sparse matrix of UMI counts with genes as rows and cells as columns.
#' The matrix written by bustools has barcodes as rows and genes as columns,
#' so we transpose it. The dimensions are checked against the barcode
#' and gene files before transposing, and the gene and barcode identifiers
#' are checked for duplicates.
#'
#' @param mtx_path Path to the `<run_id>.mtx` file.
#' @param genes_path Path to the `<run_id>.genes.txt` file.
#' @param barcodes_path Path to the `<run_id>.barcodes.txt` file.
#' @param verbose Logical. Report the dimensions of the matrix read.
#'
#' @return A `dgCMatrix` of UMI counts, genes in rows and cells in columns,
#' with the gene identifiers as row names and the cell barcodes as column
#' names. Column compressed storage is used because every downstream step
#' subsets cells.
#'
#' @noMd
#' @noRd
#'
read_bustools_matrix <- function(mtx_path, genes_path, barcodes_path,
    verbose = FALSE) {

    expected_files <- c(mtx_path, genes_path, barcodes_path)
    absent <- expected_files[!file.exists(expected_files)]
    if (length(absent) > 0) {
        stop("bustools output file(s) not found : ",
            paste(absent, collapse = ", "))
    }

    counts <- Matrix::readMM(mtx_path)
    genes <- readLines(genes_path)
    barcodes <- readLines(barcodes_path)

    # Check that the dimensions of the count matrix match the number of barcodes
    # and genes provided.
    if (nrow(counts) != length(barcodes) || ncol(counts) != length(genes)) {
        stop("the matrix read from ", mtx_path, " is ", nrow(counts), " by ",
            ncol(counts), ", but ", length(barcodes), " barcodes and ",
            length(genes), " genes were provided. bustools count writes a ",
            "barcodes by genes matrix, so the barcode file should hold ",
            nrow(counts), " entries and the gene file ", ncol(counts),
            ". Check that the genes and barcodes files have not been swapped.")
    }

    # Check for duplicated names.
    if (anyDuplicated(genes) > 0) {
        stop("the gene file ", genes_path, " contains duplicated ",
            "identifiers, for instance ", genes[anyDuplicated(genes)], ".")
    }
    if (anyDuplicated(barcodes) > 0) {
        stop("the barcode file ", barcodes_path, " contains duplicated ",
            "barcodes, for instance ", barcodes[anyDuplicated(barcodes)], ".")
    }

    counts <- Matrix::t(counts)
    dimnames(counts) <- list(genes, barcodes)

    if (isTRUE(verbose)) {
        message("Read a count matrix of ", nrow(counts), " features by ",
            ncol(counts), " cells, holding ", length(counts@x),
            " non zero counts.")
    }
    return(methods::as(counts, "CsparseMatrix"))
}

#' @title Normalise cell barcodes names
#'
#' @description Remove non-nucleotide prefix in barcode names: a sample prefix
#' separated by an underscore, a colon or a dot (`sample1_AAAC...`) and the
#' lane or GEM well suffix (`AAAC...-1`). bustools emits the bare nucleotide
#' barcode, so the stripped form is the one the count matrix uses.
#'
#' @param barcodes Character vector of cell barcodes.
#' @param strip_prefix Logical. Remove everything up to the last underscore,
#' colon or dot.
#' @param strip_suffix Logical. Remove a trailing dash followed by digits.
#' @param to_upper Logical. Uppercase the result.
#'
#' @return A character vector holding the normalized barcodes.
#'
#' @noMd
#' @noRd
#'
harmonize_barcodes <- function(barcodes, strip_prefix = TRUE,
    strip_suffix = TRUE, to_upper = TRUE) {
    barcodes <- as.character(barcodes)
    if (isTRUE(strip_prefix)) {
        barcodes <- sub("^.*[_:.]", "", barcodes)
    }
    if (isTRUE(strip_suffix)) {
        barcodes <- sub("-[0-9]+$", "", barcodes)
    }
    if (isTRUE(to_upper)) {
        barcodes <- toupper(barcodes)
    }
    return(barcodes)
}

#' @title Match a cell type annotation to the barcodes of a count matrix
#'
#' @description Match the cell-type annotation to the barcodes of the count
#' matrix.
#'
#' @param celltype_annotation A data.frame holding a `barcode` and a
#' `celltype` column.
#' @param matrix_barcodes Character vector of barcodes as used by the count
#' matrix (its column names).
#' @param verbose Logical. Report how many barcodes were matched.
#'
#' @return A data.frame with a `barcode` column holding the matrix spelling
#' and a `celltype` column, one row per annotated barcode found in the matrix.
#'
#' @noMd
#' @noRd
#'
match_celltype_annotation <- function(celltype_annotation, matrix_barcodes,
    verbose = TRUE) {
    # Check that the annotation has the two required columns.
    if (!all(c("barcode", "celltype") %in% colnames(celltype_annotation))) {
        stop("the celltype_annotation data.frame must contain a 'barcode' ",
            "and a 'celltype' column, got : ",
            paste(colnames(celltype_annotation), collapse = ", "))
    }

    # Harmonize barcodes in both the annotation and the matrix.
    annotation_key <- harmonize_barcodes(celltype_annotation$barcode)
    matrix_key <- harmonize_barcodes(matrix_barcodes, strip_prefix = FALSE,
        strip_suffix = FALSE)

    # Look for duplicated barcodes, throw a warning and keep only the first one.
    duplicated_barcodes <- duplicated(annotation_key)
    if (any(duplicated_barcodes)) {
        warning(sum(duplicated_barcodes), " duplicated barcode(s) in the ",
            "cell type annotation. Only the first occurrence of each ",
            "barcode is kept.")
        celltype_annotation <-
            celltype_annotation[!duplicated_barcodes, , drop = FALSE]
        annotation_key <- annotation_key[!duplicated_barcodes]
    }

    matched_index <- match(annotation_key, matrix_key)
    keep <- !is.na(matched_index)

    if (!any(keep)) {
        # Diagnose the most common cause of a total mismatch : the annotation
        # follows the opposite strand convention. Diagnosed only, never fixed
        # silently.
        reverse_complement_hits <- 0
        tryCatch({
            sampled <- head(annotation_key, 1000)
            reverse_complement_hits <- sum(as.character(
                reverseComplement(DNAStringSet(sampled))) %in% matrix_key)
        }, error = function(e) NULL)
        stop("none of the ", length(annotation_key), " annotated barcodes ",
            "match the ", length(matrix_key), " barcodes of the count ",
            "matrix.\n  example annotated barcode (harmonised) : ",
            annotation_key[1], "\n  example matrix barcode : ",
            matrix_key[1], "\n",
            if (reverse_complement_hits > 0) {
                paste0("  NOTE : ", reverse_complement_hits,
                    " of the first ", length(head(annotation_key, 1000)),
                    " annotated barcodes match the matrix after reverse ",
                    "complementing. The annotation likely follows the ",
                    "opposite strand convention, please reverse complement ",
                    "its barcodes before use.\n")
            },
            "  Check the sequencing technology and the barcode convention ",
            "of the annotation.")
    }

    # Report the number of matched barcodes and its fraction.
    matched_fraction <- sum(keep) / length(annotation_key)
    if (isTRUE(verbose)) {
        message("Matched ", sum(keep), " of ", length(annotation_key),
            " annotated barcodes to the ", length(matrix_key),
            " barcodes of the count matrix.")
    }
    if (matched_fraction < 0.5) {
        warning("only ", round(100 * matched_fraction, 1), "% of the ",
            "annotated barcodes were found in the count matrix.")
    }

    return(data.frame(
        barcode = matrix_barcodes[matched_index[keep]],
        celltype = as.character(celltype_annotation$celltype[keep]),
        stringsAsFactors = FALSE))
}

#' @title Associate barcode to cell type annotation for pseudobulking
#'
#' @description Returns a data.frame associating the barcodes of the count
#' matrix to the cell type annotation provided in the DropletMetadata object.
#'
#' @param droplet_metadata A DropletMetadata object.
#' @param count_matrix The count matrix, genes in rows and cells in columns,
#' as returned by `read_bustools_matrix`.
#' @param min_umi_per_barcode Numeric. When no annotation is provided, only
#' barcodes holding at least this many UMI counts are pooled. The default of
#' 0 keeps every barcode.
#' @param verbose Logical. Report how many barcodes were matched or pooled.
#'
#' @return A data.frame with a `barcode` and a `celltype` column, keyed on
#' the matrix spelling of the barcodes.
#'
#' @noMd
#' @noRd
#'
resolve_celltype_annotation <- function(droplet_metadata, count_matrix,
    min_umi_per_barcode = 0, verbose = TRUE) {
    matrix_barcodes <- colnames(count_matrix)
    if (is.null(matrix_barcodes)) {
        stop("the count matrix has no column names, cell barcodes are ",
            "expected as column names.")
    }

    if (nrow(droplet_metadata@celltype_annotation) > 0) {
        return(match_celltype_annotation(
            droplet_metadata@celltype_annotation, matrix_barcodes,
            verbose = verbose))
    }

    keep <- rep(TRUE, length(matrix_barcodes))
    if (min_umi_per_barcode > 0) {
        keep <- Matrix::colSums(count_matrix) >= min_umi_per_barcode
        if (!any(keep)) {
            stop("no barcode holds at least ", min_umi_per_barcode,
                " UMI counts.")
        }
    }
    warning("no cell type annotation was provided : the ", sum(keep),
        " barcode(s) of the count matrix are pooled into a single ",
        "pseudobulk sample named 'all_cells'. No cell calling or empty ",
        "droplet filtering was performed, so this pool also holds the ",
        "ambient RNA of the empty droplets. Provide a celltype_annotation ",
        "(or a min_umi_per_barcode threshold) to generate biologically ",
        "meaningful calls per cell type.")
    return(data.frame(
        barcode = matrix_barcodes[keep],
        celltype = "all_cells",
        stringsAsFactors = FALSE))
}
