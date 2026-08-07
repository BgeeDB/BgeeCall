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
