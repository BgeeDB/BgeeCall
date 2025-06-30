#' @title Download FASTQ files by Species and Tissue/Organ
#'
#' @description This function downloads FASTQ files from SRA for a specified species and tissue type.
#' It first queries the SRA database, filters the results by tissue type (if provided), and downloads
#' the specified number of libraries. If the SRAmetadb SQLite file is not already present, it will attempt
#' to download it using wget, which supports resuming.
#'
#' @param species_name Taxon ID of the species to download data for
#' @param tissue_keyword (Optional) Tissue or organ keyword to filter within the Attributes field
#' @param download_dir Directory where the FASTQ files should be downloaded
#' @param sratk_path Path to the SRA Toolkit installation directory
#' @param libraries_number Maximum number of libraries to download
#' @param sqlite_file Path to the SRAmetadb SQLite file
#' @noRd
download_sra_by_tissue <- function(species_name, tissue_keyword = NULL, download_dir, sratk_path, libraries_number, sqlite_file = "SRAmetadb.sqlite") {
    # Step 0: Check if the SQLite file exists; if not, download it
    if (!file.exists(sqlite_file)) {
        message("SRAmetadb SQLite file not found. Downloading with wget...")

        # Download the compressed SQLite file using wget with resume support
        download_command <- paste(
            "wget -c https://gbnci.cancer.gov/backup/SRAmetadb.sqlite.gz -O",
            paste0(sqlite_file, ".gz")
        )
        system(download_command)

        # Check if download was successful
        if (file.exists(paste0(sqlite_file, ".gz"))) {
            # Unzip the file
            message("Download completed. Decompressing the SQLite file...")
            system(paste("gunzip", paste0(sqlite_file, ".gz")))
            message("Decompression complete.")
        } else {
            stop("Failed to download SRAmetadb.sqlite file. Please check your internet connection or try downloading manually.")
        }
    }

    # Connect to the SRAmetadb SQLite file
    sra_con <- dbConnect(SQLite(), sqlite_file)
    on.exit(dbDisconnect(sra_con)) # Ensure the connection closes after function runs

    # Step 1: Query the SRA database for the species
    query <- paste0("SELECT * FROM sra WHERE taxon_id = '", species_name, "'")
    sra_data <- dbGetQuery(sra_con, query)

    # Step 2: Filter for the tissue/organ of interest within the Attributes field
    sra_data_filtered <- sra_data %>%
        filter(library_strategy == "RNA-Seq") %>%
        filter(if (!is.null(tissue_keyword)) grepl(tissue_keyword, attributes, ignore.case = TRUE) else TRUE)

    # Limit the number of libraries to the specified count
    sra_data_filtered <- head(sra_data_filtered, libraries_number)

    # Step 3: Download the FASTQ files for the filtered SRA IDs
    lapply(sra_data_filtered$experiment_name, function(sra_id) {
        run_accessions <- get_run_accessions_by_experiment(sra_id)
        library_path <- file.path(download_dir, sra_id)
        dir.create(library_path, showWarnings = FALSE)
        
        # Download each run for the experiment
        lapply(run_accessions, function(run_accession) {
            system(paste0(sratk_path, "/bin/prefetch --max-size 500G ", run_accession))
            system(paste0(sratk_path, "/bin/fasterq-dump --split-3 --outdir ", library_path, " ", run_accession))
        })
    })
}

#' @title Retrieve Run Accessions by Experiment Accession
#'
#' @description This function retrieves all run_accessions associated with a specific experiment accession.
#'
#' @param experiment_accession The experiment accession ID
#' @return A vector of run accessions for the experiment
#' @noRd
get_run_accessions_by_experiment <- function(experiment_accession) {
    query <- paste0("SELECT run_accession FROM sra WHERE experiment_accession = '", experiment_accession, "' AND run_accession IS NOT NULL")
    run_accessions <- dbGetQuery(sra_con, query)
    return(run_accessions$run_accession)
}
