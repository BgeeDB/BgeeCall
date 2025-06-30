#' @title Retrieve GTF files from ensembl, ensembl metazoa or ncbi 
#'
#' @description Downloads the GTF files of the given species for the specified ensembl release
#'
#' @param species_gtf list of genomes to download following from ensembl or path to the file containing the list of genomes
#' @param ensembl_release ensembl release from which we want to GTF files
#' @param ensembl_metazoa_release ensembl metazoa release from which we want to GTF files
#' @param outDir path to where we want to store the GTF files
#' @param from_file boolean of whether species_gtf is a file or a list of genomes
#'
#' @author Alessandro Brandulas Cammarata
#' 
#' @importFrom RCurl getURL
#' @importFrom readr read_tsv write_tsv parse_date
#' @importFrom stringr str_extract str_detect
#' @importFrom tools file_path_sans_ext
#' @importFrom curl curl_download
#' @importFrom jsonlite fromJSON
#' @importFrom dplyr select filter mutate arrange group_by summarise ungroup distinct
#' 
#' @noMd
#' @noRd
#' 
#' 
retrieve_gtf_files <- function(species_gtf=c("homo_sapiens/Homo_sapiens.GRCh38", "gallus_gallus/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b"), ensembl_release=112, ensembl_metazoa_release=59, outDir="./genomes/", from_file=FALSE, taxon_ids=FALSE) {

if (!dir.exists(outDir)) {
    dir.create(outDir, recursive = TRUE)
}

if (!dir.exists(outDir) || !file.access(outDir, 2) == 0) {
    stop(sprintf("Invalid or missing [%s]: %s\n", outDir, Sys.info()[["errno"]]))
}

if (from_file == TRUE) {
    tsv <- read_tsv(species_gtf, col_types = cols())
    species_gtf <- unique(tsv$genomeFilePath)
}

if (taxon_ids == TRUE) {
  species_gtf = find_genome_file_paths(species_gtf, ensembl_release, ensembl_metazoa_release, outDir) 
  species_gtf <- as.character(species_gtf)
}

working_path = getwd()
setwd(outDir)
print(sprintf("GTF files will be downloaded in %s", getwd()))

for (specie_gtf in sort(species_gtf)) {
  gtf_file <- basename(sprintf("%s.gtf.gz", specie_gtf))
  
  if (file.exists(gtf_file) && file.size(gtf_file) > 0) {
    next
  }
  
  ens_url <- sprintf("ftp://ftp.ensembl.org/pub/release-%d/gtf/%s.%d.gtf.gz", ensembl_release, specie_gtf, ensembl_release)
  metazoa_url <- sprintf("ftp://ftp.ensemblgenomes.org/pub/release-%d/metazoa/gtf/%s.%d.gtf.gz", ensembl_metazoa_release, specie_gtf, ensembl_metazoa_release)
  
  if (download_file(ens_url, gtf_file)) {
    cat(sprintf("Downloaded from Ensembl: %s\n", gtf_file))
  } else if (download_file(metazoa_url, gtf_file)) {
    cat(sprintf("Downloaded from Ensembl Metazoa: %s\n", gtf_file))
  } else if (str_detect(specie_gtf, "^\\w+/\\w+?_((GC[FA])_(\\d\\d\\d)(\\d\\d\\d)(\\d\\d\\d).*)$")) {
    match <- str_match(specie_gtf, "^\\w+/\\w+?_((GC[FA])_(\\d\\d\\d)(\\d\\d\\d)(\\d\\d\\d).*)$")
    ncbi_url <- sprintf("ftp://ftp.ncbi.nlm.nih.gov/genomes/all/%s/%s/%s/%s/%s/%s_genomic.gtf.gz", 
                        match[3], match[4], match[5], match[6], match[2], match[2])
    if (download_file(ncbi_url, gtf_file)) {
      cat(sprintf("Downloaded from NCBI: %s\n", gtf_file))
    } else {
      cat(sprintf("No GTF file found for %s in NCBI\n", specie_gtf))
    }
  } else {
    cat(sprintf("No GTF file found for %s in Ensembl %d, Ensembl Metazoa %d, or NCBI\n", 
                specie_gtf, ensembl_release, ensembl_metazoa_release))
  }
}
setwd(working_path)
}
#' @title Download from URL 
#'
#' @description Downloads files from URL and stores them in destination path
#'
#' @param url the url link to the object to download
#' @param dest the path where to store the data
#'
#' @author Alessandro Brandulas Cammarata
#'
#' @return Boolean of whether the retrival happened successfully
#' 
#' @importFrom RCurl getURL
#' 
#' @noMd
#' @noRd
#' 
#' 
download_file <- function(url, dest) {
  tryCatch({
    bin <- getURL(url)
    writeBin(bin, dest)
    return(TRUE)
  }, error = function(e) {
    message(sprintf("Error downloading %s: %s", url, e))
    return(FALSE)
  })
}

#' @title unzip selected file 
#'
#' @description unzips the provided file
#'
#' @param file path to the file to unzip
#'
#' @author Alessandro Brandulas Cammarata
#'
#' @return Boolean of whether the decompression happened successfully
#' 
#' 
#' @noMd
#' @noRd
#' 
#' 
unzip_file <- function(file) {
  tryCatch({
    system(sprintf("gunzip -f %s", file))
    return(TRUE)
  }, error = function(e) {
    message(sprintf("Error unzipping %s: %s", file, e))
    return(FALSE)
  })
}

#' @title Retrieve FASTA files from ensembl, ensembl metazoa or ncbi 
#'
#' @description Downloads the FASTA files of the given species for the specified ensembl release
#'
#' @param gtf_folder path to the folder containing the gtf files of the species for which we want the FASTA files
#' @param ensembl_release ensembl release from which we want to GTF files
#' @param ensembl_metazoa_release ensembl metazoa release from which we want to GTF files
#' @param outDir path to where we want to store the FASTA files
#'
#' @author Alessandro Brandulas Cammarata
#' 
#' @importFrom RCurl getURL
#' @importFrom readr read_tsv write_tsv parse_date
#' @importFrom stringr str_extract str_detect
#' @importFrom tools file_path_sans_ext
#' @importFrom curl curl_download
#' @importFrom jsonlite fromJSON
#' @importFrom dplyr select filter mutate arrange group_by summarise ungroup distinct
#' 
#' @noMd
#' @noRd
#' 
#' 
retrieve_fasta_files <- function(gtf_folder="./genomes/", ensembl_release=112, ensembl_metazoa_release=59, outDir="./genomes/") {

if (gtf_folder == "" || ensembl_release == 0 || ensembl_metazoa_release == 0 || outDir == "") {
  stop("Invalid or missing argument.")
}

if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = TRUE)
}

if (!dir.exists(outDir) || !file.access(outDir, 2) == 0) {
  stop(sprintf("Invalid or missing [%s]: %s\n", outDir, Sys.info()[["errno"]]))
}

working_path = getwd()
gtf_folder <- normalizePath(gtf_folder, winslash = "/", mustWork = FALSE)
setwd(outDir)
print(sprintf("FASTA files will be downloaded in %s", getwd()))

gtf_files <- list.files(gtf_folder, pattern = "\\.gtf\\.gz$", full.names = TRUE)
for (gtf in gtf_files) {
  if (file.info(gtf)$size == 0) {
    stop(sprintf("Problem with GTF file [%s]\n", gtf))
  }
  species_name <- tolower(strsplit(basename(gtf), ".", fixed = TRUE)[[1]][1])
  prefix <- sub("\\.gtf\\.gz$", "", basename(gtf))
  print(sprintf("Downloading %s", species_name))

  if (file.exists(sprintf("%s.genome.fa", prefix)) && file.info(sprintf("%s.genome.fa", prefix))$size > 0) {
    next
  }
  
  ens_url <- sprintf("ftp://ftp.ensembl.org/pub/release-%d/fasta/%s/dna/%s.dna.toplevel.fa.gz", ensembl_release, species_name, prefix)
  metazoa_url <- sprintf("ftp://ftp.ensemblgenomes.org/pub/release-%d/metazoa/fasta/%s/dna/%s.dna.toplevel.fa.gz", ensembl_metazoa_release, species_name, prefix)
  
  if (download_file(ens_url, sprintf("%s.genome.fa.gz", prefix))) {
    cat(sprintf("Downloaded from Ensembl: %s.genome.fa.gz\n", prefix))
    unzip_file(sprintf("%s.genome.fa.gz", prefix))
  } else if (download_file(metazoa_url, sprintf("%s.genome.fa.gz", prefix))) {
    cat(sprintf("Downloaded from Ensembl Metazoa: %s.genome.fa.gz\n", prefix))
    unzip_file(sprintf("%s.genome.fa.gz", prefix))
  } else if (str_detect(prefix, "^\\w+?_((GC[FA])_(\\d\\d\\d)(\\d\\d\\d)(\\d\\d\\d).*)$")) {
    match <- str_match(prefix, "^\\w+?_((GC[FA])_(\\d\\d\\d)(\\d\\d\\d)(\\d\\d\\d).*)$")
    ncbi_url <- sprintf("ftp://ftp.ncbi.nlm.nih.gov/genomes/all/%s/%s/%s/%s/%s/%s_genomic.fna.gz", 
                        match[3], match[4], match[5], match[6], match[2], match[2])
    if (download_file(ncbi_url, sprintf("%s.genome.fa.gz", prefix))) {
      cat(sprintf("Downloaded from NCBI: %s.genome.fa.gz\n", prefix))
      unzip_file(sprintf("%s.genome.fa.gz", prefix))
    } else {
      cat(sprintf("No genome file found for %s in NCBI\n", species_name))
    }
  } else {
    cat(sprintf("No genome file found for %s in Ensembl %d, Ensembl Metazoa %d, or NCBI\n", 
                species_name, ensembl_release, ensembl_metazoa_release))
  }
}
setwd(working_path)
}

#' @title retrieves FASTA and GTF
#'
#' @description retrieves both FASTA and GTF files for specified species
#'
#' @param species_gtf list of genomes to download following from ensembl or path to the file containing the list of genomes
#' @param ensembl_release ensembl release from which we want to GTF files
#' @param ensembl_metazoa_release ensembl metazoa release from which we want to GTF files
#' @param outDir path to where we want to store the GTF files
#' @param from_file boolean of whether species_gtf is a file or a list of genomes
#'
#' @author Alessandro Brandulas Cammarata
#' 
#' @importFrom RCurl getURL
#' @importFrom readr read_tsv write_tsv parse_date
#' @importFrom stringr str_extract str_detect
#' @importFrom tools file_path_sans_ext
#' @importFrom curl curl_download
#' @importFrom jsonlite fromJSON
#' @importFrom dplyr select filter mutate arrange group_by summarise ungroup distinct
#' 
#' @noMd
#' @noRd
#' 
#' 
retrieve_fasta_gtf <- function(species_gtf=c("homo_sapiens/Homo_sapiens.GRCh38", "gallus_gallus/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b"), ensembl_release=112, ensembl_metazoa_release=59, outDir="./genomes/", from_file=FALSE) {
    retrieve_gtf_files(species_gtf, ensembl_release, ensembl_metazoa_release, outDir, from_file)
    retrieve_fasta_files(gtf_folder=outDir, ensembl_release, ensembl_metazoa_release, outDir=outDir)
}

#' @title finds the ensembl FTP file path
#'
#' @description takes a taxon ID as input and retrieves the location of the genome and fasta files from the specified ensembl version
#'
#' @param species_gtf list of genomes to download following from ensembl or path to the file containing the list of genomes
#' @param ensembl_release ensembl release from which we want to GTF files
#' @param ensembl_metazoa_release ensembl metazoa release from which we want to GTF files
#' @param outDir path to where we want to store the GTF files
#'
#' @author Alessandro Brandulas Cammarata
#' 
#' @importFrom curl curl_download
#' @importFrom jsonlite fromJSON
#' @importFrom dplyr select filter mutate arrange group_by summarise ungroup distinct
#' @importFrom readr read_delim
#' 
#' @noMd
#' @noRd
#' 
#' 
find_genome_file_paths <- function(species_taxon_ids=c(9606, 9031), ensembl_release=112, ensembl_metazoa_release=59, outDir="./genomes/") {  
  # Helper function to download and process JSON data
  process_json <- function(url, output_filename) {
    output_file <- file.path(outDir, output_filename)
    if (!file.exists(output_file)) {
      message(paste("Downloading:", url))
      curl_download(url, destfile = output_file, quiet = TRUE)
    } else {
      message(paste("File already exists:", output_file))
    }
    tryCatch({
      fromJSON(output_file, flatten = TRUE)
    }, error = function(e) {
      warning(paste("Failed to process JSON from:", url, "Error:", e))
      return(NULL)
    })
  }
  
  # Helper function to find matching rows by taxonomy ID
  find_matching_rows <- function(species_metadata, taxon_id) {
    species_metadata %>%
      filter(organism.taxonomy_id == taxon_id) %>%
      select(assembly.assembly_name, organism.url_name, organism.name, organism.taxonomy_id, organism.strain) %>%
      rename(assembly_name = assembly.assembly_name,
             url_name = organism.url_name,
             organism_name = organism.name,
             taxonomy_id = organism.taxonomy_id,
             strain = organism.strain) %>%
      distinct()
  }
  
  # Helper function to process NCBI file
  process_ncbi <- function(ncbi_url, taxon_id) {
    output_file <- file.path(outDir, basename(ncbi_url))
    if (!file.exists(output_file)) {
      message(paste("Downloading:", ncbi_url))
      curl_download(ncbi_url, destfile = output_file, quiet = TRUE)
    } else {
      message(paste("File already exists:", output_file))
    }
    
    ncbi_data <- read_delim(output_file, delim = "\t", skip = 1, col_names = TRUE, show_col_types = FALSE, na = c("na", "NA", ""))
    
    # Rename the problematic column
    colnames(ncbi_data)[colnames(ncbi_data) == "#assembly_accession"] <- "assembly_accession"
    
    # Match the taxon_id with species_taxid and filter by refseq_category
    matching_rows <- ncbi_data %>%
      filter(species_taxid == taxon_id) %>%
      filter(refseq_category %in% c("reference genome", "representative genome")) %>%
      select(asm_name, assembly_accession, organism_name)
    
    # Return concatenated genomefilepath if a match is found
    if (nrow(matching_rows) > 0) {
      species_name_formatted <- gsub(" ", "_", tolower(matching_rows$organism_name))
      return(paste0(species_name_formatted, "/", species_name_formatted, "_", matching_rows$assembly_accession, "_", matching_rows$asm_name))
    } else {
      return(NA)  # Return NA if no matches
    }
  }
  
  # URLs for Ensembl Vertebrates, Metazoa, and NCBI
  url_ensembl <- paste0("https://ftp.ensembl.org/pub/release-", ensembl_release, "/species_metadata_EnsemblVertebrates.json")
  url_metazoa <- paste0("http://ftp.ensemblgenomes.org/pub/metazoa/release-", ensembl_metazoa_release, "/species_metadata_EnsemblMetazoa.json")
  url_ncbi <- "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"
  
  # Initialize a list to store the genome file paths
  genomefilepaths <- list()
  
  # Process the Vertebrates JSON (save as species_metadata_EnsemblVertebrates.json)
  species_metadata_vertebrates <- process_json(url_ensembl, "species_metadata_EnsemblVertebrates.json")
  
  # Process the Metazoa JSON (save as species_metadata_EnsemblMetazoa.json)
  species_metadata_metazoa <- process_json(url_metazoa, "species_metadata_EnsemblMetazoa.json")
  
  # Loop through each species_taxon_id in the list
  for (taxon_id in species_taxon_ids) {
    matching_rows <- data.frame()  # Initialize as empty
    
    if (!is.null(species_metadata_vertebrates)) {
      matching_rows <- find_matching_rows(species_metadata_vertebrates, taxon_id)
      
      if (nrow(matching_rows) > 1) {
        matching_rows <- matching_rows %>% filter(strain == "reference")
      }
    }
    
    if (nrow(matching_rows) == 0) {
      if (!is.null(species_metadata_metazoa)) {
        matching_rows <- find_matching_rows(species_metadata_metazoa, taxon_id)
        
        if (nrow(matching_rows) > 1) {
          matching_rows <- matching_rows %>% filter(strain == "reference")
        }
      }
    }
    
    # If still no match, process the NCBI file
    if (nrow(matching_rows) == 0) {
      genomefilepath <- process_ncbi(url_ncbi, taxon_id)
      genomefilepaths <- append(genomefilepaths, list(genomefilepath))  # Append as list element
    } else {
      # Create the concatenated string if matching rows are found
      genomefilepath <- paste0(matching_rows$organism_name[1], "/", matching_rows$url_name[1], ".", matching_rows$assembly_name[1])
      genomefilepaths <- append(genomefilepaths, list(genomefilepath))  # Append as list element
    }
  }
  
  return(genomefilepaths)
  
}

#' @title retrieves FASTA and GTF files from taxon ID
#'
#' @description retrieves both FASTA and GTF files for specified species
#'
#' @param species_gtf list of genomes to download following from ensembl or path to the file containing the list of genomes
#' @param ensembl_release ensembl release from which we want to GTF files
#' @param ensembl_metazoa_release ensembl metazoa release from which we want to GTF files
#' @param outDir path to where we want to store the GTF files
#'
#' @author Alessandro Brandulas Cammarata
#' 
#' @importFrom RCurl getURL
#' @importFrom readr read_tsv write_tsv parse_date
#' @importFrom stringr str_extract str_detect
#' @importFrom tools file_path_sans_ext
#' @importFrom curl curl_download
#' @importFrom jsonlite fromJSON
#' @importFrom dplyr select filter mutate arrange group_by summarise ungroup distinct
#' 
#' @noMd
#' @noRd
#' 
#' 
retrieve_fasta_gtf_from_taxonid <- function(taxon_id=c(9606, 9031), ensembl_release=112, ensembl_metazoa_release=59, outDir="./genomes/") {
    # Get the absolute path of outDir
    full_outDir <- normalizePath(outDir, mustWork = FALSE)
    
    # Print the full path
    print(paste0("Retrieving GTF and FASTA files from taxon ID at: ", full_outDir))
    
    # Call your retrieval functions
    retrieve_gtf_files(taxon_id, ensembl_release, ensembl_metazoa_release, outDir, taxon_ids=TRUE)
    retrieve_fasta_files(gtf_folder=outDir, ensembl_release, ensembl_metazoa_release, outDir=outDir)
}