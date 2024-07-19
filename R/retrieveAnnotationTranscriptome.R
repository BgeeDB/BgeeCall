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
#' @import RCurl
#' @import readr
#' @import stringr
#' @import tools
#' 
#' @noMd
#' @noRd
#' 
#' 
retrieve_gtf_files <- function(species_gtf=c("homo_sapiens/Homo_sapiens.GRCh38", "gallus_gallus/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b"), ensembl_release=112, ensembl_metazoa_release=59, outDir="./genomes/", from_file=FALSE) {

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
#' @import RCurl
#' 
#' @noMd
#' @noRd
#' 
#' 
download_file <- function(url, dest) {
  tryCatch({
    bin <- getBinaryURL(url)
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
#' @import RCurl
#' @import readr
#' @import stringr
#' @import tools
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
#' @import RCurl
#' @import readr
#' @import stringr
#' @import tools
#' 
#' @noMd
#' @noRd
#' 
#' 
retrieve_fasta_gtf <- function(species_gtf=c("homo_sapiens/Homo_sapiens.GRCh38", "gallus_gallus/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b"), ensembl_release=112, ensembl_metazoa_release=59, outDir="./genomes/", from_file=FALSE) {
    retrieve_gtf_files(species_gtf, ensembl_release, ensembl_metazoa_release, outDir, from_file)
    retrieve_fasta_files(gtf_folder=outDir, ensembl_release, ensembl_metazoa_release, outDir=outDir)
}
