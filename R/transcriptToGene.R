#' @title Generate transcript to gene mapping for intergenic
#'
#' @description Generate transcript to gene mapping for intergenic regions as 
#' used by tximport. Gene and transcript columns are identical.
#'
#' @param myBgeeMetadata A Reference Class BgeeMetadata object.
#' @param myUserMetadata A Reference Class UserMetadata object.
#'
#' @author Julien Wollbrett
#'
#' @return transcript to gene mapping for intergenic regions
#'
#' @noMd
#' @noRd
#' 
#' @examples { 
#' user <- new("UserMetadata", species_id = "6239")
#' bgee <- new("BgeeMetadata", intergenic_release = "0.1")
#' intergenic_tx2gene(bgee, user)
#' }
intergenic_tx2gene <- function(myBgeeMetadata, myUserMetadata) {
    species_path <- get_species_path(myBgeeMetadata, myUserMetadata)
    bgee_intergenic_file <- file.path(species_path, 
                                      myBgeeMetadata@fasta_intergenic_name)
    if (!file.exists(bgee_intergenic_file)) {
        if (!dir.exists(species_path)) {
            dir.create(species_path, recursive = TRUE)
        }
        download_fasta_intergenic(myBgeeMetadata, myUserMetadata, 
                                  bgee_intergenic_file)
    }
    bgee_intergenic <- readDNAStringSet(bgee_intergenic_file)
    #intergenic ID correspond to part of the header before the first space character
    all_transcripts <- as.data.frame(sub("^([^ ]+).*", "\\1", 
                                         names(bgee_intergenic)))
    # Create a second column idantical to the first one. 
    # Then each intergenic region will be consider as a gene (of the same name) 
    # for tximport
    all_transcripts[,2] <- all_transcripts[,1]
    names(all_transcripts) <- c("TXNAME", "GENEID")
    return(all_transcripts)
}


#' @title Create TxDb annotation
#'
#' @description Create TxDb annotation from gtf or gff3 annotations
#'
#' @param myAbundanceMetadata A descendant object of the Class 
#' myAbundanceMetadata.
#' @param myUserMetadata A Reference Class UserMetadata object.
#'
#' @author Julien Wollbrett
#'
#' @return TxDb annotation
#'
#' @import GenomicFeatures
#'
#' @noMd
#' @noRd
#'
create_TxDb <- function(myAbundanceMetadata, myUserMetadata) {
    # create txdb from GRanges Object
    txdb <- makeTxDbFromGRanges(myUserMetadata@annotation_object, 
                                taxonomyId = as.numeric(myUserMetadata@species_id))
    return(txdb)
}


#' @title Create transcript to gene mapping file
#'
#' @description Create transcript to gene mapping file as used by tximport. 
#' The file contains both genic and intergenic regions.
#'
#' @param myAbundanceMetadata A descendant object of the Class 
#' myAbundanceMetadata.
#' @param myBgeeMetadata A Reference Class BgeeMetadata object.
#' @param myUserMetadata A Reference Class UserMetadata object.
#'
#' @author Julien Wollbrett
#'
#' @export
#' 
#' @return path to the tx2gene file
#'
create_tx2gene <- function(myAbundanceMetadata, myBgeeMetadata, myUserMetadata) {
    # create tx2gene from annotations
    annotation_path <- get_annotation_path(myBgeeMetadata, myUserMetadata)
    tx2gene_file <- myAbundanceMetadata@tx2gene_file
    if(myAbundanceMetadata@ignoreTxVersion == TRUE) {
        tx2gene_file <- myAbundanceMetadata@tx2gene_file_without_version
    }
    tx2gene_path <- file.path(annotation_path, tx2gene_file)
    if (!file.exists(tx2gene_path)) {
        cat(paste0("Generate file ", tx2gene_file, ".\n"))
        if(!dir.exists(annotation_path)) {
            dir.create(annotation_path, recursive = TRUE)
        }
        txdb <- create_TxDb(myAbundanceMetadata, myUserMetadata)
        k <- biomaRt::keys(txdb, keytype = "TXNAME")
        tx2gene <- as.data.frame(biomaRt::select(txdb, k, "GENEID", "TXNAME"))
        intergenic_tx2gene <- intergenic_tx2gene(myBgeeMetadata, myUserMetadata)
        tx2gene <- rbind(tx2gene, intergenic_tx2gene)
        # Remove the transcript version that can be present in transcript id of 
        # gtf files
        if(myAbundanceMetadata@ignoreTxVersion) {
            cat(paste0("remove transcript version info in ", 
                       tx2gene_file, " file.\n"))
            tx2gene$TXNAME <- gsub(pattern = "\\..*", "", tx2gene$TXNAME )
        }
        write.table(x = tx2gene, file = tx2gene_path, sep = "\t", 
                    row.names = FALSE, quote = FALSE)
    }
    return(tx2gene_path)
}

#' @title Run tximport
#'
#' @description Run tximport. Will summarize abundance estimation from transcript 
#' level to gene level if `myAbundanceMetadata@txout == FALSE`. 
#' Otherwise keep abundance estimation at transcript level.
#'
#' @param myAbundanceMetadata A descendant object of the Class 
#' myAbundanceMetadata.
#' @param myBgeeMetadata A Reference Class BgeeMetadata object.
#' @param myUserMetadata A Reference Class UserMetadata object.
#' @param abundanceFile  (Optional) Path to the abundance file. NULL by default.
#' If not NULL, the file located at `abundanceFile` will be used to run tximport.
#' Otherwise (Default) the path to the abundance file is deduced fom attributes of
#' classes `BgeeMetadata`, `UserMetadata` and `AbundanceMetadata`
#' 
#' @author Julien Wollbrett
#'
#' @import rhdf5
#' @import tximport
#'
#' @export
#' 
#' @examples {
#' ah <- AnnotationHub()
#' ah_resources <- query(ah, c("Ensembl", "Caenorhabditis elegans", "84"))
#' annotation_object <- ah_resources[["AH50789"]]
#' user <- new("UserMetadata", species_id = "6239")
#' user <- setAnnotationFromObject(user, annotation_object, "WBcel235_84")
#' bgee <- new("BgeeMetadata", intergenic_release = "0.1")
#' kallisto <- new("KallistoMetadata")
#' abundance_file <- system.file("extdata", "abundance.tsv", package = "BgeeCall")
#' run_tximport(kallisto, bgee, user, abundance_file)
#' }
#'
run_tximport <- function (myAbundanceMetadata, myBgeeMetadata, myUserMetadata,
                          abundanceFile = NULL) {
    tx2gene_path <- create_tx2gene(myAbundanceMetadata, myBgeeMetadata, 
                                   myUserMetadata)
    tx2gene <- read.table(tx2gene_path,header = TRUE, sep = "\t")
    output_path <- get_tool_output_path(myAbundanceMetadata, myBgeeMetadata, 
                                        myUserMetadata)
    abundance_file <- file.path(output_path, myAbundanceMetadata@abundance_file)
    if (!file.exists(abundance_file)) {
        stop(paste0("can not generate presence/absence calls. 
                Abundance file is missing : ", abundance_file, "."))
    }
    txi <- tximport(abundance_file, type= myAbundanceMetadata@tool_name, 
                    tx2gene = tx2gene, txOut = myAbundanceMetadata@txOut, 
                    ignoreTxVersion = myAbundanceMetadata@ignoreTxVersion)
    return(txi)
}

abundance_without_intergenic <- function (myAbundanceMetadata, 
                                          myBgeeMetadata, myUserMetadata) {
    file_without_intergenic_name <- "abundance_without_intergenic.tsv"
    
    # remove intergenic from tx2gene
    tx2gene_path <- create_tx2gene(myAbundanceMetadata, 
                                   myBgeeMetadata, myUserMetadata)
    tx2gene <- read.table(tx2gene_path,header = TRUE, sep = "\t")
    tx2gene_without_intergenic <- 
        tx2gene[as.character(tx2gene$TXNAME)!=as.character(tx2gene$GENEID),]
    
    # remove intergenic from abundance file
    output_path <- get_tool_output_path(myAbundanceMetadata, 
                                        myBgeeMetadata, myUserMetadata)
    abundance_file <- file.path(output_path, "abundance.tsv")
    abundance <- read.table(abundance_file, header = TRUE, sep = "\t")
    abundance_without_intergenic <- 
        abundance[which(abundance[[myAbundanceMetadata@transcript_id_header]] 
                        %in% tx2gene_without_intergenic$TXNAME),]
    temp_abundance_file_without_intergenic <- 
        file.path(output_path, file_without_intergenic_name)
    write.table(abundance_without_intergenic, temp_abundance_file_without_intergenic, 
                sep = "\t", row.names = FALSE )
    
    # calculate corrected TPM value
    abundance_without_intergenic[ myAbundanceMetadata@abundance_header ] <- 
        countToTpm(abundance_without_intergenic[[myAbundanceMetadata@count_header]], 
                   abundance_without_intergenic[[myAbundanceMetadata@eff_length_header]])
    txi_without_intergenic <- 
        tximport(temp_abundance_file_without_intergenic, 
                 type= myAbundanceMetadata@tool_name, 
                 tx2gene = tx2gene_without_intergenic, 
                 txOut = myAbundanceMetadata@txOut, 
                 ignoreTxVersion = myAbundanceMetadata@ignoreTxVersion)
    file.remove(temp_abundance_file_without_intergenic)
    return(txi_without_intergenic)
}
