#' @title Retrieve dataframe of reference intergenic IDs
#'
#' @description Retrieve a dataframe of reference intergenic IDs. 
#' It is generated using the reference intergenic fasta file from Bgee FTP.
#'
#' @param myBgeeMetadata A Class BgeeMetadata object.
#' @param myUserMetadata A Class UserMetadata object.
#'
#' @author Julien Wollbrett
#' @author Julien Roux
#'
#' @return A dataframe containing reference intergenic ids
#' 
#' @import Biostrings
#' 
#' @noMd
#' @noRd
#'
get_ref_intergenic_ids <- function(myBgeeMetadata, 
    myUserMetadata) {
    bgee_intergenic_file <- retrieve_intergenic_path(myBgeeMetadata, myUserMetadata)
    bgee_intergenic <- readDNAStringSet(bgee_intergenic_file)
    # keep only intergenic ids from fasta file
    return(as.data.frame(sub("^([^ ]+).*", "\\1", names(bgee_intergenic))))
}

#' @title Generate presence absence
#'
#' @description Generate presence absence calls. It correponds to 
#' the last part of the generation of the expression calls workflow. 
#' It runs the last part of the workflow generating present/absent 
#' expression calls. 
#' This function should only be used by advanced user who already 
#' manually run all previous parts of the pipeline.
#' If you are not an advanced user it is safer to run the function 
#' ``generate_calls_workflow`` that run all steps of the worklow
#'
#'@seealso generate_calls_workflow
#'
#' @param myAbundanceMetadata A descendant object of the Class myAbundanceMetadata
#' (optional).
#' @param myBgeeMetadata A Class BgeeMetadata object (optional).
#' @param myUserMetadata A Class UserMetadata object.
#'
#' @author Julien Wollbrett
#' @author Julien Roux
#' @author Sara Fonseca Costa
#' 
#' @return path to the 4 output files
#'
#' @export
#' 
#' @examples {
#' # this example reuse data present in the directory 'extdata' of the package.
#' user <- new('UserMetadata', working_path = system.file('extdata', 
#' package = 'BgeeCall'), species_id = '6239', rnaseq_lib_path = system.file(
#' 'extdata', 'SRX099901_subset', package = 'BgeeCall'), 
#' annotation_name = 'WBcel235_84', simple_arborescence = TRUE)
#' calls_output <- generate_presence_absence(myUserMetadata = user)
#' 
#' #
#' }
#'
generate_presence_absence <- function(myAbundanceMetadata = new("KallistoMetadata"), 
    myBgeeMetadata = new("BgeeMetadata"), myUserMetadata) {
    
    system.file()
    # load data
    ref_intergenic <- get_ref_intergenic_ids(myBgeeMetadata, 
        myUserMetadata)
    tool_path <- get_tool_path(myAbundanceMetadata, 
        myBgeeMetadata, myUserMetadata)
    
    # use the standard output dir or the one defined by the user
    output_path <- get_tool_output_path(myAbundanceMetadata, 
                                        myBgeeMetadata, myUserMetadata)
    
    # check if calls have to be created
    if (file.exists(file.path(output_path, 
                              paste0(myAbundanceMetadata@gene_calls_file_name, myAbundanceMetadata@cutoff_type, ".tsv"))) & !myAbundanceMetadata@overwrite_calls ){
        if(isTRUE(myUserMetadata@verbose)) {
            message("no need to regenerate calls")
        }
    }
    else {
        # biotype mapping information will depend on
        # summarization at gene level or not
        biotype_mapping <- ""
        if (myAbundanceMetadata@txOut) {
            biotype_mapping <- load_transcript_to_biotype(myAbundanceMetadata, 
                myBgeeMetadata, myUserMetadata)
        } else {
            biotype_mapping <- load_gene_to_biotype(myAbundanceMetadata, 
                myBgeeMetadata, myUserMetadata)
        }
    
        # run tximport for file with intergenic regions (if
        # myAbundanceMetadata@txOut = FALSE, then tximport
        # will summarize transcript level estimates at
        # gene level)
        tximportObject <- run_tximport(myAbundanceMetadata = myAbundanceMetadata, 
            myBgeeMetadata = myBgeeMetadata, myUserMetadata = myUserMetadata)
    
        # recalculate TPM without intergenic regions and
        # run tximport (if myAbundanceMetadata@txOut =
        # FALSE, then tximport will summarize at
        # transcript level estimates at gene level)
    
        # remove intergenic
        tximportObject_without_intergenic <- abundance_without_intergenic(myAbundanceMetadata, 
            myBgeeMetadata, myUserMetadata)
    
        # transform tximportObect in order to easily process information
        abundance <- transform_tximport(tximportObject, biotype_mapping)
        # transform tximportObject without intergenic in
        # order to easily process information
        abundance_without_intergenic <- transform_tximport(tximportObject_without_intergenic, 
            biotype_mapping)
    
        # define coding and intergenic abundance subset
        selected_coding <- abundance$biotype %in% "protein_coding"
        selected_intergenic <- (abundance$type == "intergenic" & 
            abundance$id %in% biotype_mapping$id)
        
        # actual call generation depending on cutoff_type
        if(isTRUE(myUserMetadata@verbose)) {
            message("Generate present/absent expression calls using ",
                myAbundanceMetadata@cutoff_type, "cutoff")
        }
        # init variables not used in all approaches. Will be used to add a line in the cutoff info file
        # if not null
        # only used in intergenic approach
        r_cutoff <- NULL
        #only used in pvalue approach
        mean_pvalue<- NULL
        #only used in pvalue approach
        sd_pvalue <- NULL


        if(myAbundanceMetadata@cutoff_type == 'intergenic') {
            if(isTRUE(myUserMetadata@verbose)) {
                results <- calculate_abundance_cutoff(abundance, 
                    selected_coding, selected_intergenic, myAbundanceMetadata@cutoff)
            } else {
                results <- suppressMessages(calculate_abundance_cutoff(abundance, 
                    selected_coding, selected_intergenic, myAbundanceMetadata@cutoff))
            }
            abundance_cutoff <- results[1]
            r_cutoff <- results[2]
            # add presence/absence information
            abundance$call <- ifelse(abundance$abundance >= 
                                         abundance_cutoff, "present", "absent")
            # transfert presence/absence annotations to the
            # abundances without intergenic
            abundance_without_intergenic <- merge(abundance_without_intergenic, 
                                                  abundance[, c("id", "call")], 
                                                  by = "id")
        # generate calls and calculate abundance_cutoff
        }else if (myAbundanceMetadata@cutoff_type == 'pValue') {
            pvalue_generated <- generate_theoretical_pValue(counts =abundance,
                                                     myAbundanceMetadata@cutoff)
            abundance <- pvalue_generated$counts_with_pValue
            mean_pvalue <- pvalue_generated$mean
            sd_pvalue <- pvalue_generated$sd
            abundance_cutoff <- min(na.omit(abundance$abundance[abundance$pValue <= myAbundanceMetadata@cutoff]))
            # abundances without intergenic
            abundance_without_intergenic <- merge(abundance_without_intergenic, 
                                                  abundance[, c("id", "zScore", "pValue", "call")], 
                                                  by = "id")
        } else if (myAbundanceMetadata@cutoff_type == 'qValue') {
            abundance <- generate_qValue(counts =abundance,
                                         myAbundanceMetadata@cutoff)
            
            abundance_cutoff <- min(na.omit(abundance$abundance[abundance$qValue <= myAbundanceMetadata@cutoff]))

            # abundances without intergenic
            abundance_without_intergenic <- merge(abundance_without_intergenic, 
                                                  abundance[, c("id", "qValue", "call")], 
                                                  by = "id")
        } else {
            stop("unknown cutoff type : ", myAbundanceMetadata@cutoffType, ". Should be 
            \"pValue\" or \"intergenic\" or \"qValue\"")
        }
   
        # generate name of output file
        calls_file_name <- paste0(myAbundanceMetadata@gene_calls_file_name, myAbundanceMetadata@cutoff_type, ".tsv")
        cutoff_info_file_name <- paste0(myAbundanceMetadata@gene_cutoff_file_name,  myAbundanceMetadata@cutoff_type, ".tsv")
        distribution_file_name <- paste0(myAbundanceMetadata@gene_distribution_file_name,  myAbundanceMetadata@cutoff_type, ".pdf")
        if (myAbundanceMetadata@txOut) {
            calls_file_name <- paste0(myAbundanceMetadata@transcript_calls_file_name,  myAbundanceMetadata@cutoff_type, ".tsv")
            cutoff_info_file_name <- paste0(myAbundanceMetadata@transcript_cutoff_file_name,  myAbundanceMetadata@cutoff_type, ".tsv")
            distribution_file_name <- paste0(myAbundanceMetadata@transcript_distribution_file_name,  myAbundanceMetadata@cutoff_type, ".pdf")
        }
    
        # generate pdf plot
        pdf(file = file.path(output_path, distribution_file_name), 
            width = 6, height = 5)
        plot_distributions(abundance, selected_coding, 
            selected_intergenic, abundance_cutoff, myUserMetadata)
    
        # generate cutoff info file
        cutoff_info_file <- cutoff_info(counts = abundance, column = "call", abundance_cutoff = abundance_cutoff, 
                                        r_cutoff = r_cutoff, mean_pvalue=mean_pvalue, sd_pvalue=sd_pvalue, 
                                        myUserMetadata = myUserMetadata, myAbundanceMetadata = myAbundanceMetadata)
        dev.off()
        

    
        # Save calls and stats to output folder
        calls_file_path <- file.path(output_path, calls_file_name)
        cutoff_info_file_path <- file.path(output_path, 
            cutoff_info_file_name)
        write.table(abundance_without_intergenic, file = calls_file_path, 
            quote = FALSE, sep = "\t", col.names = TRUE, 
            row.names = FALSE)
        write.table(t(t(cutoff_info_file)), file = cutoff_info_file_path, 
            quote = FALSE, sep = "\t", col.names = FALSE, 
            row.names = TRUE)
        #save values of the Slots of the BgeeCall internal S4 objects
        s4_summary_df <- as.data.frame(
            generate_S4_object_properties_output(myAbundanceMetadata,
                myBgeeMetadata,
                myUserMetadata))
        s4_slots_path <- file.path(output_path, "S4_slots_summary.tsv")
        write.table(s4_summary_df, file = s4_slots_path, quote = FALSE, sep = "\t",
                    col.names = TRUE, row.names = FALSE)    
        calls_result <- list()
        calls_result$calls_tsv_path <- calls_file_path
        calls_result$cutoff_info_file_path <- cutoff_info_file_path
        calls_result$abundance_tsv <- file.path(output_path, 
            myAbundanceMetadata@abundance_file)
        calls_result$TPM_distribution_path <- file.path(output_path, 
            distribution_file_name)
        calls_result$S4_slots_summary <- s4_slots_path
        return(calls_result)
        ## t(t(cutoff_info_file)) is a solution to export a
        ## vector vertically
    }
}

#' @title Transform tximport object
#'
#' @description transform tximport object in order to easily process information
#'
#' @param tximportObject A tximport Object
#' @param biotype_mapping Mapping between gene or transcript IDs and biotypes
#'
#' @return A tximport Object with biotype and without the countsFromAbundance 
#' column
#'
#' @noMd
#' @noRd
#'
transform_tximport <- function(tximportObject, biotype_mapping) {
    tx_df <- as.data.frame(tximportObject)
    tx_df$id <- rownames(tx_df)
    tx_df$countsFromAbundance <- NULL
    abundance <- merge(tx_df, biotype_mapping, by = "id", 
        all = FALSE)
}

