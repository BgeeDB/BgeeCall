#' @title AbundanceMetadata s4 class
#'
#' @description An S4 class that is the parent class of all abundance tool 
#' Classes. 
#' It contains information needed to all abundance tools.
#' This class can be seen as an abstract class, you should never instanciate 
#' it.
#'
#' @slot txOut Similar to tximport txOut parameter. Allows to keep abundance at
#' transcript level if TRUE (default = FALSE)
#' @slot ignoreTxVersion logical used to remove transcript version in 
#' transcript ID if TRUE (default = FALSE)
#' @slot cutoff_type Defines the approach used to generate present/absent calls.
#' default value is 'pValue', allowing calls to be generated using a pValue.
#' Other possible values are 'intergenic' allowing to use a ratio of intergenic
#' sequences considered as present as a threshold, or use qValue allowing calls
#' to be generated from a qValue.
#' @slot cutoff numeric value of the cutoff used to generate the present/absent
#' calls. If value of the slot cutoff_type is 'pValue' this cutoff will correspond
#' to the highest pValue allowing to define a gene as present. If value of the slot 
#' cutoff_type is 'intergenic' this cutoff will correspond to the proportion of 
#' intergenic present divided by proportion of protein coding present. 
#' If value of the slot cutoff_type is 'qValue' this cutoff will correspond to the
#' highest qValue allowing to define a gene as present. The qValue is calculated
#' based on the proportion of intergenic/(intergenic + genic) at each unique abundance value (TPM).
#' The default value is 0.05. Be careful when modifying this value as it could have a huge impact 
#' on present/absent calls.
#' @slot full_transcriptome_file Name of the fasta file containing both 
#' transcriptomic and intergenic regions. This file is created by the pipeline.
#' You should edit this slot only if you already have such a file with a 
#' different name.
#' @slot tx2gene_file Name of the file containing the mapping between 
#' transcript IDs and gene IDs (See the tximport package vignette for more 
#' details). This file is created by the pipeline. You should edit this 
#' slot only if you already have such a file with a different name. This 
#' file must be store at get_species_path()
#' @slot tx2gene_file_without_version Name of the file containing the mapping 
#' between transcript IDs and gene IDs if ignoreTxVersion == TRUE (See the 
#' tximport package vignette for more details). This file is created by the 
#' pipeline. You should edit this slot only if you already have such a file 
#' with a different name. This file must be store at get_species_path()
#' @slot gene2biotype_file Name of the file containing the mapping between 
#' gene IDs and biotypes. This file is created by the pipeline. You should 
#' edit this slot only if you already have such a file with a different name.
#' @slot tool_name Name of the tool that will be use to generate transcript 
#' abundance estimation. All descendant of this class have to define a value 
#' for this slot (in the prototype section)
#' @slot abundance_file Name of the transcript-level abundance file. All 
#' descendant of this class have to define a value for this slot (in the 
#' prototype section)
#' @slot read_size_kmer_threshold read size of the library below which 
#' transcript index is created using a smaller kmer size
#' @slot transcript_id_header Name of the header of the column that contains 
#' transcript ID
#' @slot count_header Name of the header of the column that contains count
#' @slot abundance_header  Name of the header of the column that contains 
#' abundance
#' @slot eff_length_header Name of the header of the column that contains 
#' effective length
#' @slot transcript_calls_file_name default name of file containing all 
#' transcript ids and calls (if calls created at transcript level)
#' @slot gene_calls_file_name default name of file containing all gene ids and 
#' calls (if calls created at gene level)
#' @slot transcript_cutoff_file_name default name of file containing summary 
#' of cutoff used to generate transcript expression calls (if calls created 
#' at transcript level)
#' @slot gene_cutoff_file_name default name of file containing summary of 
#' cutoff used to generate gene expression calls (if calls created at gene 
#' level)
#' @slot transcript_distribution_file_name default name of density plot file 
#' containing TPM distribution of all transcripts  (if calls created at 
#' transcript level)
#' @slot gene_distribution_file_name default name of density plot file 
#' containing TPM distribution of all genes (if calls created at gene level)
#' 
#' 
AbundanceMetadata <- setClass(
    # Set the name for the class
    Class = "AbundanceMetadata",
    
    # Define the slots
    representation = representation(
        txOut = "logical",
        ignoreTxVersion = "logical",
        cutoff_type = "character",
        cutoff = "numeric",
        full_transcriptome_file = "character",
        tx2gene_file = "character",
        tx2gene_file_without_version = "character",
        tx2biotype_file = "character",
        tx2biotype_file_without_tx_version = "character",
        gene2biotype_file = "character",
        tool_name = "character",
        abundance_file = "character",
        abundance_file_without_tx_version = "character",
        read_size_kmer_threshold = "numeric",
        transcript_id_header = "character",
        count_header = "character",
        abundance_header = "character",
        eff_length_header = "character",
        transcript_calls_file_name = "character",
        gene_calls_file_name = "character",
        transcript_cutoff_file_name = "character",
        gene_cutoff_file_name = "character",
        transcript_distribution_file_name = "character",
        gene_distribution_file_name = "character"
    ),
    
    # Set the default values for the slots.
    prototype = prototype(
        txOut = FALSE,
        ignoreTxVersion = FALSE,
        cutoff_type = 'pValue',
        cutoff = 0.05,
        full_transcriptome_file = "transcriptome_with_intergenic.fa",
        tx2gene_file = "tx2gene.tsv",
        tx2gene_file_without_version = "tx2gene_without_tx_version.tsv",
        tx2biotype_file = "tx2biotype.tsv",
        tx2biotype_file_without_tx_version = "tx2biotype_without_tx_version.tsv",
        gene2biotype_file = "gene2biotype.tsv",
        transcript_calls_file_name = "transcript_level_abundance+calls",
        gene_calls_file_name = "gene_level_abundance+calls",
        transcript_cutoff_file_name = "transcript_cutoff_info_file",
        gene_cutoff_file_name = "gene_cutoff_info_file",
        transcript_distribution_file_name = 
            "transcript_TPM_genic_intergenic+cutoff",
        gene_distribution_file_name = "gene_TPM_genic_intergenic+cutoff"
    )
    
)
