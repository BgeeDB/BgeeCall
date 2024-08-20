#' @title subsetting the annotation file
#'
#' @description Retrieves a specific part of the annotation field from the GTF file
#'
#' @param split_annotation the whole field we want to subset, for example:gene_id "FBgn0264003"; gene_name "mir-5613"; gene_source "FlyBase"; gene_biotype "pre_miRNA" 
#' @param field_name the field we want to filter the annotation with
#'
#' @author Alessandro Brandulas Cammarata
#' @author Julien Wollbrett
#' @author Julien Roux
#' @author Marta Rosikiewicz
#'
#' @return the subsetted data based on the provided field, for example: gene_id: FBgn0264003 
#' 
#' @noMd
#' @noRd
#' 
get_annot_value <- function(split_annotation, field_name){
  ## find the right field
  field_all <- split_annotation[grep(field_name, split_annotation, fixed = T)];

  ## split the field
  field_value <- strsplit(field_all, ' ', fixed = T)[[1]][2];

  ## remove the last ';' if necessairy
  field_value <- sub(';', '', field_value,fixed = T)

  return(field_value)
}

#' @title Removes N from intergenic regions 
#'
#' @description Removes intergenic regions which contain too many Ns instead of ACTG
#'
#' @param chr_number the id of the chromosome/contig
#' @param chr_intergenic_regions intergenic regions defined using gtf file
#' @param max_block_size threshold on maximum size of a block of N. If a sequence contains a block of N bigger or equals to this threshold we remove the block and split the
##	 sequence in 2. It will potentially create 2 intergenic sequences (if the block is not at the beginning nor the end of the sequence).
#' @param max_proportion_N maximum proportion of N in a sequence. If the proportion of N in a sequence is bigger than this threshold we do not keep the sequence.
#' @param min_intergenic_length minimum length of an intergenic region we want to keep. We remove all intergenic sequences smaller than this minimmum length.
#'
#' @author Alessandro Brandulas Cammarata
#' @author Julien Wollbrett
#' @author Julien Roux
#' @author Marta Rosikiewicz
#' 
#' @noMd
#' @noRd

remove_Ns_from_intergenic <- function (chr_number, chr_intergenic_regions, max_block_size = 31, max_proportion_N = 0.05, min_intergenic_length = 999) {
  ## Intergenic regions after removing unwanted N
  intergenic_regions_without_N <- matrix(ncol = 6, nrow = 0)
  colnames(intergenic_regions_without_N) <- c("chr", "start", "end", "upstream/downstream","strand", "sequence")

  # For each intergenic sequence
  for (line in 1 : nrow(chr_intergenic_regions)) {

    # absolut start and end positions in the chromosome
    chr_start <- as.numeric(chr_intergenic_regions[line,"start"])
    chr_end <- as.numeric(chr_intergenic_regions[line,"end"])

    # start position in the intergenic sequence. Allows to retrieve the sequence
    intergenic_start <- 1

    # retrieve intergenic sequence from fasta file using subseq function from BioStrings library
    intergenic_sequence <- chr_intergenic_regions[line,"sequence"]

    # if less N than the minimum size of a block we add the intergenic sequence to corrected intergenic regions
    if (max_block_size >=  count_number_of_occurences("N", intergenic_sequence)) {
      intergenic_regions_without_N <- rbind(intergenic_regions_without_N, chr_intergenic_regions[line,])

      # if no block of N and proportion of N lower than threshold then we add the sequence to corrected intergenic regions
    } else if (count_number_of_occurences(strrep("N", max_block_size), intergenic_sequence) ==  0
               & proportion_of_N(intergenic_sequence) <=  max_proportion_N) {
      intergenic_regions_without_N <- rbind(intergenic_regions_without_N, chr_intergenic_regions[line,])

    # Now need to parse the sequence in order to find potential block of N
    } else {

      # for each bp of the sequence
      seq_position <- 0
      while (seq_position + 1 <=  nchar(intergenic_sequence)) {
        seq_position <- seq_position + 1

        # start a new block of N if the current bp is a N
        if (substr(intergenic_sequence, seq_position, seq_position) ==  "N") {
          block_size = 1

          # continue to increase size of the block of N until we are at a real bp position or at the end of the sequence
          while (seq_position + 1 <=  nchar(intergenic_sequence) && substr(intergenic_sequence,seq_position + 1,seq_position + 1) ==  "N") {
            block_size <- block_size + 1
            seq_position <- seq_position + 1
          }

          # check the size of the block of N to know if it should be removed
          if (block_size >=  max_block_size) {

            #absolute position in the chromosome
            current_chr_stop <- as.numeric(chr_intergenic_regions[line,"start"]) + seq_position - (block_size + 1)
            # position in the intergenic sequence (the first one in this addition correspond to the start position of the sequence)
            current_intergenic_stop <- seq_position - block_size
            current_size <- (current_chr_stop - chr_start) + 1
            # check size and proportion of N of the subsequence
            if (current_size >=  min_intergenic_length && proportion_of_N(substr(intergenic_sequence, intergenic_start, current_intergenic_stop)) <=  max_proportion_N) {
              intergenic_regions_without_N <- rbind(intergenic_regions_without_N,
                                                    cbind(chr_number, chr_start, current_chr_stop, substr(intergenic_sequence, intergenic_start, current_intergenic_stop))[1,])
            }
            intergenic_start <- 1 + seq_position
            chr_start <- as.numeric(chr_intergenic_regions[line,"start"]) + seq_position
          }
        }
      }
      # test if it remains one intergenic region to add
      last_portion_size <- (seq_position - intergenic_start) + 1
      # check size and proportion of N of the subsequence
      if (last_portion_size >=  min_intergenic_length && proportion_of_N(substr(intergenic_sequence, intergenic_start, seq_position)) <=  max_proportion_N) {
        intergenic_regions_without_N <- rbind(intergenic_regions_without_N,
                                              cbind(chr_number, chr_start, chr_end, substr(intergenic_sequence, intergenic_start, seq_position))[1,])
      }
    }

  }
  return(intergenic_regions_without_N)
}

#' @title Calculate the percentage of N in a sequence 
#'
#' @description Calculate the percentage of N in a sequence
#'
#' @param sequence dna sequence to analyze
#'
#' @author Alessandro Brandulas Cammarata
#' @author Julien Wollbrett
#' @author Julien Roux
#' @author Marta Rosikiewicz
#' 
#' @noMd
#' @noRd
proportion_of_N <- function (sequence) {
  return (countPattern("N", sequence) / nchar(sequence))
}

#' @title Counts the number of occurence of a pattern
#'
#' @description counts the number of times one string is present in another string (e.g count_number_of_occurences("N", "ATGNNCTN") ==  3)
#'
#' @param pattern pattern to find in the string
#' @param string string in which the pattern could be present
#'
#' @author Alessandro Brandulas Cammarata
#' @author Julien Wollbrett
#' @author Julien Roux
#' @author Marta Rosikiewicz
#' 
#' @noMd
#' @noRd
count_number_of_occurences <- function (pattern, string) {
  return (lengths(regmatches(string, gregexpr(pattern, string))))
}

#' @title generates the initial intergenic regions
#'
#' @description generates intergenic regions from the species GTF and FASTA file by taking intergenic regions upstream and downstream of each gene
#'
#' @param gene_gtf_path full path to input gene gtf file
#' @param genome_fasta_path full path to input genome fasta file
#' @param N_block_size number of successive N from which it is considered as a block of N (and removed)
#' @param N_proportion higher proportion of N than this threshold results in removing the sequence
#' @param output_gtf_path full path to output folder + base name for output files
#' @param download_gtf_fasta boolean on whether to download the gtf and fasta from ensembl or
#' @param Rout_path path to  output .Rout file (optional)
#'
#' @author Alessandro Brandulas Cammarata
#' @author Julien Wollbrett
#' @author Julien Roux
#' @author Marta Rosikiewicz
#' 
#' @import GenomicFeatures
#' @import Biostrings
#'
#' 
#' @noMd
#' @noRd
#' 
generate_initial_intergenic_regions <- function(gene_gtf_path="./genomes/Homo_sapiens.GRCh38.gtf.gz", genome_fasta_path="./genomes/Homo_sapiens.GRCh38.genome.fa", N_block_size=31, N_proportion=0.05, output_gtf_path="./intergenic_regions/", Rout_path="./") {
print(sprintf("Generating intergenic regions for genome: %s", gene_gtf_path))
## reading in fasta file
cat("Reading FASTA file...\n")
ref_fasta_genome <- readDNAStringSet(genome_fasta_path)

# header of chromosomes/contigs has to be changed in order to map name of chr/contig present in gtf file.
# Only the name of chs/contig has to be kept
# regex allowing to retrieve name of the chs/contig (keep all characters before the first space)
chr_contig_regex <- "^([^ ]+).+"
# change header of each sequence
names(ref_fasta_genome) <- sub(chr_contig_regex, "\\1", names(ref_fasta_genome))

## reading in gtf file (gzipped, no need to uncompress)
cat("Reading GTF file...\n")
gene_gtf <- as.matrix(read.table(file = gzfile(gene_gtf_path, "r"), sep = "\t", strip.white = TRUE, as.is = TRUE, colClasses = "character", comment.char = '#'))
## Header lines starting with # are not read

## selecting exon lines
gene_gtf_exon <- gene_gtf[gene_gtf[,3] == "exon",]
## selecting gene lines
gene_gtf <- gene_gtf[gene_gtf[,3] == "gene",]

##  selecting only genes from assembled chromosome sequence (potentially useful code to select only fully asembled genome sequence)
##  chromosomes_all <- c(1:100, "X", "Y", "MT")
##  gene_gtf_exon <- gene_gtf_exon[gene_gtf_exon[,1] %in% chromosomes_all,]

## getting gene start, stop, chromosome, strand, biotype for each gene (from the exon GTF)
## GTF format reports 1-based coordinates (start and end)
cat("Extracting gene informations...\n")
## splitting the annotation field using single space "; " as a pattern
split_annotation_list <- strsplit(gene_gtf_exon[,9], "; ",  fixed  = T)

## getting the vector of the gene IDs (1 for every exon)
gene_ids <- sapply(split_annotation_list, function(x){ get_annot_value(x, 'gene_id') })

## getting the vector of the transcript IDs (1 for every exon)
transcript_ids <- sapply(split_annotation_list, function(x){ get_annot_value(x, 'transcript_id') })

## getting the table with mappings between gene IDs and transcript (for tximport)
tx2gene_ids <- unique(cbind(transcript_ids, gene_ids), MARGIN = 1)

## getting the vector of gene_biotypes: splitting gene gtf
split_annotation_list <- strsplit(gene_gtf[,9], "; ",  fixed  = T)
gene_biotypes <- cbind(
  sapply(split_annotation_list, function(x){ get_annot_value(x, 'gene_id') }),
  sapply(split_annotation_list, function(x){ get_annot_value(x, 'gene_biotype') }),
  "genic")

## getting the chromosome, start and the end of the gene
## For start, take the minimum exon start position
## For end, take the maximum exon end position
gene_start <- sapply(split(as.numeric(gene_gtf_exon[,4]), gene_ids), function(x){ sort(as.numeric(x))[1] })
gene_stop <- sapply(split(as.numeric(gene_gtf_exon[,5]), gene_ids), function(x){ rev(sort(as.numeric(x)))[1] })
gene_chr <- sapply(split(gene_gtf_exon[,1], gene_ids), function(x){ x[1] })

## chromosome/contig names from given gtf files
chromosomes <- unique(gene_gtf[,1])
## removing patch contigs before selecting intergenic regions
chromosomes <- chromosomes[grep('PATCH', chromosomes, invert = TRUE, ignore.case = TRUE)]

###################################################################
## Matrix summarizing impact of N removal
summary_N_removal <- matrix(ncol = 4, nrow = 2, data = 0)
colnames(summary_N_removal) <- c("intergenic_regions", "total_bp", "total_N", "proportion_N")
rownames(summary_N_removal) <- c("before", "after")

# To avoid scientific notation we change the value of option "scipen" to 999. At the end of the script we will change to its initial value
scipen_initial_value <- getOption("scipen")
options(scipen = 999)

## getting the vector of the gene IDs (1 for every gene)
gene_ids <- sapply(split_annotation_list, function(x){ get_annot_value(x, 'gene_id') }) 

intergenic_chr = c()
intergenic_starts = c()
intergenics_end = c()
intergenic_name = c()
intergenic_strand = c()

#Defining regions at 500bp upstream and downstream of genes
for(gene_nrow in 1:nrow(gene_gtf)){
  intergenic_chr =  append(intergenic_chr, as.character(gene_gtf[gene_nrow, 1]))
  intergenic_chr =  append(intergenic_chr, as.character(gene_gtf[gene_nrow, 1]))
  intergenic_strand = append(intergenic_strand, as.character(gene_gtf[gene_nrow, 7]))
  intergenic_strand = append(intergenic_strand, as.character(gene_gtf[gene_nrow, 7]))
  if(gene_gtf[gene_nrow, 7] == "+"){
    #Getting intergenic region information for intergenic upstream of forward strand gene
    intergenic_starts = append(intergenic_starts, max(0, as.numeric(gene_gtf[gene_nrow, 4]) - 1500))
    intergenics_end = append(intergenics_end, max(0, as.numeric(gene_gtf[gene_nrow, 4]) - 500))
    intergenic_name = append(intergenic_name, paste0("upstream_", gene_ids[gene_nrow]))
    
    #Getting information for downstream intergenic region
    intergenic_starts = append(intergenic_starts, max(0, as.numeric(gene_gtf[gene_nrow, 5]) + 500))
    intergenics_end = append(intergenics_end, max(0, as.numeric(gene_gtf[gene_nrow, 5]) + 1500))
    intergenic_name = append(intergenic_name, paste0("downstream_", gene_ids[gene_nrow]))
  } else {
    #Getting intergenic region information for intergenic upstream of forward strand gene
    intergenic_starts = append(intergenic_starts, max(0, as.numeric(gene_gtf[gene_nrow,4]) - 1500))
    intergenics_end = append(intergenics_end, max(0, as.numeric(gene_gtf[gene_nrow, 4]) - 500))
    intergenic_name = append(intergenic_name, paste0("downstream_", gene_ids[gene_nrow]))
    
    #Getting information for downstream intergenic region
    intergenic_starts = append(intergenic_starts, max(0, as.numeric(gene_gtf[gene_nrow, 5]) + 500))
    intergenics_end = append(intergenics_end, max(0, as.numeric(gene_gtf[gene_nrow, 5]) + 1500))
    intergenic_name = append(intergenic_name, paste0("upstream_", gene_ids[gene_nrow]))
  }
  
}

## Select the set of intergenic regions
cat("Selecting set of intergenic regions...\n")

## This object will include the coordinates of the selected intergenic regions
final_intergenic_regions <- matrix(c(intergenic_chr, intergenic_starts, intergenics_end, intergenic_name, intergenic_strand), ncol = 5)
colnames(final_intergenic_regions) <- c("chr", "start", "end", "upstream/downstream", "strand")

##Final dataset after filtering
Reference_intergenic <- matrix(ncol = 6, nrow = 0)
colnames(Reference_intergenic) <- c("chr", "start", "end", "upstream/downstream", "strand", "sequence")

#Subsetting intergenic regions to regions overlapping with gene proximity
for(chr in chromosomes){
  print(chr)
  if(( sum(gene_chr == as.character(chr)) == 0 )){ print(paste0("skipped", chr)); next }
  gene_IR <- IRanges(start = gene_start[gene_chr == as.character(chr)], end = gene_stop[gene_chr == as.character(chr)])
  intergenic_all_IR <- slice(coverage(gene_IR), lower = 0, upper = 0, rangesOnly = TRUE)
  Max_chromosomal_distance = max(end(intergenic_all_IR)[length(intergenic_all_IR)], end(gene_IR)[length(gene_IR)])
  intergenic_gene1k_IR <- IRanges(start = intergenic_starts[intergenic_chr == as.character(chr)], end = intergenics_end[intergenic_chr == as.character(chr)])
  intergenic_ref_IR <- intersect(intergenic_all_IR, intergenic_gene1k_IR)
  intergenic_ref_IR <- intergenic_ref_IR[intergenic_ref_IR@width >= 1001]
  #print(intergenic_ref_IR@width)
  chr_intergenic_regions = final_intergenic_regions[final_intergenic_regions[,"chr"] == chr & final_intergenic_regions[, "start"] %in% intergenic_ref_IR@start & as.numeric(final_intergenic_regions[, "end"]) <= as.numeric(Max_chromosomal_distance), ]
  
  if(is.null(nrow(chr_intergenic_regions))){
    print(paste0("skipped", chr)); next
  } else if(length(duplicated(chr_intergenic_regions[, "start"])) != 0){
    chr_intergenic_regions = chr_intergenic_regions[!duplicated(chr_intergenic_regions[, "start"]), ]
  }
  
  if(( nrow(chr_intergenic_regions) == 0 || is.null(nrow(chr_intergenic_regions)) )){ print(paste0("skipped ", chr));next }
  
  chr_sequence <- as.character(ref_fasta_genome[[chr]])
  sequence = apply(chr_intergenic_regions, 1, function(x) substr(chr_sequence, x["start"], x["end"]))
  chr_intergenic_regions = cbind(chr_intergenic_regions, sequence)
  
  # Keep information of number of N, bp and number of intergenic regions before removing blocks of N
  summary_N_removal["before","total_N"] <- summary_N_removal["before","total_N"] + sum(apply(chr_intergenic_regions, 1, function(x) count_number_of_occurences("N", x["sequence"])))
  summary_N_removal["before","intergenic_regions"] <- summary_N_removal["before","intergenic_regions"] + nrow(chr_intergenic_regions)
  summary_N_removal["before","total_bp"] <- summary_N_removal["after","total_bp"] + sum(as.numeric(chr_intergenic_regions[,"end"]) - as.numeric(chr_intergenic_regions[,"start"]) + 1)
  
  
  # Remove blocks of N and intergenic regions with big proportion of N
  chr_intergenic_regions_after_N_removal <- remove_Ns_from_intergenic(chr, chr_intergenic_regions, as.numeric(N_block_size), as.numeric(N_proportion))
  
  # Keep information of number of N, bp and number of intergenic regions before removing blocks of N
  summary_N_removal["after","total_N"] <- summary_N_removal["after","total_N"] + sum(apply(chr_intergenic_regions_after_N_removal, 1, function(x) count_number_of_occurences("N", x["sequence"])))
  summary_N_removal["after","intergenic_regions"] <- summary_N_removal["after","intergenic_regions"] + nrow(chr_intergenic_regions_after_N_removal)
  summary_N_removal["after","total_bp"] <- summary_N_removal["after","total_bp"] + sum(as.numeric(chr_intergenic_regions_after_N_removal[,"end"]) - as.numeric(chr_intergenic_regions_after_N_removal[,"start"]) + 1)
  Reference_intergenic <- rbind(Reference_intergenic, chr_intergenic_regions_after_N_removal[,c(1:6)])
}
options(scipen = scipen_initial_value)

## preparing intergenic gtf data
cat("\nPreparing intergenic GTF data...\n")
intergenic_regions_gtf <- matrix(ncol = 9, nrow = nrow(Reference_intergenic))
intergenic_regions_gtf[,1] <- Reference_intergenic[,1]
intergenic_regions_gtf[,2] <- "intergenic"
intergenic_regions_gtf[,3] <- "exon"
intergenic_regions_gtf[,c(4,5)] <- Reference_intergenic[,c(2,3)]
intergenic_regions_gtf[,c(6,8)] <- "."

intergenic_regions_gtf[,7] <- Reference_intergenic[, 5]

## intergenic_id - chr "_" start "_" stop
intergenic_id <- Reference_intergenic[, 4]
gene_id <- paste("gene_id ", paste(intergenic_id, ";", sep = ""), sep = "")
transcript_id <- paste("transcript_id ", paste(intergenic_id, ";", sep = ""), sep = "")
intergenic_regions_gtf[,9] <- apply(cbind(gene_id, transcript_id), 1, function(x){ paste(x, collapse = " ") })

## Caluclate proportion of N before and after N removal
summary_N_removal[,"proportion_N"] <- summary_N_removal[,"total_N"] / summary_N_removal[,"total_bp"]

####################################
output_file_path <- file.path(output_gtf_path, basename(substring(gene_gtf_path, 1, (nchar(gene_gtf_path) -7))))

if (!dir.exists(output_gtf_path)) {
  dir.create(output_gtf_path, recursive = TRUE)
}


## Output:
##Reference intergenic object
message("write file : ", paste(output_file_path, "_interegenic", sep = ""))

write.table(x = Reference_intergenic,
            file = paste(output_file_path, "_interegenic", sep = ""),
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

##Intergenic only in fasta format
#We need to convert the Reference_intergenic object into a fasta file
convert_to_fasta(Reference_intergenic, paste(output_file_path, "_fasta_interegenic.fa", sep = ""))

##Intergenic only
message("write file : ", paste(output_file_path, "gtf_interegenic", sep = ""))

write.table(x = intergenic_regions_gtf,
            file = paste(output_file_path, "_gtf_intergenic", sep = ""),
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

## GTF file with both genic exons and intergenic regions
message("write file : ", paste(output_file_path, ".gtf_all", sep = ""))
write.table(x = rbind(gene_gtf_exon, intergenic_regions_gtf),
            file = paste(output_file_path, ".gtf_all", sep = ""),
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

## GTF file with only genic exons
message("write file : ", paste(output_file_path, ".gtf_transcriptome", sep = ""))
write.table(x = gene_gtf_exon,
            file = paste(output_file_path, ".gtf_transcriptome", sep = ""),
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

## Table summaryzing the step of block of Ns removal
message("write file : ", paste(output_file_path, ".Nremoval", sep = ""))
write.table(x = summary_N_removal,
            file = paste(output_file_path, ".Nremoval", sep = ""),
            sep = "\t",
            row.names = TRUE,
            col.names = TRUE,
            quote = FALSE)

## Table of mapping between transcript_id and gene_id
intergenic_tx2gene_ids <- cbind(intergenic_id, intergenic_id)
all_tx2gene_ids <- rbind(tx2gene_ids, intergenic_tx2gene_ids)
message("write file : ", paste(output_file_path, ".tx2gene", sep = ""))
write.table(x = all_tx2gene_ids,
            file = paste(output_file_path, ".tx2gene", sep = ""),
            sep = "\t",
            row.names = FALSE,
            col.names = c("TXNAME", "GENEID"),
            quote = FALSE)

## Table of mapping between gene_id and both biotype and type (genic or intergenic)
intergenic_gene_biotypes <- cbind(
  intergenic_id,
  NA,
  "intergenic")
gene_biotypes <- rbind(gene_biotypes, intergenic_gene_biotypes)
message("write file : ", paste(output_file_path, ".gene2biotype", sep = ""))
write.table(x = gene_biotypes,
            file = paste(output_file_path, ".gene2biotype", sep = ""),
            sep = "\t",
            row.names = FALSE,
            col.names = c("id", "biotype", "type"),
            quote = FALSE)

}


#' @title downloads GTF and FASTA file from ensembl and generates intergenic regions from them
#'
#' @description allows to generate the initial close to gene intergenic regions needed to make the presence/absence calls by first download the GTF and FASTA files from ensembl
#'
#' @param species_gtf list of genomes to download following from ensembl or path to the file containing the list of genomes
#' @param ensembl_release ensembl release from which we want to GTF files
#' @param ensembl_metazoa_release ensembl metazoa release from which we want to GTF files
#' @param gtf_dir path to where we want to store the GTF files
#' @param from_file boolean of whether species_gtf is a file or a list of genomes
#' @param intergenic_dir path to where the initial intergenic regions will be stored
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
generate_intergenic_with_ensembl <- function(species_gtf = c("homo_sapiens/Homo_sapiens.GRCh38", "gallus_gallus/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b"), ensembl_release = 112, ensembl_metazoa_release = 59, gtf_dir = "./genomes/", from_file = FALSE, intergenic_dir = "./intergenic_regions") {

    retrieve_fasta_gtf(species_gtf = species_gtf, ensembl_release = ensembl_release, ensembl_metazoa_release = ensembl_metazoa_release, outDir = gtf_dir, from_file = from_file)

    if (!dir.exists(intergenic_dir)) {
        dir.create(intergenic_dir, recursive = TRUE)
    }

    for (gtf in species_gtf) {

        gtf <- sub(".*/", "", gtf)
        fasta = paste0(gtf, ".genome.fa")
        gtf <- paste0(gtf, ".gtf.gz")

        generate_initial_intergenic_regions(gene_gtf_path=paste0(gtf_dir, "/", gtf), genome_fasta_path=paste0(gtf_dir, "/", fasta))
    }
}


#' @title Changed to format of the Reference_intergenic object into a fasta file
#'
#' @description We take the Reference_intergenic object and shift the columns or fuse the column so that the format is the same as a fasta file
#'
#' @param Reference_intergenic object containing the intergenic regions
#'
#' @author Alessandro Brandulas Cammarata
#' 
#' @import dplyr
#' @import readr
#' 
#' @noMd
#' @noRd
#' 
convert_to_fasta <- function(Ref, output_file) {
  # Initialize empty data frames and list
  upstream <- data.frame()
  downstream <- data.frame()
  fasta_seq <- list()
  
  # Open the output file for writing
  out_f <- file(output_file, 'w')
  
  # Initialize variables for looping
  i <- 1
  odd <- FALSE
  
  # Loop through each line in the input data frame
  for (line in 1:nrow(Ref)) {
    new_header <- paste0('>', Ref[i, "upstream/downstream"], '    ', Ref[i, "chr"], Ref[i, "strand"], '    ', Ref[i, "start"], '-', Ref[i, "end"])
    
    writeLines(new_header, out_f)
    fasta_seq[[i]] <- Ref[i, "sequence"]
    writeLines(Ref[i, "sequence"], out_f)
    
    i <- i + 1
  }
  
  # Close the output file
  close(out_f)
}

#' @title Generates reference intergenic regions by removing outlier intergenics
#'
#' @description We sum the TPM values of intergenic regions in specified directory containing multiple libraries and find a threshold to remove outlier intergenic regions based on protein coding genes
#'
#' @param intergenic_regions_path path to the intergenic regions
#' @param tx2gene_path path to the tx2gene file of the species
#' @param gene2biotype_path path to the gene2biotype file of the species
#' @param TPM_path path to the abundance files of the libraries
#' @param output_path path to the output folder
#'
#' @author Alessandro Brandulas Cammarata
#' 
#' @import dplyr
#' @import MASS
#' @import ggplot2
#' 
#' @noMd
#' @noRd
#'
#' 

find_reference_intergenic_regions <- function(intergenic_regions_path, tx2gene_path, gene2biotype_path, TPM_path, output_path="./reference_intergenic_regions", species="my_species") {
  print("Generating reference intergenic regions...")

  # Create output folder if it does not exist
  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
  }

  numLibs = 0

  # List all files in the directory and its subdirectories which contain the abundance files
  libraries <- list.files(path = TPM_path, 
                    pattern = "abundance.tsv", 
                    full.names = TRUE, 
                    recursive = TRUE)

  # Read tx2gene file
  tx2gene <- read.table(file = file.path(tx2gene_path), header = TRUE, sep = "\t")
  colnames(tx2gene) <- c("target_id", "gene_id")
  # Read gene2biotype file
  gene2biotype <- read.table(file = file.path(gene2biotype_path), header = TRUE, sep = "\t")
  # Merge tx2gene and gene2biotype by gene_id
  tx2gene <- merge(tx2gene, gene2biotype, by.x = "gene_id", by.y = "id")
  # Create template of data.frame that will contain summed counts and ponderated mean of eff_length
  firstLibId <- libraries[1]
  summed <- read.table(file = firstLibId, header = TRUE, sep = "\t")
  summed <- merge(x = summed , y = tx2gene, by = "target_id")
  #sanity check that transcript are always in the same order in abundance files
  summed <- summed[order(summed$target_id),]
  ## detect number of transcripts and number of libraries to create matrix with proper size
  number_transcripts <- nrow(summed)
  numberLib <- length(libraries)
  # create data frame at proper size in order not to add new columns for each library
  effec_length_info <- as.data.frame(matrix(nrow=number_transcripts, ncol=numberLib))
  raw_counts_info <- as.data.frame(matrix(nrow=number_transcripts, ncol=numberLib))

  for(libraryId in libraries){
    ## For each library use the file at transcriptID level
    file <- libraryId
    if (file.exists(file)){
      cat("  Reading the ", libraryId, "\n")
      numLibs = numLibs + 1
      ## read transcript level data for each library
      kallisto_transcript_current <- read.table(file, h=T, sep="\t")
      # Remove all the transcripts that are not in the tx2gene file
      kallisto_transcript_current <- kallisto_transcript_current %>% filter(target_id %in% tx2gene$target_id)
      #sanity check that transcript are always in the same order in abundance files
      kallisto_transcript_current <- kallisto_transcript_current[order(kallisto_transcript_current$target_id),]
      raw_counts_info[,numLibs] <- kallisto_transcript_current$est_counts
      effec_length_info[,numLibs] <- kallisto_transcript_current$eff_length
    } else {
      stop("  No data found from", libraryId, "\n")
    }
  }
  
  if (numLibs == 1){
    ## provide output directly from gene level
    file2 <- paste0(output_path, "/abundance_gene_level+fpkm+intergenic.tsv")
    summed <- read.table(file2, h=T, sep="\t")
    # add mapping of gene and 
  } else {
    ## sum the raw counts for all libraries
    summed$est_counts <- apply(raw_counts_info, 1, sum)
    
    ## calculate effect_length for each transcript by using weighted.mean...
    calculate_effect_length <- function(counts, effec_length){
      myWeightedMean <- c()
      ## transpose the matrix (this means each column is a transcript ID and each row is a library)
      rawcounts <- t(counts)
      raw_effeclength <- t(effec_length)
      
      for (i in 1:ncol(raw_effeclength)) {
        if (sum(rawcounts[,i]) == 0 ){
          ## provide the same weight for a transcript in case the est_count is always zero
          rawcounts[,i] <- ifelse(rawcounts[,i] == 0, 1)
          weightedMean <-  weighted.mean(x=raw_effeclength[,i], w=rawcounts[,i])
          myWeightedMean[i] <- weightedMean
        } else {
          weightedMean <-  weighted.mean(x=raw_effeclength[,i], w=rawcounts[,i])
          myWeightedMean[i] <- weightedMean
        }
      }
      return(myWeightedMean)
    }
    ## each row is the weighted.mean of eff_length for each transcriptID
    summed$eff_length <- as.numeric(as.character(calculate_effect_length(raw_counts_info, effec_length_info)))
    
    # calculate effect_length for each transcript by using weighted.mean...
    ## re-calculate TPM and FPKM after collect the weighted.mean of eff_length for each transcriptID and after sum the est_count for each transcriptID
    estCount_to_tpm <- function(est_count, effec_length){
      rate <- log(est_count) - log(effec_length)
      denom <- log(sum(exp(rate)))
      exp(rate - denom + log(1e6))
    }
    summed$tpm <- as.numeric(as.character(estCount_to_tpm(summed$est_counts, summed$eff_length)))
    
    estCount_to_fpkm <- function(est_count, effec_length){
      N <- sum(est_count)
      exp( log(est_count) + log(1e9) - log(effec_length) - log(N) )
    }
    summed$fpkm <- as.numeric(as.character(estCount_to_fpkm(summed$est_counts, summed$eff_length)))
    
    ## select intergenic regions
    intergenic_regions <- summed[summed$type == "intergenic", c("target_id", "est_counts", "tpm", "fpkm", "type", "biotype")]
    names(intergenic_regions)[1] <- "gene_id"
    ## Sum est_counts, TPM and FPKM for genic regions
    sumGenic <- summed[summed$type == "genic", c("gene_id", "est_counts", "tpm", "fpkm", "type", "biotype")]
    sumGenic <- aggregate(sumGenic[,2:4], list(sumGenic$gene_id), sum)
    names(sumGenic)[1] <- "gene_id"
    sumGenic$type <-  rep("genic", times=length(sumGenic$gene_id))
    select_genic <- dplyr::filter(summed, type == "genic")
    select_biotype_genic = distinct(select_genic, gene_id, .keep_all= TRUE)
    sumGenic$biotype <- select_biotype_genic$biotype
    ## Final Table with genic and intergenic information
    summed <- rbind(sumGenic, intergenic_regions)
  }
  
  ## if no library wass found for this species
  if (numLibs == 0){
    cat("  No library found")
    next
  }
  
  ## Else:
  cat("  Plotting density of aggregated data\n")
  ## Density plot of summed data
  pdf(file = paste0(output_path, "/distribution_TPM_genic_intergenic_sum_", species, ".pdf"), width = 6, height = 5)
  ## par(mar=c(5,6,1,1)) ## bottom, left, top and right margins
  ## density of log2(TPM) of summed data
  dens <- density(log2(na.omit(summed$tpm) + 10^-6))
  ## Subgroups densities. Visualization trick: we add an invisible set of points at x=-30, to make densities comparable
  ## genic regions
  dens_genic <- density(c(rep(-30, times=sum(summed$type != "genic")), log2(summed$tpm[summed$type == "genic"] + 10^-6)))
  ## protein-coding genes only (had to take care of NAs strange behavior)
  dens_coding <- density(c(rep(-30, times=sum(!summed$biotype %in% "protein_coding")), log2(summed$tpm[summed$biotype %in% "protein_coding"] + 10^-6)))
  ## intergenic
  dens_intergenic <- density(c(rep(-30, times=sum(summed$type != "intergenic")), log2(summed$tpm[summed$type == "intergenic"] + 10^-6)))
  ## Plot whole distribution
  plot(dens, ylim=c(0, max(c(dens$y, dens_genic$y[dens_genic$x > -15], dens_coding$y[dens_coding$x > -15], dens_intergenic$y[dens_intergenic$x > -15]))*1.1), xlim=c(-23, 21), lwd=2, main=paste0("species summed TPM Density plot", " (", numLibs, " libraries)"), bty="n", axes=T, xlab="log2(TPM + 10^-6)")
  ## Add subgroups distributions (genic, intergenic, etc):
  ## genic
  lines(dens_genic, col="firebrick3", lwd=2)
  ## protein-coding genes
  lines(dens_coding, col="firebrick3", lwd=2, lty=2)
  ## intergenic
  lines(dens_intergenic, col="dodgerblue3", lwd=2)
  ## legend
  legend("topleft", c(paste0("all (", length(summed[,1]),")"), paste0("genic (", sum(summed$type == "genic"), ")"), paste0("coding (", sum(summed$biotype %in% "protein_coding"), ")"), paste0("intergenic (", sum(summed$type == "intergenic"), ")")), lwd=2, col=c("black", "firebrick3", "firebrick3", "dodgerblue3"), lty=c(1, 1, 2, 1), bty="n")
  dev.off()
  
  ## Redo plot for log2(FPKMs) (probably not proportional anymore)
  pdf(file = paste0(output_path, "/distribution_FPKM_genic_intergenic_sum_", species, ".pdf"), width = 6, height = 5)
  ## par(mar=c(5,6,1,1)) ## bottom, left, top and right margins
  ## FPKM density of summed data
  dens <- density(log2(na.omit(summed$fpkm) + 10^-6))
  ## Subgroups densities. Visualization trick: we add an invisible set of points at x=-30, to make densities comparable
  ## genic regions
  dens_genic <- density(c(rep(-30, times=sum(summed$type != "genic")), log2(summed$fpkm[summed$type == "genic"] + 10^-6)))
  ## protein-coding genes only (had to take care of NAs strange behavior)
  dens_coding <- density(c(rep(-30, times=sum(!summed$biotype %in% "protein_coding")), log2(summed$fpkm[summed$biotype %in% "protein_coding"] + 10^-6)))
  ## intergenic
  dens_intergenic <- density(c(rep(-30, times=sum(summed$type != "intergenic")), log2(summed$fpkm[summed$type == "intergenic"] + 10^-6)))
  ## Plot whole distribution
  plot(dens, ylim=c(0, max(c(dens$y, dens_genic$y[dens_genic$x > -15], dens_coding$y[dens_coding$x > -15], dens_intergenic$y[dens_intergenic$x > -15]))*1.1), xlim=c(-23, 21), lwd=2, main=paste0("species summed FPKM Density plot", " (", numLibs, " libraries)"), bty="n", axes=T, xlab="log2(FPKM + 10^-6)")
  ## Add subgroups distributions (genic, intergenic, etc):
  ## genic
  lines(dens_genic, col="firebrick3", lwd=2)
  ## protein-coding genes
  lines(dens_coding, col="firebrick3", lwd=2, lty=2)
  ## intergenic
  lines(dens_intergenic, col="dodgerblue3", lwd=2)
  ## legend
  legend("topleft", c(paste0("all (", length(summed[,1]),")"), paste0("genic (", sum(summed$type == "genic"), ")"), paste0("coding (", sum(summed$biotype %in% "protein_coding"), ")"), paste0("intergenic (", sum(summed$type == "intergenic"), ")")), lwd=2, col=c("black", "firebrick3", "firebrick3", "dodgerblue3"), lty=c(1, 1, 2, 1), bty="n")
  dev.off()
  
  ## Redo plot for log(read counts)
  pdf(file = paste0(output_path, "/distribution_counts_genic_intergenic_sum_", species, ".pdf"), width = 6, height = 5)
  ## par(mar=c(5,6,1,1)) ## bottom, left, top and right margins
  ## density of read counts on summed data
  dens <- density(log2(na.omit(summed$est_counts) + 1))
  ## Subgroups densities. Visualization trick: we add an invisible set of points at x=-10, to make densities comparable
  ## genic regions
  dens_genic <- density(c(rep(-10, times=sum(summed$type != "genic")), log2(summed$est_counts[summed$type == "genic"] + 1)))
  ## protein-coding genes only (had to take care of NAs strange behavior)
  dens_coding <- density(c(rep(-10, times=sum(!summed$biotype %in% "protein_coding")), log2(summed$est_counts[summed$biotype %in% "protein_coding"] + 1)))
  ## intergenic
  dens_intergenic <- density(c(rep(-10, times=sum(summed$type != "intergenic")), log2(summed$est_counts[summed$type == "intergenic"] + 1)))
  ## Plot whole distribution
  plot(dens, ylim=c(0, max(c(dens$y, dens_genic$y[dens_genic$x > -15], dens_coding$y[dens_coding$x > -15], dens_intergenic$y[dens_intergenic$x > -15]))*1.1), xlim=c(-1, 20), lwd=2, main=paste0("species summed read counts Density plot", " (", numLibs, " libraries)"), bty="n", axes=T, xlab="log2(read counts + 1)")
  ## Add subgroups distributions (genic, intergenic, etc):
  ## genic
  lines(dens_genic, col="firebrick3", lwd=2)
  ## protein-coding genes
  lines(dens_coding, col="firebrick3", lwd=2, lty=2)
  ## intergenic
  lines(dens_intergenic, col="dodgerblue3", lwd=2)
  ## legend
  legend("topleft", c(paste0("all (", length(summed[,1]),")"), paste0("genic (", sum(summed$type == "genic"), ")"), paste0("coding (", sum(summed$biotype %in% "protein_coding"), ")"), paste0("intergenic (", sum(summed$type == "intergenic"), ")")), lwd=2, col=c("black", "firebrick3", "firebrick3", "dodgerblue3"), lty=c(1, 1, 2, 1), bty="n")
  dev.off()
  
  
  ## Remove outlier intergenic regions that overlap too much with coding regions
  cat("  Definition of Reference intergenic regions\n")
  
  ## Focus on regions with enough signal (remove TPM = 0 or very small)
  summed_filtered <- summed[summed$tpm > 10^-6, ]
  
  ## open PDF device
  pdf(file = paste0(output_path, "/distribution_TPM_genic_intergenic_sum_deconvolution_", species, ".pdf"), width = 6, height = 5)
  
  set.seed(123)
  ## Coding regions
  # Filter for protein coding genes and log transform the TPM values
  pcoding_data <- log2(summed_filtered$tpm[summed_filtered$biotype %in% "protein_coding"])
  # Plot the TPM distribution of protein coding genes using a density plot
  ggplot(data.frame(log2_TPM = pcoding_data), aes(x = log2_TPM)) +
    geom_density() +
    labs(x = "log2(TPM) - protein-coding genes", y = "Density") +
    theme_minimal()
  
  ## Intergenic regions
  # Filter for intergenic regions and log transform the TPM values
  intergenic_data <- log2(summed_filtered$tpm[summed_filtered$type %in% "intergenic"])
  # Plot the TPM distribution of protein coding genes using a density plot
  ggplot(data.frame(log2_TPM = intergenic_data), aes(x = log2_TPM)) +
    geom_density() +
    labs(x = "log2(TPM) - protein-coding genes", y = "Density") +
    theme_minimal()
  
  ## Plot the density of the original data, and the density of regions classified to different gaussians
  cat("  Plotting density of deconvoluted genic and intergenic regions\n")
  dens <- density(log2(summed_filtered$tpm))
  ## Plot whole distribution
  plot(dens, ylim=c(0, max(dens$y)*1.1), xlim=c(-20, 20), lwd=2, main=paste0("species summed TPM Density plot", " (", numLibs, " libraries)"), bty="n", axes=T, xlab="log2(TPM)")
  
  ## protein-coding genes only (had to take care of NAs strange behavior)
  dens_coding <- density(log2(summed_filtered$tpm[summed_filtered$biotype %in% "protein_coding"]))
  ## Normalize density for number of observations
  dens_coding$y <- dens_coding$y * sum(summed_filtered$biotype %in% "protein_coding") / length(summed_filtered$tpm)
  lines(dens_coding, col="firebrick3", lwd=2, lty=2)
  
  
  ## intergenic
  dens_intergenic <- density(log2(summed_filtered$tpm[summed_filtered$type == "intergenic"]))
  dens_intergenic$y <- dens_intergenic$y * sum(summed_filtered$type == "intergenic") / length(summed_filtered$tpm)
  lines(dens_intergenic, col="dodgerblue3", lwd=2)

  ## legend
  legend("topleft", c(paste0("all (", length(summed_filtered[,1]),")"), paste0("coding (", sum(summed_filtered$biotype %in% "protein_coding"), ")"), paste0("intergenic (", sum(summed_filtered$type == "intergenic"), ")")), lwd=2, col=c("black", "firebrick3", "dodgerblue3"), lty=c(1, 2, 1), bty="n")
  
  
  ## Misclassification of intergenic regions based on protein coding overlap
  ## Removes intergenic regions which TPM is more likely to belong to a coding regions than an intergenic region (PDF coding >= PDF intergenic)
  
  #Points where the Kernel Density estimation will be evaluated (we use -15, 15 since by emperical evidence all RNAseq TPM data fall within that range)
  # Fit KDE for each distribution
  kde1 <- density(intergenic_data, n=2048, from=-15, to=15)
  kde2 <- density(pcoding_data, n=2048, from=-15, to=15)
  
  # Save the pdf of the kde of intergenic and protein coding into their own variables
  pdf1 <- kde1$y
  pdf2 <- kde2$y
  
  # Find the crossover point where PDF of protein-coding is always higher than PDF of intergenics
  crossover_index <- which(pdf2 > pdf1)[1]
  # Throw an error if no crossover point is found
  if (crossover_index == 0) {
    stop("No crossover point found where pdf2 is higher than pdf1.")
  }
  crossover_value <- kde1$x[crossover_index]
  
  # Identify outlier intergenic regions based on the crossover value
  mislabelled_values <- intergenic_data[intergenic_data >= crossover_value]
  
  print(paste("Crossover value where PDF of dist2 becomes higher than PDF of dist1:", crossover_value))
  print(paste("Number of mislabelled intergenic regions identified by KDE crossover point:", length(mislabelled_values)))
  #print(mislabelled_values)
  
  # Prepare data for visualization
  data <- data.frame(value = c(intergenic_data, pcoding_data),
                     distribution = factor(c(rep('intergenic PDF', length(intergenic_data)), rep('Protein Coding PDF', length(pcoding_data)))))
  
  mislabelled_df <- data.frame(value = mislabelled_values, y = rep(0, length(mislabelled_values)))
  
  ggplot() +
    geom_density(data = data, aes(x = value, fill = distribution), alpha = 0.6, adjust = 1.5) +
    geom_vline(aes(xintercept = crossover_value, color = 'Threshold of exclusion'), linetype = 'dashed', size = 1) +
    #geom_point(data = mislabelled_df, aes(x = value, y = y), color = 'red', size = 3, alpha = 0.8) +
    scale_color_manual(name = '', values = c('Threshold of exclusion' = 'black')) +
    labs(title = 'Intergenic regions exclusion based on protein coding overlap', x = 'Value', y = 'Density') +
    theme_minimal() +
    scale_fill_manual(values = c('intergenic PDF' = 'blue', 'Protein Coding PDF' = 'green')) +
    theme(
      legend.title = element_blank(),
      legend.position = c(0.2, 0.8),  
      legend.text = element_text(size = 7),  
      legend.background = element_rect(fill = alpha('white', 0.6)) 
    )
  
   dev.off()
  
  ###########
  
  ## Export file with summed data and classification of intergenic and coding regions
  cat("  Exporting aggregated data and classification of coding and intergenic regions\n")
  ## Add new column to summed object
  summed$classification <- NA
  summed$classification[summed$tpm > 10^-6 & summed$biotype %in% "protein_coding"] <- "coding"
  summed$classification[(summed$tpm < 10^-6 | summed$tpm >= 2**crossover_value) & summed$type == "intergenic"] <- "Outlier_intergenic"
  summed$classification[summed$tpm > 10^-6 & summed$tpm < 2**crossover_value & summed$type == "intergenic"] <- "Reference_intergenic"
  write.table(summed, file = paste0(output_path, "/sum_abundance_gene_level+fpkm+intergenic+classification_", species, ".tsv"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
  
  rm(summed)
  rm(summed_filtered)
  rm(numLibs)
}
#' @title Generates reference intergenic regions by removing outlier intergenics
#'
#' @description We transform the fasta file containing all intergenic regions to the reference intergenic regions by removing outliers
#'
#' @param intergenic_fasta intergenic regions fasta file
#' @param summed_tpm_file path to the file containing the summed reads for intergenic regions
#' @param output_file path to the output reference intergenic fasta file with removed outliers
#'
#' @author Alessandro Brandulas Cammarata
#' 
#' @import dplyr
#' @import stringr
#' 
#' @noMd
#' @noRd
#'
#' 
generate_fasta_intergenic_regions <- function(intergenic_fasta, summed_tpm_file, output_file) {
  ## Read the summed TPM file
  summed_tpm <- read.table(file = summed_tpm_file, header = TRUE, sep = "\t")

  # Filter for intergenic regions where classification matches the pattern 'Reference_intergenic'
  reference_intergenic <- summed_tpm %>% filter(str_detect(classification, 'Reference_intergenic')) %>% pull(gene_id)

  # Open the input FASTA file and output file
  fasta_output = file(output_file, "w")
  fasta_input = file(intergenic_fasta, "r")

  # Initialize the variable describing whether the current sequence is in the reference intergenic regions
  core <- FALSE

  # Loop through each line in the input file
  while (TRUE) {
    line <- readLines(fasta_input, n = 1)
    if (length(line) == 0) break # End of file

    # Check if the line is a header line 
    if (startsWith(line, ">")) {
      # Extract the transcript ID from the line
      enst_id <- str_extract(line, '[a-z]+_[A-Za-z]+\\d+')
      print(enst_id)
      
      # Check if the extracted transcript ID is in the reference intergenic regions
      if (enst_id %in% reference_intergenic) {
        # Write the header line to the output file
        writeLines(line, fasta_output)
        core <- TRUE
      } else {
        core <- FALSE
      }
    } else {
      # If it's not a header line, write it to the output file if core is TRUE
      if (core) {
        writeLines(line, fasta_output)
        core <- FALSE
      }
    }
  }
}


#' @title Generates reference intergenic regions by removing outlier intergenics
#'
#' @description We sum the TPM values of intergenic regions in specified directory containing multiple libraries and find a threshold to remove outlier intergenic regions based on protein coding genes
#' we then filter the fasta file containing all intergenic regions to the reference intergenic regions by removing outliers
#'
#' @param intergenic_regions_path path to the intergenic regions
#' @param tx2gene_path path to the tx2gene file of the species
#' @param gene2biotype_path path to the gene2biotype file of the species
#' @param TPM_path path to the abundance files of the libraries
#' @param output_path path to the output folder
#'
#' @author Alessandro Brandulas Cammarata
#' 
#' @import dplyr
#' @import MASS
#' @import ggplot2
#' 
#' @noMd
#' @noRd
#'
#' 

generate_reference_intergenic_regions <- function(intergenic_regions_path, tx2gene_path, gene2biotype_path, TPM_path, output_path="./reference_intergenic_regions", species="my_species") {

  find_reference_intergenic_regions(intergenic_regions_path, tx2gene_path, gene2biotype_path, TPM_path, output_path, species)

  # Extract base name from intergenic_regions_path and append "_fasta_intergenic.fa"
  base_name <- basename(intergenic_regions_path)
  fasta_file_name <- paste0(base_name, "_fasta_intergenic.fa")

  # Construct the full path for the fasta file
  fasta_intergenic <- file.path(output_path, fasta_file_name)

  generate_fasta_intergenic_regions(intergenic_fasta = intergenic_regions_path, summed_tpm_file = paste0(output_path, "/sum_abundance_gene_level+fpkm+intergenic+classification_", species, ".tsv"), paste0(output_path, "/", substr(base_name, 1, nchar(base_name) - 3), "_reference.fa"))
}