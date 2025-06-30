#' @importFrom stats density p.adjust pnorm qnorm sd median
#' @importFrom dplyr %>% arrange bind_rows distinct mutate n pull rename summarise ungroup
#' @importFrom data.table data.table fread fwrite rbindlist setDT setkey setnames
#' @importFrom Biostrings DNAString DNAStringSet readDNAStringSet writeXStringSet width subseq reverseComplement
#' @importFrom IRanges IRanges coverage slice
#' @importFrom GenomicFeatures makeTxDbFromGFF exonsBy transcriptsBy
#' @importFrom rtracklayer import export
#' @importFrom biomaRt useEnsembl getBM
#' @importFrom readr parse_date read_tsv write_tsv
#' @importFrom sjmisc str_contains
#' @importFrom tximport tximport summarizeToGene
#' @importFrom methods new slot "slot<-" is
#' @importFrom ggplot2 ggplot aes geom_density geom_vline scale_color_manual labs theme_minimal scale_fill_manual theme element_blank element_text element_rect
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics legend lines par plot
#' @importFrom utils read.table write.table
#' @importFrom rslurm slurm_apply get_job_status get_slurm_out
#' @importFrom RCurl getURL
#' @importFrom stringr str_extract str_detect str_match
#' @importFrom tools file_path_sans_ext
#' @importFrom curl curl_download
#' @importFrom jsonlite fromJSON

NULL 