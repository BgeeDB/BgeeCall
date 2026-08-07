#' @importFrom stats aggregate approxfun complete.cases density integrate
#' @importFrom stats median na.omit p.adjust pnorm qnorm quantile sd
#' @importFrom stats weighted.mean
#' @importFrom dplyr %>% arrange bind_rows distinct mutate n pull rename
#' @importFrom dplyr summarise ungroup
#' @importFrom data.table data.table fread fwrite rbindlist setDT setkey
#' @importFrom data.table setnames
#' @importFrom Biostrings DNAString DNAStringSet readDNAStringSet
#' @importFrom Biostrings reverseComplement subseq width writeXStringSet
#' @importFrom IRanges IRanges coverage end slice
#' @importFrom GenomicFeatures exonsBy transcriptsBy
#' @importFrom txdbmaker makeTxDbFromGRanges
#' @importFrom rtracklayer import export
#' @importFrom biomaRt useEnsembl getBM
#' @importFrom readr cols parse_date read_tsv write_tsv
#' @importFrom sjmisc str_contains
#' @importFrom tximport tximport summarizeToGene
#' @importFrom methods new slot "slot<-" is
#' @importFrom ggplot2 aes element_blank element_rect element_text
#' @importFrom ggplot2 geom_density geom_vline ggplot labs scale_color_manual
#' @importFrom ggplot2 scale_fill_manual theme theme_minimal
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics arrows axis legend lines mtext par plot
#' @importFrom utils download.file head packageVersion read.table untar
#' @importFrom utils unzip write.table
#' @importFrom rslurm slurm_apply get_job_status get_slurm_out
#' @importFrom RCurl getURL
#' @importFrom stringr str_extract str_detect str_match
#' @importFrom tools file_path_sans_ext
#' @importFrom curl curl_download
#' @importFrom jsonlite fromJSON

NULL
