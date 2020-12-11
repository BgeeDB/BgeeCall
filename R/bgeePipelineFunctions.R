# these functions should be part of
# presenceAbsence.R they mainly correspond to
# duplicated code between the pipeline and this
# package.  They do not need to be exported We
# decided to move them in a different file to
# easily update them when the pipeline is updated.
# XXX We should maybe centralise these functions
# somewhere outside of the package and the bgee
# pipeline

#' @title Calculate TPM cutoff
#'
#' @description This function calculate the TPM cutoff.
#' This cutoff will correspond to the minimal value of TPM for which the ratio
#' of genes and intergenic regions is equal to `intergenic_cutoff`
#' (by default = 0.05)
#' or lower (first test if at least 1 TPM value has this property):
#'
#' @param counts TPM information for both genic and intergenic regions
#' @param selected_coding TPM information for both protein_coding regions
#' @param selected_intergenic TPM information for intergenic regions
#' @param intergenic_cutoff value of the ratio between prop of
#' intergenic present and protein coding present (called r in this function)
#'
#' @return the TPM cutoff and the value of r (proportion of intergenic regions
#' considered as present )
#'
#' @author Julien Roux
#' @author Julien Wollbrett
#'
#' @noMd
#' @noRd
#'
calculate_abundance_cutoff <- function(counts,
    selected_coding,
    selected_intergenic,
    intergenic_cutoff = 0.05) {
    ## r = (number of intergenic regions with TPM values
    ## higher than x * number of coding regions) /
    ## (number of coding regions with TPM values higher
    ## than x * number of intergenic regions) = 0.05
    ## What is value of x (cutoff)? calculate the
    ## distribution of r for a range of TPMs, then
    ## select the closest value to `cutoff`
    ## Counting how many intergenic regions have equal
    ## or higher value of TPM for every value of TPM For
    ## each gene's TPM (sorted), calculate r
    summed_intergenic <-
        vapply(unique(sort(counts$abundance[selected_coding])),
        function(x) {
            return(sum(counts$abundance[selected_intergenic] >=
            x))
        }, numeric(1))
    ## It is not necessary to do the same for coding
    ## regions: for a sorted vector the number of
    ## position + 1.  Here, it is a bit trickier since
    ## we do not consider all coding TPM values, but the
    ## unique ones, so we use the rle function to know
    ## the lengths of the runs of similar values and sum
    ## them
    summed_coding <-
        c(0, cumsum(rle(sort(
            counts$abundance[selected_coding]
        ))$lengths))
    summed_coding <- summed_coding[-(length(summed_coding))]
    summed_coding <- sum(selected_coding) - summed_coding
    
    ## Now we can calculate r
    r <-
        (summed_intergenic / sum(selected_intergenic)) / (summed_coding / sum(selected_coding))
    
    prot_coding_present <- 100 - (intergenic_cutoff * 100)
    
    ## Select the minimal value of TPM for which the
    ## ratio of genes and intergenic regions is equal to
    ## `intergenic_cutoff` or lower (first test if at
    ## least 1 TPM value has this property):
    if (sum(r < intergenic_cutoff) == 0) {
        TPM_cutoff <- sort(unique(counts$abundance[selected_coding]))[which(r ==  min(r))[1]]
        r_cutoff <- min(r)
        warning(
            "There is no TPM cutoff for which ",
            prot_coding_present,
            "% of the expressed genes would be coding. TPM cutoff is fixed at
the first value with maximum coding/intergenic ratio. r=",
            r_cutoff,
            "at TPM=",
            TPM_cutoff,
            "\n"
        )
    } else {
        TPM_cutoff <- sort(unique(counts$abundance[selected_coding]))[which(r < intergenic_cutoff)[1]]
        r_cutoff <- intergenic_cutoff
        message("TPM cutoff for which ",prot_coding_present,
            "% of the expressed genes are coding found at TPM = ", TPM_cutoff)
    }
    return(c(TPM_cutoff, r_cutoff))
}

#' @title plot distribution of TPMs
#'
#' @description Plotting of the distribution of TPMs for coding and intergenic
#' regions + cutoff
#'
#' @param counts TPM information for both genic and intergenic regions
#' @param selected_coding TPM information for protein_coding regions
#' @param selected_intergenic TPM information for intergenic regions
#' @param cutoff TPM cutoff below which calls are considered as absent
#' @param myUserMetadata A Reference Class UserMetadata object.
#' This object has to be edited before running kallisto @seealso UserMetadata.R
#'
#' @author Julien Roux
#' @author Julien Wollbrett
#'
#' @noMd
#' @noRd
#'
#' @return plot distribution of TPMs for coding and intergenic regions + cutoff
#'
plot_distributions <- function(counts,
    selected_coding, selected_intergenic,
    cutoff, myUserMetadata) {
    ## Plotting of the distribution of TPMs for coding
    ## and intergenic regions + cutoff Note: this code
    ## is largely similar to plotting section in
    ## rna_seq_analysis.R
    
    par(mar = c(5, 6, 1, 1))  ## bottom, left, top and right margins
    dens <- density(log2(na.omit(counts$abundance) + 10 ^ -6))
    
    ## protein-coding genes only (had to take care of
    ## NAs strange behavior)
    dens_coding <- density(log2(counts$abundance[selected_coding] +
        10 ^ -6))
    ## Normalize density for number of observations
    dens_coding$y <-
        dens_coding$y * sum(selected_coding) / length(counts$abundance)
    
    ## intergenic
    dens_intergenic <-
        density(log2(counts$abundance[selected_intergenic] + 10 ^ -6))
    dens_intergenic$y <-
        dens_intergenic$y * sum(selected_intergenic) / length(counts$abundance)
    
    ## Plot whole distribution
    title <- basename(myUserMetadata@rnaseq_lib_path)
    plot(
        dens,
        ylim = c(0, max(dens$y) * 1.1),
        xlim = c(-23, 21),
        lwd = 2,
        main = title,
        bty = "n",
        axes = FALSE,
        xlab = ""
    )
    axis(2, las = 1)
    
    axis(
        1,
        at = seq(-30 , 30, by = 10),
        line = 0,
        mgp = c(3, 0.5, 0),
        cex.axis = 0.8
    )
    mtext(
        expression(log[2]('TPM'+10 ^ -6)),
        1,
        adj = 1,
        padj = 0,
        line = 0.2,
        at = par("usr")[1],
        col = "black",
        cex = 0.8
    )
    
    
    ## Plot the TPM cutoff abline(v=cutoff, col='gray',
    ## lty=1, lwd=2)
    arrows(
        log2(cutoff + 1e-05),
        par("usr")[3],
        log2(cutoff + 1e-05),
        par("usr")[4] / 2,
        col = "gray",
        lty = 1,
        lwd = 2,
        angle = 160,
        length = 0.1
    )
    
    ## Add subgroups distributions (coding, intergenic,
    ## etc): protein-coding genes
    lines(dens_coding,
        col = "firebrick3",
        lwd = 2,
        lty = 2)
    ## intergenic
    lines(dens_intergenic,
        col = "dodgerblue3",
        lwd = 2,
        lty = 2)
    
    ## legend
    legend(
        "topright",
        c(
            "all",
            "selected protein-coding genes",
            "selected intergenic regions"
        ),
        lwd = 2,
        col = c("black",
                "firebrick3", "dodgerblue3"),
        lty = c(1, 2,
                2),
        bty = "n"
    )
    return()
}

# this function is not exactly the same than in the
# bgee pipeline.  We removed the max_intergenic
# variable. In the pipeline this variable was used
# to define the reference intergenic regions.  We
# do not need it anymore in this package as we
# precomputed the list of reference intergenic
# regions we also remove TPM_final_cutoff,
# FPKM_cutoff, FPKM_final_cutoff

#' @title Cutoff information
#'
#' @description Calculate summary statistics to export in cutoff info file
#'
#' @param counts TPM information for both genic and intergenic regions
#' @param column Name of the column containing presence/absence information
#' @param TPM_cutoff TPM cutoff below which calls are considered as absent
#' @param r_cutoff Proportion of intergenic regions considered as present
#' @param mean_pvalue mean information. Provided only when pvalue approach is used
#' @param sd_pvalue standard deviation. Provided only when pvalue approach is used 
#' @param myUserMetadata A Reference Class UserMetadata object.
#' @param myAbundanceMetadata A Reference Class AbundanceMetadata object.
#'
#' @return summary statistics to export in cutoff info file
#'
#' @author Julien Roux
#' @author Julien Wollbrett
#'
#' @noMd
#' @noRd
#'
cutoff_info <- function(counts, column, abundance_cutoff, r_cutoff, mean_pvalue=NULL, 
    sd_pvalue=NULL, myUserMetadata, myAbundanceMetadata) {
    ## Calculate summary statistics to export in cutoff
    ## info file
    genic_present <- sum(counts[[column]][counts$type ==
        "genic"] == "present") / sum(counts$type == "genic") * 100
    number_genic_present <- sum(counts[[column]][counts$type ==
        "genic"] == "present")
    
    coding_present <- sum(counts[[column]][counts$biotype %in%
        "protein_coding"] == "present") / sum(counts$biotype %in%
        "protein_coding") * 100
    number_coding_present <-
        sum(counts[[column]][counts$biotype %in%
        "protein_coding"] == "present")
    
    intergenic_present <- sum(counts[[column]][counts$type ==
        "intergenic"] == "present") / sum(counts$type ==
        "intergenic") * 100
    number_intergenic_present <-
        sum(counts[[column]][counts$type == "intergenic"] == "present")
    
    ## Export cutoff_info_file
    to_export <- c(
        basename(myUserMetadata@rnaseq_lib_path),
        abundance_cutoff,
        genic_present,
        number_genic_present,
        sum(counts$type == "genic"),
        coding_present,
        number_coding_present,
        sum(counts$biotype %in% "protein_coding"),
        intergenic_present,
        number_intergenic_present,
        sum(counts$type == "intergenic")
    )
    names(to_export) <- c(
        "libraryId",
        "cutoffTPM",
        "proportionGenicPresent",
        "numberGenicPresent",
        "numberGenic",
        "proportionCodingPresent",
        "numberPresentCoding",
        "numberCoding",
        "proportionIntergenicPresent",
        "numberIntergenicPresent",
        "numberIntergenic"
    )
    if(myAbundanceMetadata@cutoff_type == "intergenic") {
        to_export <- c(to_export, r_cutoff)
        names(to_export)[length(to_export)] <- "ratioIntergenicCodingPresent"
    } else if(myAbundanceMetadata@cutoff_type == "pValue") {
        to_export <- c(to_export, myAbundanceMetadata@cutoff, mean_pvalue, sd_pvalue)
        names(to_export)[length(to_export)-2] <- "pValueCutoff"
        names(to_export)[length(to_export)-1] <- "meanIntergenic"
        names(to_export)[length(to_export)] <- "sdIntergenic"
    } else if(myAbundanceMetadata@cutoff_type == "qValue") {
        to_export <- c(to_export, myAbundanceMetadata@cutoff)
        names(to_export)[length(to_export)] <- "qValueCutoff"
    } else {
        stop("unknown cutoff type : ", myAbundanceMetadata@cutoff_type, ". Should be 
        \"pValue\" or \"intergenic\" or \"qValue\"")
    }
    return(to_export)
}

#' @title Generate theoretical pValue
#'
#' @description Uses reference intergenic to calculate a pValue per gene id
#' Then uses the pValue cutoff provided by user to generate present/absent calls
#'
#' @param counts A list of estimated counts
#' @param pValueCutoff the pValue cutoff
#'
#' @return counts with zscore, pvalue and calls, but also the mean and the sd of ref. intergenic
#' 
#' @author Sara Fonseca Costa
#'
#' @noMd
#' @noRd
#' 
#' @import dplyr
#'
generate_theoretical_pValue <- function(counts, pValueCutoff) {
    ## select values with TPM > 0 (because we will use log2 scale)
    selected_intergenic <- filter(counts, abundance > 0 & type == "intergenic")

    ## select genic region from the library
    selected_count <- filter(counts, abundance > 0)  

    ## calculate z-score for each gene_id using the reference intergenic 
    selected_count$zScore <- (log2(selected_count$abundance) - mean(log2(selected_intergenic$abundance))) / sd(log2(selected_intergenic$abundance))
    ## calculate p-values for each gene_id
    selected_count$pValue <- pnorm(selected_count$zScore, lower.tail = FALSE)
    counts_with_pValue <- merge(counts, selected_count[, c("id", "zScore", "pValue")], 
                                by = "id", all.x=TRUE)
    
    counts_with_pValue$call <- ifelse((counts_with_pValue$pValue > pValueCutoff | 
                                           is.na(counts_with_pValue$pValue)), "absent", "present")
    mean <- 2^(mean(log2(selected_intergenic$abundance)))
    sd <- 2^(sd(log2(selected_intergenic$abundance)))
    return(list(counts_with_pValue = counts_with_pValue, mean = mean, sd = sd))
}

#' @title Generate qValue
#'
#' @description Calculate a qValue based on the ratio
#' intergenic/(intergenic+genic) for each gene id.
#' Then uses the qValue cutoff provided by user to generate present/absent calls
#'
#' @param counts A list of estimated counts
#' @param qValueCutoff the qValue cutoff
#'
#' @return counts with qvalue and calls
#' 
#' @author Sara Fonseca Costa
#'
#' @noMd
#' @noRd
#' 
#' @import dplyr
#'
generate_qValue <- function(counts, qValueCutoff) {
    
    ## select genic and intergenic regions
    selected_genic <- counts$type %in% "genic"
    selected_intergenic <- counts$type %in% "intergenic"
    
    ## collect the density distributions of genic and intergenic
    dens_genic <- density(log2(counts$abundance[selected_genic] + 10^-6))
    dens_intergenic <- density(log2(counts$abundance[selected_intergenic] + 10^-6))
    
    ## perform the linear interpolation for each type
    genicRegion <- approxfun(dens_genic$x, dens_genic$y)
    intergenicRegion <- approxfun(dens_intergenic$x, dens_intergenic$y)
    ## numerical integration
    numInt_genicRegion <- integrate(genicRegion, min(dens_genic$x), max(dens_genic$x), subdivisions=1000, rel.tol = .Machine$double.eps^0.01)$value
    numInt_intergenicRegion <- integrate(intergenicRegion, min(dens_intergenic$x), max(dens_intergenic$x), subdivisions=1000, rel.tol = .Machine$double.eps^0.01)$value
    
    ## for each abundance value (TPM) collect the genic and intergenic linear interpolation
    interpolationInfo <- c()
    for (i in 1:nrow(counts)) {
        log2TpmValue <- log2(counts$abundance[i])
        genicY <- genicRegion(log2TpmValue)
        intergenicY <- intergenicRegion(log2TpmValue)
        ## create an info table
        geneInfo <- c(log2TpmValue, genicY,intergenicY)
        interpolationInfo <- rbind(interpolationInfo, geneInfo)
    }
    interpolationInfo <- data.frame(interpolationInfo)
    colnames(interpolationInfo) <- c("log2TpmValue", "genicY","intergenicY")
    interpolationInfo$geneId <- counts$gene_id
    
    ## Function to calculate q-Value per unique TPM
    calculate_qValue <- function(log2Tpm, dens_genic, dens_intergenic, interpolationGenic, interpolationIntergenic, numInt_genic, numInt_intergenic){
        unscaled_genic <- integrate(interpolationGenic, log2Tpm, max(dens_genic$x), subdivisions=1000, stop.on.error = FALSE)$value
        scaled_genic <- unscaled_genic / numInt_genic
        unscaled_intergenic <- integrate(interpolationIntergenic, log2Tpm, max(dens_intergenic$x), subdivisions=1000, stop.on.error = FALSE)$value
        scaled_intergenic <- unscaled_intergenic / numInt_intergenic
        ## calculate qValue for target gene
        qValue <- scaled_intergenic / (scaled_intergenic + scaled_genic)
        return(qValue)
    }
    
    ## Use the calculate_qValue function to calculate the q-Value for each unique log2 TPM value and perform the calls
    qValueInfo <- c() 
    for (i in 1:nrow(interpolationInfo)) {
        
        log2TpmValue <- interpolationInfo$log2TpmValue[i]
        genicY_info <- interpolationInfo$genicY[i]
        intergenicY_info <- interpolationInfo$intergenicY[i]
        
        ## attribute qValue = 1 to inf log2TPM values and to values with Na to genicY_info and intergenicY_info (peak on left side of the plot)
        if ( log2TpmValue == "-Inf" | log2TpmValue != "-Inf" & is.na(genicY_info) == TRUE & is.na(intergenicY_info) == TRUE){
            qValue <- "1"
            qValueInfo <- rbind(qValueInfo,qValue)
            
        } else if (log2TpmValue != "-Inf" & (is.na(intergenicY_info) == TRUE & is.na(genicY_info) == FALSE) | (is.na(intergenicY_info) == FALSE & is.na(genicY_info) == TRUE)){
            
            ## retrieve log2TPM value where was possible to calculate the linear interpolation for genic and intergenic
            maxNumIntegration <- filter(interpolationInfo, genicY != "NaN" & intergenicY != "NaN" )
            
            ## if log2TPM value is a negative value (means peak more at left side)
            ## we attribute qValue based on the minimum where was possible to calculate the linear interpolation for both of the types
            if(log2TpmValue < 0){
                
                log2TpmValue <- min(maxNumIntegration$log2TpmValue)
                qValue <- calculate_qValue(log2TpmValue, dens_genic, dens_intergenic, genicRegion,
                                           intergenicRegion, numInt_genicRegion, numInt_intergenicRegion)
                qValueInfo <- rbind(qValueInfo,qValue)       
                
            } else {
                
                ## if log2TPM value is positive, this represent the values in the tail on the right side of the plot.
                ## we attribute qValue for this cases based on the last log2TPM value where was possible to calculate the linear interpolation for both of the types
                log2TpmValue <- max(maxNumIntegration$log2TpmValue)
                qValue <- calculate_qValue(log2TpmValue, dens_genic, dens_intergenic, genicRegion,
                                           intergenicRegion, numInt_genicRegion, numInt_intergenicRegion)
                qValueInfo <- rbind(qValueInfo,qValue)    
            }
        } else {
            
            ## calculate q-value using the linear interpolation info for a particular log2TPM value
            qValue <- calculate_qValue(log2TpmValue, dens_genic, dens_intergenic, genicRegion,
                                       intergenicRegion, numInt_genicRegion, numInt_intergenicRegion)
            qValueInfo <- rbind(qValueInfo,qValue)     
        }
    }
    ## add qValue columns to the final table
    counts$qValue <- as.numeric(qValueInfo[,1])
    counts$call <- ifelse(counts$qValue <= qValueCutoff , "present", "absent")
    return(counts)
}

#' @title Recalculate TPMs
#'
#' @description Recalculate TPMs using functions from
#' ```https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/```
#'
#' @param counts A list of estimated counts
#' @param effLen A list of effective length
#'
#' @return A list of recalculated TPMs
#'
#' @noMd
#' @noRd
#'
countToTpm <- function(counts, effLen) {
    rate <- log(counts) - log(effLen)
    denom <- log(sum(exp(rate)))
    exp(rate - denom + log(1e+06))
}
