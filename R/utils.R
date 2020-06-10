# Usefull functions potentially used everywhere in
# the package.  These functions are only used
# inside of the package. They are not exported.

#' @title get Operating System
#'

#' @description Function used to detect the OS in which the package is run.
#' Return 'linux', 'osx', or 'windows', depending on the OS
#'
#' @noMd
#' @noRd
#'
get_os <- function() {
    sysinf <- Sys.info()
    if (!is.null(sysinf)) {
        os <- sysinf["sysname"]
        if (os == "Darwin") {
            os <- "osx"
        } else {
            os <- .Platform$OS.type
            if (os == "unix") {
                if (grepl("^darwin", R.version$os)) {
                    os <- "osx"
                }
                if (grepl("linux-gnu", R.version$os)) {
                    os <- "linux"
                }
            } else if (os != "windows") {
                stop(paste0("Unrecognized Operating System : ",
                            os, "!!\n"))
            }
        }
        return(tolower(os))
    }
    stop("Sys.info() function returned a null value")
}

#' @title Path to bgee release directory
#'
#' @description helper function to get the path to the bgee release directory
#'
#' @param myBgeeMetadata A Reference Class BgeeMetadata object.
#' @param myUserMetadata A Reference Class UserMetadata object.
#'
#' @return path to the file containing details about intergenic releases
#'
#' @noMd
#' @noRd
#'
get_intergenic_release_path <- function(myBgeeMetadata,
                                        myUserMetadata) {
    return(file.path(
        myUserMetadata@working_path,
        paste0(
            myBgeeMetadata@intergenic_prefix,
            myBgeeMetadata@intergenic_release
        )
    ))
}

#' @title Path to species directory
#'
#' @description helper function to get the path to the species directory
#'
#' @noMd
#' @noRd
#'
get_species_path <- function(myBgeeMetadata, myUserMetadata) {
    if (nchar(myUserMetadata@species_id) == 0) {
        stop("the object of the UserMetadata class must contains a species_id")
    }
    return(file.path(
        get_intergenic_release_path(myBgeeMetadata,
                                    myUserMetadata),
        myUserMetadata@species_id
    ))
}

#' @title Path to transcriptome directory
#'
#' @description helper function to get the path to one transcriptome with
#' intergenic regions of one species
#'
#' @noMd
#' @noRd
#'
get_transcriptome_path <- function(myBgeeMetadata,
    myUserMetadata) {
    return(file.path(
        get_species_path(myBgeeMetadata,
            myUserMetadata),
        paste0(
            "transcriptome_",
            gsub("\\.",
            "_", myUserMetadata@transcriptome_name)
        )
    ))
}

#' @title Path to annotation directory
#'
#' @description helper function to get the path to annotation directory of one
#' species. This annotation directory will contain both gene2biotype and tx2gene
#' files.
#'
#' @noMd
#' @noRd
#'
get_annotation_path <- function(myBgeeMetadata, myUserMetadata) {
    return(file.path(
        get_species_path(myBgeeMetadata,
            myUserMetadata),
        paste0(
            "annotation_",
            gsub("\\.",
            "_", myUserMetadata@annotation_name)
        )
    ))
}

#' @title Path to abundance tool directory
#'
#' @description helper function to get the path to the abundance tool directory
#'
#' @noMd
#' @noRd
#'
get_tool_path <- function(myAbundanceMetadata,
    myBgeeMetadata,
    myUserMetadata) {
    return(file.path(
        get_species_path(myBgeeMetadata,
            myUserMetadata),
        myAbundanceMetadata@tool_name
    ))
}

#' @title Path to the index created directory
#'
#' @description helper function to get the path to transcriptome specific
#' directory of one abundance tool.
#'
#' @noMd
#' @noRd
#'
get_tool_transcriptome_path <- function(myAbundanceMetadata,
    myBgeeMetadata,
    myUserMetadata) {
    return(file.path(
        get_tool_path(myAbundanceMetadata,
            myBgeeMetadata, myUserMetadata),
        paste0(
            "transcriptome_",
            gsub("\\.", "_", myUserMetadata@transcriptome_name)
        )
    ))
}

#' @title Path to the output directory of one abundance tool
#'
#' @description helper function to get the path to the output directory of
#' one abundance tool.
#'
#' @noMd
#' @noRd
#'
get_tool_output_path <- function(myAbundanceMetadata,
    myBgeeMetadata,
    myUserMetadata) {
    if (is.na(myUserMetadata@output_dir) ||
        myUserMetadata@output_dir == "") {
        if (myUserMetadata@simple_arborescence) {
            return(
                file.path(
                    get_intergenic_release_path(myBgeeMetadata,
                        myUserMetadata),
                    "all_results",
                    get_output_dir(myUserMetadata)
                )
            )
        }
        return(file.path(
            get_tool_transcriptome_path(myAbundanceMetadata,
                myBgeeMetadata, myUserMetadata),
            paste0(
                "annotation_",
                gsub("\\.", "_", myUserMetadata@annotation_name)
            ),
            get_output_dir(myUserMetadata)
        ))
    } else {
        return(myUserMetadata@output_dir)
    }
}

#' @title Download fasta intergenic
#'
#' @description Check if reference intergenic fasta file has already
#' been downloaded. If not the file is downloaded from Bgee FTP or from
#' the community repository depending on myBgeeMetadata@intergenic_release.
#' if myBgeeMetadata@intergenic_release == "community" then reference intergenic
#' wil be downloaded from the Zenodo community repository. Otherwise Reference
#' intergenic sequences will be downloaded from the official Bgee FTP.
#' Be careful when using reference intergenic sequences generated by the community
#' as the Bgee team do not deeply review them.
#'
#'
#' @param myBgeeMetadata A Reference Class BgeeMetadata object (optional)
#' @param myUserMetadata A Reference Class UserMetadata object.
#' @param intergenic_file path where intergenic file will be saved
#'
#' @export
#'
#' @examples {
#' bgee_intergenic_file <- file.path(getwd(), 'intergenic.fasta')
#' userMetadata <- new('UserMetadata', species_id = '7227')
#' }
#'
#' @return download fasta intergenic from Bgee FTP or from the Zenodo community and save it
#' locally
#'
download_fasta_intergenic <-
    function(myBgeeMetadata = new("BgeeMetadata"),
        myUserMetadata,
        intergenic_file) {
        if (myBgeeMetadata@intergenic_release == "community") {
            intergenic_url <-
                retrieve_community_ref_intergenic_url(myUserMetadata@species_id)
        } else {
            all_releases <- myBgeeMetadata@all_releases
            intergenic_url <- as.character(all_releases$referenceIntergenicFastaURL[
                all_releases$release == myBgeeMetadata@intergenic_release])
            intergenic_url <-
                gsub("SPECIES_ID",
                myUserMetadata@species_id,
                intergenic_url)
        }
        success <- download.file(url = intergenic_url,
            destfile = intergenic_file)
        if (success != 0) {
            stop(
                "ERROR: Downloading reference intergenic sequences from FTP was not successful."
            )
        }
    }

#' @title Retireve name of fastq files
#'
#' @description Check the presence of fastq files in the fastq directory.
#' * If no myUserMetadata@run_ids are provided it will find all distinct run ids.
#'   These run will be considered as technical replicates and merged together.
#'   With technical replicates it is not possible to have a combination of
#'   single-end AND paired-end run. Then if the first run is detected as paired-end
#'   (presence of 2 fastq files with names finishing with _1 and _2) all fastq
#'   files will be considered as paired-end fastq files. If a mix of single_end
#'   and paired_end run are detected the function will return an error
#' * If myUserMetadata@run_ids are provided they will be considered as technical
#'   replicates and run in the same transcript expression estimation analyse.
#'   For the same reason than described in previous section, it is not possible
#'   to combine single-end and paired-end runs.
#'
#' @return A character containing the name of all fastq files
#'
#' @noMd
#' @noRd
#'
get_merged_fastq_file_names <- function(myUserMetadata) {
    message("debug problem with fastq files... rnaseq lib path: ", myUserMetadata@rnaseq_lib_path)
    fastq_files <- get_fastq_files(myUserMetadata)
    message("debug problem with fastq files... all fastq_files: ", fastq_files)
    
    
    # filter list of fastq files if run_ids are
    # provided\024
    if (length(myUserMetadata@run_ids) != 0) {
        fastq_files <- unique(grep(
            paste(myUserMetadata@run_ids,
                collapse = "|"),
            fastq_files,
            value = TRUE
        ))
    }
    fastq_files_names <- ""
    if (is_pair_end(fastq_files)) {
        first_files <- sort(grep("_1", fastq_files,
            value = TRUE))
        second_files <- sort(grep("_2", fastq_files,
            value = TRUE))
        if (length(first_files) != length(second_files)) {
            stop(
                paste0(
                    "Can not run a paired-end expression estimation if not same
                    number of file finishing with _1 and _2... In library ",
                    basename(myUserMetadata@rnaseq_lib_path)
                )
            )
        }
        if (length(first_files) + length(second_files) !=
            length(fastq_files)) {
            stop(
                paste0(
                    "Can not run a paired-end expression estimation if not all
                    fastq file names end with _1 or _2... In library ",
                    basename(myUserMetadata@rnaseq_lib_path)
                )
            )
        }
        for (i in seq_along(first_files)) {
            run_1 <- sub("^([^_]+).*", "\\1", first_files[i],
                perl = TRUE)
            run_2 <- sub("^([^_]+).*", "\\1", second_files[i],
                perl = TRUE)
            if (run_1 == run_2) {
                # combine all fastq_files in a character like A_1
                # A_2 B_1 B_2 ...
                fastq_files_names = paste(
                    fastq_files_names,
                    file.path(
                        myUserMetadata@rnaseq_lib_path,
                        first_files[i]
                    ),
                    file.path(
                        myUserMetadata@rnaseq_lib_path,
                        second_files[i]
                    ),
                    sep = " "
                )
            }
        }
    } else {
        if (length(grep("_1", fastq_files, value = TRUE)) !=
            0 || length(grep("_2", fastq_files, value = TRUE)) !=
            0) {
            stop(
                paste0(
                    "Looks like a combination of single-end and paired-end
                    (file name end with _1 or _2) fastq files for library ",
                    basename(myUserMetadata@rnaseq_lib_path),
                    ".\n"
                )
            )
        }
        for (i in seq_along(fastq_files)) {
            fastq_files_names = paste(
                fastq_files_names,
                file.path(myUserMetadata@rnaseq_lib_path,
                    fastq_files[i]),
                sep = " "
            )
        }
    }
    return(fastq_files_names)
}

#' @title get fastq files
#'
#' @description retrieve all fastq files names present in the
#' myUserMetadata@rnaseq_lib_path directory
#'
#' @noMd
#' @noRd
#'
get_fastq_files <- function(myUserMetadata) {
    # all files of the library directory
    library_files <-
        list.files(path = myUserMetadata@rnaseq_lib_path)
    fastq_files <- ""
    i <- 1
    for (library_file in library_files) {
        if (grepl(".fq$", library_file) || grepl(".fq.gz$",
            library_file) || grepl(".fastq.gz$", library_file) ||
            grepl(".fastq$", library_file)) {
            fastq_files[i] <- library_file
            i <- i + 1
        }
    }
    return(fastq_files)
}

#' @title is paired-end
#'
#' @description check is the first element of a vector of fastq files names
#' correspond to a paired-end file. Return a boolean
#'
#' @noMd
#' @noRd
#'
is_pair_end <- function(fastq_files) {
    return(grepl("_1.", fastq_files[1]) || grepl("_2.",
        fastq_files[1]))
}

#' @title Retrieve name of output directory
#'
#' @description Retireve name of output directory depending on
#' myUserMetadata@rnaseq_lib_path and myUserMetadata@run_ids
#' if myUserMetadata@run_ids is empty then the output directory will correspond
#' to name of the last directory of myUserMetadata@rnaseq_lib_path.
#' Otherwise the name of the output directory will be the concatenation of the
#' name of the last directory of myUserMetadata@rnaseq_lib_path and all
#' myUserMetadata@run_ids
#'
#' @param myUserMetadata Reference Class UserMetadata object.
#'
#' @return Name of the output directory
#'
#' @noMd
#' @noRd
#'
get_output_dir <- function(myUserMetadata) {
    if (length(myUserMetadata@rnaseq_lib_path) == 0) {
        stop(
            "No fastq path provided. Please edit `rnaseq_lib_path`
            attribute of UserMetadata class"
        )
    }
    if (length(myUserMetadata@run_ids) == 0) {
        return(basename(myUserMetadata@rnaseq_lib_path))
    } else {
        return(paste0(
            basename(myUserMetadata@rnaseq_lib_path),
            "_",
            paste(myUserMetadata@run_ids, collapse = "_")
        ))
    }
}

#' @title Retrieve path to kallisto directory
#'
#' @description Retireve path to kallisto directory. This path depends
#' on the OS
#'
#' @param myAbundanceMetadata Reference Class AbundanceMetadata object.
#' @param myUserMetadata Reference Class UserMetadata object.
#'
#' @return path to kallisto dir
#'
#' @noMd
#' @noRd
#'
get_kallisto_dir_path <- function(myAbundanceMetadata,
    myUserMetadata) {
    os_version <- get_os()
    if (os_version == "linux") {
        return(
            file.path(
                myUserMetadata@working_path,
                myAbundanceMetadata@kallisto_linux_dir
            )
        )
    } else if (os_version == "osx") {
        return(
            file.path(
                myUserMetadata@working_path,
                myAbundanceMetadata@kallisto_osx_dir
            )
        )
    } else if (os_version == "windows") {
        return(
            file.path(
                myUserMetadata@working_path,
                myAbundanceMetadata@kallisto_windows_dir
            )
        )
    } else {
        stop(
            "can not access to kallisto dir for this operating system.
            Please install it by yourself and change the value of the
            parameter `download_kallisto` of the AbundanceMetadata
            object to FALSE"
        )
    }
}

#' @title Retrieve path to kallisto program
#'
#' @description Retireve path to kallisto program file. This path depends
#' on the OS
#'
#' @param myAbundanceMetadata Reference Class AbundanceMetadata object.
#' @param myUserMetadata Reference Class UserMetadata object.
#'
#' @return path to kallisto program
#'
#' @noMd
#' @noRd
#'
get_kallisto_program_path <- function(myAbundanceMetadata,
    myUserMetadata) {
    if (!myAbundanceMetadata@download_kallisto) {
        return(myAbundanceMetadata@tool_name)
    }
    os_version <- get_os()
    kallisto_dir <- get_kallisto_dir_path(myAbundanceMetadata,
        myUserMetadata)
    if (os_version == "linux" || os_version == "osx") {
        return(file.path(
            kallisto_dir,
            myAbundanceMetadata@unix_kallisto_name
        ))
    } else if (os_version == "windows") {
        return(file.path(
            kallisto_dir,
            myAbundanceMetadata@windows_kallisto_name
        ))
    } else {
        stop("can not access to kallisto program file for Your operating system.")
    }
}

#' @title Remove transcript version from abundance file
#'
#' @description Remove the transcript version that can be present in transcript
#' id of transcriptome files.
#' Removing the transcript version means detecting a dot in transcript id and
#' removing the dot and all caracters following it. This function will not affect
#' the transcript ID of intergenic regions because the name these intergenic
#' sequences can contain dot that do not correspond to transcript version (e.g
#' KN149689.1_73093_93092).
#' This function has been develop because the ignoreTxVersion parameter of tximport
#' do not allow to remove the dot only in a subset of transcript IDs.
#' It is called if the logical 'ignoreTxVersion' attribut of the AbundanceMetadata
#' class is set to TRUE.
#'
#' @noMd
#' @noRd
#'
removeTxVersionFromAbundance <- function(myAbundanceMetadata,
    myBgeeMetadata,
    myUserMetadata) {
    output_path <- get_tool_output_path(myAbundanceMetadata,
        myBgeeMetadata, myUserMetadata)
    abundance_path <-
        file.path(output_path, myAbundanceMetadata@abundance_file)
    abundance_path_without_tx_version <-
        file.path(output_path,
            myAbundanceMetadata@abundance_file_without_tx_version)
    abundance_data <- read.table(file = abundance_path,
        header = TRUE, sep = "\t")
    intergenic_ids <-
        get_intergenic_ids(myBgeeMetadata, myUserMetadata)
    if (myAbundanceMetadata@tool_name == "kallisto") {
        ref_intergenic_transcripts <-
            abundance_data[abundance_data$target_id
                %in% intergenic_ids$intergenic_ids, ]
        gtf_transcripts <-
            abundance_data[!(abundance_data$target_id
                %in% intergenic_ids$intergenic_ids), ]
        gtf_transcripts$target_id <- gsub(pattern = "\\..*",
            "", gtf_transcripts$target_id)
        abundance_data <-
            rbind(gtf_transcripts, ref_intergenic_transcripts)
    } else {
        stop(
            paste0(
                "Removing transcript version for tool ",
                myAbundanceMetadata$tool_name,
                " is not implemented."
            )
        )
    }
    write.table(
        x = abundance_data,
        file = abundance_path_without_tx_version,
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE
    )
}

#' Generate a summary of the Slots of the S4 objects used to
#' generate present / absent calls
#'
#' @noMd
#' @noRd
#'
generate_S4_object_properties_output <-
    function(myAbundanceMetadata,
        myBgeeMetadata,
        myUserMetadata) {
        output <- matrix(nrow = 13, ncol = 2)
        colnames(output) <- c("Slot name", "Slot value")
        output[1, ] <- c("AbundanceMetadata@tool_name",
            myAbundanceMetadata@tool_name)
        output[2, ] <- c("AbundanceMetadata@txOut",
            myAbundanceMetadata@txOut)
        output[3, ] <- c("AbundanceMetadata@ignoreTxVersion",
            myAbundanceMetadata@ignoreTxVersion)
        output[4, ] <- c("AbundanceMetadata@cutoff",
            myAbundanceMetadata@cutoff)
        output[5, ] <- c(
            "AbundanceMetadata@read_size_kmer_threshold",
            myAbundanceMetadata@read_size_kmer_threshold
        )
        output[6, ] <- c("BgeeMetadata@intergenic_release",
            myBgeeMetadata@intergenic_release)
        output[7, ] <- c("UserMetadata@species_id",
            myUserMetadata@species_id)
        output[8, ] <- c("UserMetadata@reads_size",
            myUserMetadata@reads_size)
        output[9, ] <- c("UserMetadata@rnaseq_lib_path",
            myUserMetadata@rnaseq_lib_path)
        output[10, ] <- c("UserMetadata@transcriptome_name",
            myUserMetadata@transcriptome_name)
        output[11, ] <- c("UserMetadata@annotation_name",
            myUserMetadata@annotation_name)
        output[12, ] <- c("UserMetadata@simple_arborescence",
            myUserMetadata@simple_arborescence)
        output[13, ] <- c(
            "output_dir",
            get_tool_output_path(myAbundanceMetadata,
                myBgeeMetadata, myUserMetadata)
        )
        return(output)
    }

#' potentially download reference intergenic sequences and retrieve path to the
#' corresponding fasta file
#'
#' @noMd
#' @noRd
#'
retrieve_intergenic_path <-
    function(myBgeeMetadata, myUserMetadata) {
        bgee_intergenic_file <- file.path(
            get_species_path(myBgeeMetadata,
                myUserMetadata),
            myBgeeMetadata@fasta_intergenic_name
        )
        if (!file.exists(bgee_intergenic_file)) {
            # Check if custom reference intergenic path has to be used
            if (!(myBgeeMetadata@custom_intergenic_path == "")) {
                if(!file.exists(myBgeeMetadata@custom_intergenic_path)) {
                    stop("File ", myBgeeMetadata@custom_intergenic_path, 
                         " selected as custom intergenic does not exist")
                }
                if (myBgeeMetadata@intergenic_release != "custom") {
                    stop(
                        "You selected a custom intergenic path (BgeeMetadata@custom_intergenic_path)
                    but the intergenic release (BgeeMetadata@intergenic_release) was not defined
                    as `custom`. In order to use custom reference intergenic sequences please
                    provide the path to your custom_intergenic_path AND update the
                    intergenic_release to `custom`."
                    )
                }
                bgee_intergenic_file <-
                    myBgeeMetadata@custom_intergenic_path
                return(bgee_intergenic_file)
            } else {
                if (myBgeeMetadata@intergenic_release == "custom") {
                    stop(
                        "You selected a `custom`` intergenic release (BgeeMetadata@intergenic_release)
                    but the intergenic path (BgeeMetadata@custom_intergenic_path) was not defined. In
                    order to use custom reference intergenic sequences please both provide the path
                    to your custom_intergenic_path AND update the intergenic_release to `custom`."
                    )
                }
            }
            if (!dir.exists(dirname(bgee_intergenic_file))) {
                dir.create(dirname(bgee_intergenic_file), recursive = TRUE)
            }
            download_fasta_intergenic(myBgeeMetadata,
                myUserMetadata, bgee_intergenic_file)
        }
        return(bgee_intergenic_file)
    }

#' @title Retireve intergenic IDs
#'
#' @description Check if an intergenic IDs file is already present for the
#' wanted intergenic_release and species ID.
#' - if not create the file
#' - else read the file
#'
#' @return A data frame containing the ID of all intergenic sequences
#'
#' @noMd
#' @noRd
#'
get_intergenic_ids <- function(myBgeeMetadata, myUserMetadata) {
    species_path <- get_species_path(myBgeeMetadata, myUserMetadata)
    intergenic_ids_file_name <- "intergenic_ids.txt"
    intergenic_ids_file <-
        file.path(species_path, intergenic_ids_file_name)
    # generate file if it does not exist
    intergenic_ids <- ""
    if (!file.exists(intergenic_ids_file)) {
        bgee_intergenic_file <-
            retrieve_intergenic_path(myBgeeMetadata, myUserMetadata)
        bgee_intergenic <- readDNAStringSet(bgee_intergenic_file)
        # intergenic ID correspond to part of the header
        # before the first space character
        intergenic_ids <- as.data.frame(sub("^([^ ]+) .*",
                                            "\\1", names(bgee_intergenic)))
        colnames(intergenic_ids) <- "intergenic_ids"
        write.table(
            x = intergenic_ids,
            file = intergenic_ids_file,
            quote = FALSE,
            row.names = TRUE
        )
    } else {
        intergenic_ids <-
            read.table(file = intergenic_ids_file, header = TRUE)
    }
    return(intergenic_ids)
}

get_abundance_file_path <- function(myAbundanceMetadata,
    myBgeeMetadata,
    myUserMetadata) {
    output_path <- get_tool_output_path(myAbundanceMetadata,
        myBgeeMetadata, myUserMetadata)
    if (myAbundanceMetadata@ignoreTxVersion) {
        abundance_file <- file.path(output_path,
            myAbundanceMetadata@abundance_file_without_tx_version)
        if (!file.exists(abundance_file)) {
            removeTxVersionFromAbundance(myAbundanceMetadata,
                myBgeeMetadata,
                myUserMetadata)
        }
    } else {
        abundance_file <-
            file.path(output_path, myAbundanceMetadata@abundance_file)
    }
    return(abundance_file)
}

quiet <- function(x) { 
    sink(tempfile()) 
    on.exit(sink()) 
    invisible(force(x)) 
} 