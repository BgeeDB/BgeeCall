# BgeeCall, a R package for automatic RNA-Seq present/absent gene expression calls generation
[![Bioc](http://www.bioconductor.org/shields/years-in-bioc/BgeeCall.svg)](https://www.bioconductor.org/packages/devel/bioc/html/BgeeCall.html#since)
[![Bluesky](https://img.shields.io/badge/-Bluesky-3686f7?style=social&label=Follow%20bgeedb&logo=bluesky&logoColor=blue&labelColor=white&domain=https%3A%2F%2Fbsky.app)](https://bsky.app/profile/bgee.org)
[![Mastodon](https://img.shields.io/mastodon/follow/109308703977124988?style=social&label=Follow%20%40bgeedb&domain=https%3A%2F%2Fgenomic.social)](https://genomic.social/%40bgeedb)

`BgeeCall` is a collection of functions that uses [Bgee](https://bgee.org/) expertise to create gene expression present/absent calls

The `BgeeCall` package allows to: 

* Generate calls of presence/absence of expression at the gene level. You can generate these calls for any RNA-Seq samples as long as the species is available in Bgee (information available from the `list_bgee_species()` function).
* Download reference intergenic sequences for species available in Bgee.
* Generate calls of presence/absence of expression at the transcript level (beta version)

If you find a bug or have any issues with `BgeeCall` please write a bug report in our [GitHub issues manager](https://github.com/BgeeDB/BgeeCall/issues).

## How present/absent calls are generated

In Bgee present/absent gene expression calls for RNA-seq are generated using a threshold specific of each RNA-Seq library, calculated using reads mapped to reference intergenic regions. This is unlike the more usual use of an arbitrary threshold below which a gene is not considered as expressed (e.g log2(TPM) = 1). 

### Bgee database
Bgee is a database to retrieve and compare gene expression patterns in multiple animal species, produced from multiple data types (RNA-Seq, Affymetrix, in situ hybridization, and EST data). It notably integrates RNA-Seq libraries for 29 species. 

### Reference intergenic regions
Reference intergenic regions are defined in the [Bgee RNA-Seq pipeline](https://github.com/BgeeDB/bgee_pipeline/tree/master/pipeline/RNA_Seq).
Candidate intergenic regions are defined using gene annotation data. For each species, over all available libraries, reads are mapped to these intergenic regions with [kallisto](https://github.com/pachterlab/kallisto), as well as to genes. This "intergenic expression" is deconvoluted to distinguish reference intergenic from non annotated genes, which have higher expression. Reference intergenic regions are then defined as intergenic regions with low expression level over all RNA-Seq libraries, relative to genes. This step allows not to consider regions wrongly considered as intergenic because of potential gene annotation quality problem as intergenic. For more information please refer to the [Bgee RNA-Seq pipeline](https://github.com/BgeeDB/bgee_pipeline/tree/master/pipeline/RNA_Seq).

### Threshold of present/absent
BgeeCall pipeline allows to download reference intergenic regions resulting from the expertise of the Bgee team.
Moreover BgeeCall allows to use these reference intergenic regions to automatically generate gene expression calls for your own RNA-Seq libraries as long as the species is integrated to [Bgee](https://bgee.org/)
The present/absent abundance threshold is calculated for each library using the formula :

            proportion of ref intergenic present / proportion of protein coding present = 0.05

## Installation
In R:
``` {r, message = FALSE, warning = FALSE}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("BgeeCall")
```

## How to use the BgeeCall package

BgeeCall is highly tunable. Do not hesitate to have a look at the reference manual to have a precise descripton of all slots of the 4 main S4 classes (AbundanceMetadata, KallistoMetadata, BgeeMetadata and UserMetadata) or of all available functions. If you do not find any answer to your question do not hesitate to contact us.


### Load the package
``` {r, message = FALSE, warning = FALSE}
library(BgeeCall)
```


### Quick start

With the BgeeCall package it is easy to generate present/absent gene expression calls.
The most time comsuming task of this calls generation is the generation of the kallisto transcriptome index.
As the time needed for this step depend on the size of the transcriptome, we choose C. elegans as an example, because it is the smallest transcriptome file among all species available on Bgee.
To generate these calls you will need :

- a transcriptome
- genome annotations
- your RNA-Seq reads in fastq files

For this vignette we created a toy fastq file example based on the SRX099901 library using the [ShortRead](http://bioconductor.org/packages/release/bioc/html/ShortRead.html) R package

``` {r, eval=FALSE}
library("ShortRead")
# keep 48.000 reads
sampler <- FastqSampler(file.path("absolute_path","/SRX099901/SRR350955.fastq.gz"), 48000)
set.seed(1); SRR350955 <- yield(sampler)
writeFastq(object = SRR350955, 
          file =file.path( "absolute_path","SRX099901_subset", "SRR350955_subset.fastq.gz"),
          mode = "w", full = FALSE, compress = TRUE)
```

In this example we used the Bioconductor AnnotationHub to load transcriptome and gene annotations but you can load them from wherever you want.
``` {r, message = FALSE, warning = FALSE}
ah <- AnnotationHub()
ah_resources <- query(ah, c("Ensembl", "Caenorhabditis elegans", "84"))
annotation_object <- ah_resources[["AH50789"]]
transcriptome_object <- rtracklayer::import.2bit(ah_resources[["AH50453"]])
```

Once you have access to transcriptome, gene annotations and your RNA-Seq library, an object of class `UserMetadata` has to be created.
``` {r, message = FALSE, warning = FALSE}
# create an object of class UserMetadata and specify the species ID
user_BgeeCall <- new("UserMetadata", species_id = "6239")
# import annotation and transcriptome in the user_BgeeCall object
# it is possible to import them using an S4 object (GRanges, DNAStringSet) or a file (gtf, fasta)
user_BgeeCall <- setAnnotationFromObject(user_BgeeCall, annotation_object, "WBcel235_84")
user_BgeeCall <- setTranscriptomeFromObject(user_BgeeCall, transcriptome_object, "WBcel235")
# provide path to the directory of your RNA-Seq library containing all the fastq files
user_BgeeCall <- setRNASeqLibPath(user_BgeeCall, 
                                  system.file("extdata", "SRX099901_subset", package = "BgeeCall"))
```

And that's it... You can run the generation of your present/absent gene expression calls
``` {r, message = FALSE, warning = FALSE}
calls_output <- run_from_object(myUserMetadata = user_BgeeCall)
```
Each analyze generates 5 files and return path to each one of them.

* calls_tsv_path : path to main tsv file with TPM, count, length, biotype, type, and presence/absence of expression summarized at gene level (or at the [transcript level](#transcript_level) if it was requested)
``` {r, message = FALSE, warning = FALSE}
head(read.table(calls_output$calls_tsv_path, header = TRUE), n = 5)
```
* cutoff_info_file_path : path to tsv summary of the analyze containing the proportion of gene, protein coding and intergenic defined as expressed. It also contains the library ID and the present/absent TPM threshold
``` {r, message = FALSE, warning = FALSE}
read.table(calls_output$cutoff_info_file_path)
```
* abundance_tsv : path to tsv kallisto quant output file
``` {r, message = FALSE, warning = FALSE}
head(read.table(calls_output$abundance_tsv, header = TRUE), n = 5)
```
* TPM_distribution_path : path to plot in pdf reprensting density distribution of TPM values for all sequences, protein coding sequences, and intergenic sequences. The grey line corresponds to TPM threshold used to generate present/absent calls.
``` {r, eval = FALSE}
openPDF(calls_output$TPM_distribution_path)
```
* S4_slots_summary : path to tsv file containing a summary of values used for the most important slots of the three S4 classes (UserMetadata, KallistoMetadata, and BgeeMetadata).
``` {r, message = FALSE, warning = FALSE}
read.table(S4_slots_summary, header = TRUE)
```

###  Generate present/absent calls for more than one RNA-Seq library

You will potentialy be also interested to generate present/absent calls on different RNA-Seq libraries, potentially on different species, or using  The main function `generate_presence_absence()` allows to generate present/absent calls from a UserMetadata object but also from a data frame or a tsv file depending on the arguments of the function you use. Please choose one of the three following arguments :
- userMetadata : Allows to generate present/absent calls for one RNA-Seq library using one object of the class UserMetadata.  
- userDataFrame : Provide a dataframe where each row correspond to one present/absent call generation. It allows to generate present/absent calls on different libraries, species, transcriptome, genome annotations, etc.
- userFile : Similar to userDataFrame except that the information are stored in a tsv file. A template of this file called `userMetadataTemplate.tsv` is available at the root of the package.

Columns of the dataframe or the tsv file are :

- species_id : The NCBI ID of the species.
- run_ids : The runs of the RNA-Seq library you want to use for the generation of the calls.	Allows to generate expression calls for a subset of the runs of one RNA-Seq library as described in [Generate calls for a subset of RNA-Seq runs](#run_ids). If not interested by this option, leave the column empty.
- reads_size : The size of the reads of your RNA-Seq library.
- rnaseq_lib_path : Path to the directory containing all fastq files generated for this library. This directory can only contains single-end runs or paired-end runs.
- transcriptome_path : path to the transcriptome file.
- annotation_path : path to the genome annotation file. Works with GTF of GFF3 files.
- working_path : path to the working directory where results will be stored. Using the same working directory for different RNA-Seq libraries of the same species will allow to reuse previously generated data like the custom transcriptome index (generated from both transcriptome and reference intergenic sequences). By default the working path is defined by the `getwd()` function and correspond to the working directory of your R session. If not interested by this option, leave the column empty.
- output_directory : Both species results and RNA-Seq libraries results are by default stored at the same place using the value of the `working_path` column. However, this column allows you to define a different output_directory for RNA-Seq results. For instance it allows you to save calls information directly in the RNA-Seq directory. If not interested by this option, leave the column empty.

Once the file has been fill in expression calls can be generated with :
``` {r, eval=FALSE}
calls_output <- generate_calls_workflow(userFile = "path_to_your_file.tsv")
```

### Parallized generation of present/absent calls on a cluster

BgeeCall already implement everything you need to generate calls on a cluster if it uses slurm as queuing system.
The same TSV file as described in the previous section will be used as input. In addition to tuning options available when running
BgeeCall on your computer, It is possible to modify how slurm jobs are submitted. More information available in section [Modify slurm options](#slurm_options)
In order to optimize parallelization calls will be generated in 2 steps.  
\itemize {
  \item Generate data at species level (e.g trancriptome with intergenic sequences, kallisto indexes)
  \item Generate expression calls for each RNA-Seq libraries
``` {r, eval = FALSE}
# generate kallisto indexes
generate_slurm_indexes(userFile = "path_to_your_file.tsv")
# generate expression calls
generate_slurm_calls(userFile = "path_to_your_file.tsv")
```

### Reference intergenic sequences

#### Releases of reference intergenic sequences

Different releases of reference intergenic sequences are available. It is possible to list all these releases :
```{r}
list_intergenic_release()
```
It is then possible to choose one specific release to create a `BgeeMetadata` object. Always use the setter method ` setIntergenicRelease()` when changing the release of an already existing BgeeMetadata object.
```{r}
# create BgeeMetadata object and define one reference intergenic release
bgee <- new("BgeeMetadata", intergenic_release = "0.1")
# change the reference intergenic release of your BgeeMetadata object
bgee <- setIntergenicRelease(bgee, "0.2")
```
By default the reference intergenic release used when a `BgeeMetadata` object is created is the last stable one created by the Bgee team.

#### Core reference intergenic from Bgee

`Core reference intergenic` releases  are created by the Bgee team when a lot of new RNA-Seq libraries have been manually curated for already existing species and/or for new species. These releases are the only ones with a release number (e.g "0.1"). Each of these releases contains reference intergenic sequences for a list of species.
Bgee reference intergenic sequences have been generated using Bgee team expertise. The RNA-Seq libraries were manually curated as healthy and wild type. Quality Control have been done along all steps of generation of these sequences. Reference intergenic sequences have been selected from all potential intergenic regions (see [Bgee pipeline](https://github.com/BgeeDB/bgee_pipeline/tree/master/pipeline/RNA_Seq)).
BgeeCall allows to generate gene expression call from Bgee reference intergenic sequences for any RNA-Seq libraries as long as these sequences have been generated by the Bgee team. 
A tsv file containing all species available for current release of reference intergenic is available [here](https://bgee.org/ftp/intergenic/current/species_info.tsv). This file also contains a column describing the number of RNA-Seq libraries used to generated the reference intergenic sequences of each species.
It is also possible to list in R all species for which Bgee reference intergenic sequences have been created :
```{r}
list_bgee_ref_intergenic_species(myBgeeMetadata = bgee)
```

#### Community reference intergenic

If you want to use BgeeCall on a species for which Bgee does not provide reference intergenic sequences you have the possibility to create them by yourself and share them with the Bgee community by following all steps of [this tutorial](https://github.com/BgeeDB/reference_intergenic_standalone#protocol-to-generate-reference-intergenic-sequences). Do not forget that the number of RNA-Seq libraries is a key point to the generation of precise reference intergenic sequences.
It is possible to list in R all species for which reference intergenic sequences have been created by the community using the following code
```{r}
list_community_ref_intergenic_species()
```
If reference intergenic sequences of the species you are interested in are available only from the community release it is then possible to use this release to generate your present/absent calls
```{r, eval=FALSE}
# create a BgeeMetadata object using the community release
bgee <- new("BgeeMetadata", release = "community")
calls_output <- generate_calls_workflow(bgeeMetadata = bgee, userMetadata = user_BgeeCall)
```

#### Your own reference intergenic 

If you generated your own reference intergenic sequences follwowing [this tuorial](https://github.com/BgeeDB/reference_intergenic_standalone#protocol-to-generate-reference-intergenic-sequences) but did not share them for the moment (do not forget to do it...), it is also possible to use BgeeCall with a file containing the sequences. In this case you need to select the custom release and provide the path to the file containing reference intergenic sequences :
```{r, eval=FALSE}
bgee <- new("BgeeMetadata", release = "community", custom_intergenic_path = "path/to/custom/ref_intergenic.fa.gz")
calls_output <- generate_calls_workflow(bgeeMetadata = bgee, userMetadata = user_BgeeCall)
```
### <a name="transcript_level"></a>Generate present/absent calls at transcript level (beta version)

kallisto generates TPMs at the transcript level. In the Bgee pipeline we summarize this expression at the gene level to calculate our present/absent calls.
In `BgeeCall` it is now possible to generate present/absent calls at the transcript level. Be careful when using this feature as it has not been tested for the moment.
To generate such calls you only have to create one object of the class `KallistoMetadata` and edit the value of one attribute
``` {r, eval=FALSE}
kallisto <- new("KallistoMetadata", txOut = TRUE)
calls_output <- run_from_object(myAbundanceMetadata = kallisto, myUserMetadata = user_BgeeCall)
```

### Tune how to use kallisto

#### Download or reuse your own kallisto
By default `BgeeCall` will download the version 0.45 of kallisto and will use it to quantify abundance of transcripts. It will only be used by this package and will have no impact on your potential already existing version of kallisto.
However, You can use a version already installed on your conputer or your cluster. It can be useful if you prefer to use an older version of kallisto.
To do that, you only have to create one object of the class `KallistoMetadata` and edit the value of one attribute
``` {r, eval=FALSE}
kallisto <- new("KallistoMetadata", download_kallisto = FALSE)
calls_output <- run_from_object(myAbundanceMetadata = kallisto, myUserMetadata = user_BgeeCall)
```
To download and install kallisto please follow the instructions here : http://pachterlab.github.io/kallisto/download

#### Edit kallisto quant attributes
By default kallisto is run with the same parameters that we use in the RNA-Seq Bgee pipeline:

* single end : "-t 1 --single -l 180 -s 30 --bias"
* paired end : "-t 1 --bias"

It is possible to modify them and use your favourite kallisto parameters
``` {r, message = FALSE, warning = FALSE}
kallisto <- new("KallistoMetadata", single_end_parameters = "-t 3 --single -l 150 -s 30", pair_end_parameters = "-t 2 -b --seed 36")
calls_output <- run_from_object(myAbundanceMetadata = kallisto, myUserMetadata = user_BgeeCall)
```

#### Choose between two kmer size
By default 2 indexes with 2 different kmer sizes can be used by `BgeeCall`
The default kmer size of kallisto (31) is used for libraries with reads length equal or larger than 50 bp.
A kmer size of 15 is used for libraries with reads length smaller than 50 bp.
We decided not to allow to tune kmers size because the generation of the index is time consuming and index generation takes even more time with small kmers size (< 15 bp). However it is possible to modify the threshold of read length allowing to choose between default and small kmer size.

``` {r, message = FALSE, warning = FALSE}
# libraries with reads smaller than 70bp will use the index with kmer size = 15
kallisto <- new("KallistoMetadata", read_size_kmer_threshold = 70)
calls_output <- run_from_object(myAbundanceMetadata = kallisto, myUserMetadata = user_BgeeCall)
```

### <a name="run_ids"></a>Generate calls for a subset of RNA-Seq runs
By default gene expression calls are generated using all runs of the RNA-Seq library. It is possible to select only a subset of these runs.
``` {r}
# RNA-Seq run SRR350955_subsetof from the RNA-Seq library will be used to generate the calls
user_BgeeCall@run_ids <- c("SRR350955_subset")
calls_output <- run_from_object(myUserMetadata = user_BgeeCall)
```
When run IDs are selected, the name output directory combine the library ID and all selected run IDs. In our example the expression calls will be stored in the directory `SRX099901_SRR1_SRR2`.

### Modify present/absent threshold

By default the threshold of present/absent is calculated with the formula :

    proportion of ref intergenic present / proportion of protein coding present = 0.05

This 0.05 corresponds to the ratio used in the Bgee pipeline. However it is possible to edit this value. Be careful when editing this value as it has a big impact on your present absent.

```{r}
kallisto <- new("KallistoMetadata", cutoff = 0.1)
```

### Generate calls for same library using different transcriptome or annotation versions
By default the arborescence of directories created by `BgeeCall` is as simple as possible. the results will be created using the path `working_path/intergenic_release/all_results/libraryId`. Generating present/absent gene expression calls for the same RNA-Seq library using different transcriptome or annotation versions using this arborescence will overwrite previous results. 
The `UserMetadata` class has an attribute `simple_arborescence` that is `TRUE` by default. If `FALSE`, a complexe arborescence of directories containing the name of the annotation and transcriptome files will be created. This complex arborescence will then allow to generate present/absent calls for the same library using different version of transcriptome or annotaiton.
``` {r, eval=FALSE}
user_BgeeCall@run_ids <- ""
user_BgeeCall@simple_arborescence <- FALSE
calls_output <- run_from_object(myUserMetadata = user_BgeeCall)
```

### Change directory where calls are saved
By default directories used to save present/absent calls are subdirectories of `UserMetadata@working_path`. However it is possible to select the directory where you want the calls to be generated.
```{r}
user_BgeeCall@output_dir <- "path/to/calls/for/this/library/"
```
This output directory will only contains results generated at the RNA-Seq library level. All data generated at species level are still stored using the `UserMetadata@working_path`. They can then still be reused to generate calls from other libraries of the same species. 

### <a name="slurm_options"></a>Modify slurm options

Two functions are available to run BgeeCall on a slurm queuing system. Parameters described below are available for both of them.

#### Number of jobs

The full idea of using a cluster is to parallelize your jobs. By default 10 jobs are run at the same time. It is possible to modify this number with the parameter `nodes`.

```{r, eval=FALSE}
# run 50 jobs in parallel
generate_slurm_indexes(userFile = "path/to/file.tsv", nodes = 50)
```

#### Do not submit the jobs

In order to be able to check files automatically generated to run the jobs it is possible to generate these files without submiting your jobs.
More information on created files are available on the vignette of the (https://cran.r-project.org/web/packages/rslurm/vignettes/rslurm.html)[rslurm package]. 

```{r, eval=FALSE}
# create temporary files but do not submit the jobs
generate_slurm_indexes(userFile = "path/to/file.tsv", submit = FALSE)
```

#### Modify slurm options

A bash scirpt is automatically created to run the jobs. This script contains default slurm options (array, cpus-per-task, job-name, output).
All other slurm options recognized by the sbatch command can be updated b creating a named list where name correspond to long name of options (e.g do not use 'p' but 'partition'). 

```{r, eval=FALSE}
# add slurm options to the sbatch script
slurm_options_index <- list(account = "account", time = "2:00:00", partition = "partition", mem = "30G")
generate_slurm_indexes(userFile = "path/to/file.tsv", slurm_options = slurm_options_index)
```

#### Add modules to your environment

In some cluster programs are not loaded by default. The modules parameter allows to load them by adding one line in the sbatch script.
This option has been implemented to add modules but could potentially be used to add any custom line of code in the sbatch script.

```{r, eval=FALSE}
# load R 3.6.1 and kallisto in a cluster environment where software has to loaded manually
modules <- c("module add R/3.6.1;", "module add kallisto;")
generate_slurm_indexes(userFile = "path/to/file.tsv", modules = modules)
```

#### Modify BgeeCall objects

By default except for columns present in the tsv file all other slots of the 3 BgeeCall classes will use default values.
In order to tune these parameters it is possible to create the objects and pass them to the slurm functions.
\strong{ Important :} When generating these objects it is mandatory to keep the same name as in the example below.

```{r, eval=FALSE}
# create BgeeCall objects and use them to generate indexes
kallistoMetadata <- new("KallistoMetadata", download_kallisto=TRUE)
userMetadata <- new("UserMetadata", working_path = "/path/to/working/dir")
bgeeMetadata <- new("BgeeMetadata", intergenic_release = "0.1")
generate_slurm_indexes(userFile = "path/to/file.tsv", kallistoMetadata = kallistoMetadata, bgeeMetadata = bgeeMetadata, userMetadata = userMetadata)
```

