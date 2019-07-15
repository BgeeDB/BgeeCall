`BgeeCall` is a collection of functions that uses [Bgee](https://bgee.org/) expertise to create gene expression present/absent calls

The `BgeeCall` package allows to: 

* Generate calls of presence/absence of expression at the gene level. You can generate these calls for any RNA-Seq samples as long as the species is available in Bgee (information available from the `list_bgee_species()` function).
* Download reference intergenic sequences for species available in Bgee.
* Generate calls of presence/absence of expression at the transcript level (beta version)

If you find a bug or have any issues with `BgeeCall` please write a bug report in our GitHub issues manager available at (URL).

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
BiocManager::install("BgeeCall", version = "3.9")
```

## How to use the BgeeCall package

BgeeCall is highly tunable. Do not hesitate to have a look at the reference manual to have a precise descripton of all slots of the 4 main S4 classes (AbundanceMetadata, KallistoMetadata, BgeeMetadata and UserMetadata) or of all available functions.


### Load the package
``` {r, message = FALSE, warning = FALSE}
library(BgeeCall)
```


### Quick start

With the BgeeCall package it is easy to generate present/absent gene expression calls.
The most time comsuming task of this calls generation is the generation of the kallisto transcriptome index.
As the time needed for this step depend on the size of the transcriptome, we choose, as an example, the smallest transcriptome file
among all species available on Bgee (C. elegans).
To generate these calls you will need :

- a transcriptome
- gene annotations
- your RNA-Seq reads in fastq files

For this vignette we created a toy fastq file example based on the SRX099901 library using the [ShortRead](http://bioconductor.org/packages/release/bioc/html/ShortRead.html) R package

``` {r, eval=FALSE}
library("ShortRead")
# keep 48.000 reads
sampler <- FastqSampler(file.path("absolute_path","/SRX099901/SRR350955.fastq.gz"), 48000)
set.seed(1); SRR350955 <- yield(sampler)
writeFastq(object = SRR350955, file =file.path( "absolute_path","SRX099901_subset", "SRR350955_subset.fastq.gz"), mode = "w", full = FALSE, compress = TRUE)
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
# provide path to the directory of your RNA-Seq library
user_BgeeCall <- setRNASeqLibPath(user_BgeeCall, 
                                  system.file("extdata", "SRX099901_subset", package = "BgeeCall"))
```

And that's it... You can run the generation of your present/absent gene expression calls
``` {r, message = FALSE, warning = FALSE}
calls_output <- run_from_object(myUserMetadata = user_BgeeCall)
```
Each analyze generates 4 files and return path to each one of them.

* calls_tsv_path : path to main tsv file with TPM, count, length, biotype, type, and presence/absence of expression summarized at gene level (or at the [transcript level](#transcript_level) if it was requested)
``` {r, message = FALSE, warning = FALSE}
head.DataTable(x = read.table(calls_output$calls_tsv_path, header = TRUE), n = 5)
```
* cutoff_info_file_path : path to tsv summary of the analyze containing the proportion of gene, protein coding and intergenic defined as expressed. It also contains the library ID and the present/absent TPM threshold
``` {r, message = FALSE, warning = FALSE}
read.table(calls_output$cutoff_info_file_path)
```
* abundance_tsv : path to tsv kallisto quant output file
``` {r, message = FALSE, warning = FALSE}
head.DataTable(x = read.table(calls_output$abundance_tsv, header = TRUE), n = 5)
calls_output$TPM_distribution_path
calls_output$abundance_tsv
```
* TPM_distribution_path : path to plot in pdf reprensting density distribution of TPM values for all sequences, protein coding sequences, and intergenic sequences. The grey line corresponds to TPM threshold used to generate present/absent calls.
``` {r, eval = FALSE}
openPDF(calls_output$TPM_distribution_path)
```

###  Generate present/absent calls for more than one RNA-Seq library

The function `run_from_object()` is perfect to generate calls for one library. You will potentialy be also interested to run more than one call generation at the same time. It is possible to do that by using the ` run_from_file()` or the `run_from_dataframe()` functions.
With these functions you will be able to run calls generation for different:

- RNA-Seq libraries
- transcriptome
- gene annotation
- runs from the same RNA-Seq library
- species as long as they are part of Bgee

A template of the file usable as input of the function `run_from_file()` is available at the root directory of the package with the name `userMetadataTemplate.tsv`.
In this template each column correspond to one parameter used to generate gene expression calls. Each line will correspond to one expression calls generation analyze.
It is not mandatory to add a value to the `run_ids` column except if you want to generate expression calls for a subset of the runs of one RNA-Seq library as described in [Generate calls for a subset of RNA-Seq runs](#run_ids)
Once the file has been fill in expression calls can be generated with :
``` {r, eval=FALSE}
run_from_file(userMetadataFile = "path_to_your_file.tsv")
```

### list species available on Bgee

BgeeCall allows to generate gene expression call for any RNA-Seq libraries as long as the species is present in Bgee. To see all species in the last version of Bgee run :

```{r}
list_bgee_species()
```

### list reference intergenic releases

Different releases of Bgee reference intergenic sequences are available. It is possible to list all these releases :
```{r}
list_intergenic_release()
```
It is then possible to choose one specific release to create a `BgeeMetadata` object.
```{r}
bgee <- new("BgeeMetadata", intergenic_release = "0.1")
```
By default the intergenic used when a `BgeeMetadata`object is created is the last created one.

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