% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transcriptToGene.R
\name{run_tximport}
\alias{run_tximport}
\title{Run tximport}
\usage{
run_tximport(
  myAbundanceMetadata = new("KallistoMetadata"),
  myBgeeMetadata = new("BgeeMetadata"),
  myUserMetadata,
  abundanceFile = ""
)
}
\arguments{
\item{myAbundanceMetadata}{A descendant object of the Class 
myAbundanceMetadata.}

\item{myBgeeMetadata}{A Reference Class BgeeMetadata object.}

\item{myUserMetadata}{A Reference Class UserMetadata object.}

\item{abundanceFile}{(Optional) Path to the abundance file. NULL by default.
If not NULL, the file located at `abundanceFile` will be used to run tximport.
Otherwise (Default) the path to the abundance file is deduced fom attributes of
classes `BgeeMetadata`, `UserMetadata` and `AbundanceMetadata`}
}
\value{
a tximport object
}
\description{
Run tximport. Will summarize abundance estimation from transcript 
level to gene level if `myAbundanceMetadata@txout == FALSE`. 
Otherwise keep abundance estimation at transcript level.
}
\examples{
{
user <- new("UserMetadata", working_path = system.file("extdata", 
    package = "BgeeCall"), species_id = "6239", 
  rnaseq_lib_path = system.file("extdata", 
    "SRX099901_subset", package = "BgeeCall"), 
  annotation_name = "WBcel235_84", simple_arborescence = TRUE)
abundance_file <- system.file('extdata', 'abundance.tsv', package = 'BgeeCall')
tx_import <- run_tximport(myUserMetadata = user, 
abundanceFile = abundance_file)
}

}
\author{
Julien Wollbrett
}
