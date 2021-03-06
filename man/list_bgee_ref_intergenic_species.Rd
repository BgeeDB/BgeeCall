% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/retrieveFromBgee.R
\name{list_bgee_ref_intergenic_species}
\alias{list_bgee_ref_intergenic_species}
\title{List species having Bgee reference intergenic sequences}
\usage{
list_bgee_ref_intergenic_species(myBgeeMetadata = NULL, release = NULL)
}
\arguments{
\item{myBgeeMetadata}{A Reference Class BgeeMetadata object}

\item{release}{A Bgee reference intergenic release name}
}
\value{
list all species having reference intergenic sequences available in the selected
release
}
\description{
Return information related to species having Bgee reference intergenic 
sequences available for the selected Bgee intergenic release: 
\itemize{
 \item{speciesId}{ the NCBI species ID of the species}
 \item{specieName}{ scientific species name}
 \item{numberOfLibraries}{ number of libraries used to generate these reference intergenic
sequences}
 \item{genomeVersion}{ version of the genome used to generate the reference intergenic sequences}
}
If a BgeeMetadata object is provided this function retrieve the list of species using 
BgeeMetadata@intergenic_release.
If only a `release` is provided it will use it to retrieve the list of species.
If none of them are provided the default Bgee reference intergenic release will be used.
}
\examples{
{
bgee <- new("BgeeMetadata")
list_bgee_ref_intergenic_species(myBgeeMetadata = bgee)
list_bgee_ref_intergenic_species(release = '0.2')
}

}
\author{
Julien Wollbrett
}
