% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UserMetadata.R
\docType{methods}
\name{setAnnotationFromFile}
\alias{setAnnotationFromFile}
\alias{setAnnotationFromFile,UserMetadata,character,missing-method}
\alias{setAnnotationFromFile,userMetadata,character,character}
\alias{setAnnotationFromFile,UserMetadata,character,character-method}
\title{Set annotation_object of one UserMetadata object}
\usage{
setAnnotationFromFile(userObject, annotationPath, annotationName)

\S4method{setAnnotationFromFile}{UserMetadata,character,missing}(userObject, annotationPath, annotationName)

\S4method{setAnnotationFromFile}{UserMetadata,character,character}(userObject, annotationPath, annotationName)
}
\arguments{
\item{userObject}{The UserMetadata object}

\item{annotationPath}{Absolute path to the annotation file}

\item{annotationName}{(optional) Name of the annotation. Will be used to
create folders.}
}
\value{
An object of the class UserMetadata
}
\description{
Method of the class UserMetadata. Set annotation_object of
one UserMetadata object  by providing the path to a fasta transcriptome file.
}
\details{
If no annotationName is provided the name of the annotation file
will be used to create folders.
}
\examples{
{
# path to gtf annotation file
annotation_file <- system.file("extdata", "annotation.gtf", package = "BgeeCall")
user <- new("UserMetadata")
user <- setAnnotationFromFile(user, annotation_file,
                             "annotation_name")
}

}
