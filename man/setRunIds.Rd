% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UserMetadata.R
\docType{methods}
\name{setRunIds}
\alias{setRunIds}
\alias{setRunIds,UserMetadata,character-method}
\alias{setRunIds,userMetadata,character}
\title{`run_ids` Setter}
\usage{
setRunIds(userObject, runIds)

\S4method{setRunIds}{UserMetadata,character}(userObject, runIds)
}
\arguments{
\item{userObject}{The UserMetadata object}

\item{runIds}{id of all wanted runs}
}
\value{
An object of the class UserMetadata
}
\description{
Method of the class UserMetadata. Set run_ids of
one UserMetadata object by providing the id of all wanted runs
}
\examples{
{
user <- new("UserMetadata")
user <- setRunIds(user, c("RUN_1", "RUN_2"))
}

}
