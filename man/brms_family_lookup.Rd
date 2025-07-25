% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lookups.R
\name{brms_family_lookup}
\alias{brms_family_lookup}
\title{Lookup function for brms families via string identifier}
\usage{
brms_family_lookup(family, link = NULL)
}
\arguments{
\item{family}{String identifier of the family.}

\item{link}{Link to be passed to the family function.}
}
\value{
A brmsfamily object matching the string identifier and using the link
}
\description{
Lookup function for brms families via string identifier
}
\examples{
brms_family_lookup("weibull", "softplus")
}
