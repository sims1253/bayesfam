% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lookups.R
\name{aux_family_parameters_lookup}
\alias{aux_family_parameters_lookup}
\title{Lookup function for the names of the auxiliary parameters of a likelihood}
\usage{
aux_family_parameters_lookup(family)
}
\arguments{
\item{family}{The identifier string of a family.}
}
\value{
Character vector of auxiliary parameter names.
}
\description{
Lookup function for the names of the auxiliary parameters of a likelihood
}
\examples{
aux_family_parameters_lookup("beta")
}
