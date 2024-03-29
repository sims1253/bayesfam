% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test-helper.R
\name{isSingleString}
\alias{isSingleString}
\title{Check, if the input is a single string}
\usage{
isSingleString(input)
}
\arguments{
\item{input}{String argument}
}
\value{
Is a string and only one string
}
\description{
Check, if the input is a single string
}
\examples{
bayesfam:::isSingleString("abc") # should be TRUE
bayesfam:::isSingleString(c("abc", "def")) # should be FALSE, not a single string
}
