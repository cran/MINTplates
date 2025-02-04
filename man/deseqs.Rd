% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MINTplates.R
\name{deseqs}
\alias{deseqs}
\title{Decode license-plates}
\usage{
deseqs(plates)
}
\arguments{
\item{plates}{The license plates being decoded.}
}
\value{
The sequences they decodes to.
}
\description{
Decode the license-plates using the lookup table.
}
\examples{
plates=c('tRF-18-BS6PDFD2','tRF-20-51K36D26')
deseqs(plates)
}
\references{
Pliatsika V, Loher P, Telonis AG, Rigoutsos I. Bioinformatics (2016) <doi: 10.1093/bioinformatics/btw194>
}
