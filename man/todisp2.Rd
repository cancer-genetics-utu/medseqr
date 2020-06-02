\name{todisp2}
\alias{todisp2}
\title{
  Convert Human ensembl IDs to HGNC symbols.
}
\description{
  \command{todisp2} converts Human \emph{ensembl gene id}s to (\emph{HGNC gene symbol}s).
}
\usage{
  todisp2(ensg, lab, biomart=TRUE)
}
\arguments{
  \item{ensg}{
    character value or vector, IDs to be converted
  }
  \item{lab}{
    data.frame with one column holding the IDs that \option{ensg} should be converted to and with row names corresponding to \option{ensg} for merging
  }
  \item{biomart}{
    logical, should biomaRt be used to convert \option{ensg}? If \emph{FALSE} will fall back to \command{convertId2}.
  }
}
\value{
  A character vector of the same length as \option{ensg}, with input ids replaced by converted ids.
  If no replacement was found the input id is returned.
}
\keyword{ utilities }
