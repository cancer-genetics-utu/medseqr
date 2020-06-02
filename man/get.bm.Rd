\name{get.bm}
\alias{get.bm}
\title{
  Retrieve additional annotations from biomart
}
\description{
  \command{get.bm} is a wrapper for \command{getBM} from the \emph{biomaRt} package.
  It sets all needed variables and performs the query.
}
\usage{
  get.bm(values, biom.data.set="hsapiens_gene_ensembl", biom.mart=c("ensembl", "snp", "funcgen", "vega", "pride"),
         host="www.ensembl.org", biom.filter="ensembl_gene_id", biom.attributes=c("ensembl_gene_id","hgnc_symbol","description"))
}
\arguments{
  \item{values}{
    character, ids to be converted, passed to the'values' argument of \command{getBM}
  }
  \item{biom.data.set}{
    character of length one, biomart data set to use
  }
  \item{biom.mart}{
    character vector, biomart to use (uses the first element of the vector), defaults to "ensembl"
  }
  \item{host}{
    character of length one, host URL
  }
  \item{biom.filter}{
    character of length one, name of biomart filter, i.e., type of query ids, defaults to "ensembl_gene_id"
  }
  \item{biom.attributes}{
    character vector, biomart atrributes, i.e., type of desired result(s); make sure query id type is included!
  }
}
\value{
  A character vector with the converted ids.
}
\seealso{
  \command{\link[biomaRt]{getBM}}
}
\keyword{ utilities }
