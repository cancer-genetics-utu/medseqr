\name{convert.bm}
\alias{convert.bm}
\title{
  Retrieve additional annotations from biomart
}
\description{
  \command{convert.bm} is a wrapper for \command{get.bm} which in turn makes use of \command{getBM} from the \emph{biomaRt} package.
  It takes a matrix or data frame with the IDs to be converted in one column or as row names as input and returns an annotated data frame.
}
\usage{
  convert.bm(dat, id="ID", biom.data.set="hsapiens_gene_ensembl", biom.mart=c("ensembl", "snp", "funcgen", "vega", "pride"),
             host="www.ensembl.org", biom.filter="ensembl_gene_id", biom.attributes=c("ensembl_gene_id","hgnc_symbol","description"),
             rm.dups=FALSE)
}
\arguments{
  \item{dat}{
    matrix or data frame with the ids to be converted in a column or as row names
  }
  \item{id}{
    character, name of the column with the ids to be converted, special name "rownames" will use the row names
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
  \item{rm.dups}{
    logical, should duplicated input IDs (\option{biom.filter}) be removed from the result?
  }
}
\value{
  A data frame with the retrieved information.
}
\seealso{
  \command{\link[biomaRt]{getBM}}
}
\keyword{ utilities }
