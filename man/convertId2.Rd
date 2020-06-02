\name{convertId2}
\alias{convertId2}
\title{
  Convert gene IDs between various types
}
\description{
  \command{convertId2} converts IDs of Human or Mouse genes between
  \bold{EG} (\emph{Entrez Gene}), \bold{ENSEMBL} (\emph{Ensembl ID}) and \bold{SYMBOL}
  (\emph{HGNC gene symbol}).
}
\usage{
  convertId2(id, species = "Human")
}
\arguments{
  \item{id}{
    character value or vector, IDs to be converted
  }
  \item{species}{
    character, Species name, determines the annotation package to be used
  }
}
\value{
  A named character vector.
}
\note{
  The function was initially written by Henrik Edgren and later adapted by Vidal Fey.
}
\keyword{ utilities }
