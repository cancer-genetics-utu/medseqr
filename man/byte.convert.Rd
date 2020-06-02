\name{byte.convert}
\alias{byte.convert}
\title{
  Convert file sizes between different units
}
\description{
  \command{byte.convert} converts the size of a given file (or any file-size-like character) to another unit.
  By default the input is converted to Bytes.
}
\usage{
  byte.convert(x, to=c("b", "k", "m", "g", "t"), digits=2)
}
\arguments{
  \item{x}{
    numeric or character, a file-size-like vector to be converted
  }
  \item{to}{
    character or NULL, unit to convert to. If \code{NULL} then the next unit relative to the unit of \option{x} in the order of \option{to} is used,
    e.g., if \option{x} is "256K" it will be converted to Megabytes.
  }
  \item{digits}{
    length-one numeric, number of digits to round the result to, passed to \command{round}
  }
}
\value{
  A numeric or character vector, depending on the target unit
}
\keyword{ utilities }
