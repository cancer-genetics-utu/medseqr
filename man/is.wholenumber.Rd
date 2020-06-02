\name{is.wholenumber}
\alias{is.wholenumber}
\title{
  Test if a given number is a 'whole number'
}
\description{
  \command{is.wholenumber} uses information on the numerical
  characteristics of the machine R is running on to test if the input is a whole number
}
\usage{
  is.wholenumber(x, tol = .Machine$double.eps^0.5)
}
\arguments{
  \item{x}{
    numeric, a number to be tested
  }
  \item{tol}{
    double, tolerance to be tested against
  }
}
\value{
  A logical, TRUE or FALSE
}
\seealso{
  \command{\link[base]{.Machine}}
}
\keyword{ utilities }
