\name{diamond.plot}
\alias{diamond.plot}
\title{Define a fixed aspect ratio plot}
\description{
  Calls \code{\link{plot}} to create a plot space whose dimensions
  are defined by the parameters, and whose aspect ratio is fixed
  with respect to them.
}
\usage{
  diamond.plot (cl)
}
\arguments{
  \item{cl}{a list of coordinates with components \code{x} and
    \code{y}, such as returned by \code{\link{diamond.edges}}}
}
\details{
  This function calls \code{\link{plot}} with \code{type="n"}, and it
  maintains the aspect ratio of the input coordinates in the output
  graphic by setting \code{par(pin=)}, thus the map does not scale
  dynamically in an on-screen window.
}
\author{Denis White, \email{white.denis@epa.gov} }
\examples{
  base <- diamond.base ()
  diamond.plot (diamond.edges (base, 1))
  lines (diamond.edges (base, 4))
  lines (diamond.edges (base, 1), lwd=2)
  title ("Diamond Partition at Level 4", line=0) }
\keyword{hplot}
