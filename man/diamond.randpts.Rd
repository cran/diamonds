\name{diamond.randpts}
\alias{diamond.randpts}
\title{Coordinates of random points}
\description{
  Creates \code{x} and \code{y} coordinates of random points in
  each diamond element of a diamond partition, suitable
  for \code{\link{points}}.
}
\usage{
  diamond.randpts (b, d)
}
\arguments{
  \item{b}{list of coordinates of vertices of a diamond, as from 
    \code{\link{diamond.base}}}
  \item{d}{depth of recursion of partition; if \code{d=1},
    then no partition}
}
\value{
  List with \code{x} and \code{y} components.
}
\author{Denis White, \email{white.denis@epa.gov} }
\examples{
  base <- diamond.base ()
  diamond.plot (diamond.edges (base, 1))
  lines (diamond.edges (base, 4), lwd=1, col="blue")
  lines (diamond.edges (base, 1), lwd=2, col="black")
  points (diamond.randpts (base, 4), pch=19, col="red")
  title ("Random Points at Level 4", line=0)
}
\keyword{dplot}
