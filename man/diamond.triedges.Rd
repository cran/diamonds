\name{diamond.triedges}
\alias{diamond.triedges}
\title{Coordinates of edges of triangles}
\description{
  Creates \code{x} and \code{y} coordinates of bisecting edges of
  diamonds in partition, providing for two equilateral triangles,
  with output suitable for \code{\link{lines}}.
}
\usage{
  diamond.triedges (b, d)
}
\arguments{
  \item{b}{list of coordinates of vertices of a diamond, as from 
    \code{\link{diamond.base}}}
  \item{d}{depth of recursion of partition; if \code{d=1},
    then no partition}
}
\details{
  Both this function and \code{\link{diamond.edges}} are necessary
  to show triangle partitions.  See example.
}
\value{
  List with \code{x} and \code{y} components.
}
\author{Denis White, \email{white.denis@epa.gov} }
\examples{
  base <- diamond.base ()
  diamond.plot (diamond.edges (base, 1))
  lines (diamond.edges (base, 4))
  lines (diamond.triedges (base, 4))
  title ("Triangles at Level 4", line=0)
}
\keyword{dplot}
