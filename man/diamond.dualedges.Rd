\name{diamond.dualedges}
\alias{diamond.dualedges}
\title{Coordinates of edges of hexagons}
\description{
  Creates \code{x} and \code{y} coordinates of edges of
  the hexagon duals of a diamond partition, suitable for
  \code{\link{lines}}.
}
\usage{
  diamond.dualedges (b, d)
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
  lines (diamond.dualedges (base, 4), col="blue")
  lines (diamond.edges (base, 1), lwd=2, col="black")
  title ("Dual Hexagons at Level 4", line=0)
}
\keyword{dplot}
