\name{diamond.dualcents}
\alias{diamond.dualcents}
\title{Coordinates of centers of hexagons}
\description{
  Creates \code{x} and \code{y} coordinates of center points of
  the hexagon duals of a diamond partition, suitable for
  \code{\link{points}}.
}
\usage{
  diamond.dualcents (b, d)
}
\arguments{
  \item{b}{list of coordinates of vertices of a diamond, as from 
    \code{\link{diamond.base}}}
  \item{d}{depth of recursion of partition; if \code{d=1},
    then no partition}
}
\details{
  Centers are generated for "left-side" hexagons, assuming that
  diamonds are joined as in octahedron or icosahedron.  See reference.
}
\value{
  List with \code{x} and \code{y} components.
}
\references{
  White, D., (2000) Global grids from recursive diamond subdivisions
  of the surface of an octahedron or icosahedron, \emph{Environmental
  Monitoring and Assessment}, \bold{64}(1), 93-103.
}
\author{Denis White, \email{white.denis@epa.gov} }
\examples{
  base <- diamond.base ()
  diamond.plot (diamond.edges (base, 1))
  lines (diamond.dualedges (base, 4), col="blue")
  points (diamond.dualcents (base, 4), pch=19, col="red")
  lines (diamond.edges (base, 1), lwd=2, col="black")
  title ("Centers of Dual Hexagons at Level 4", line=0)
}
\keyword{dplot}
