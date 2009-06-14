\name{diamond.base}
\alias{diamond.base}
\title{Create a base diamond}
\description{
  Generates vertices of a portrait format diamond with center at
  \code{(x, y)} and short axis \code{d}.
}
\usage{
  diamond.base (x=0, y=0, a=1)
}
\arguments{
  \item{x}{x coordinate of diamond center}
  \item{y}{y coordinate of diamond center}
  \item{a}{length of short axis of diamond}
}
\details{
  This function is normally called first to establish a base diamond
  from which recursive partitions are derived.  The partitions are
  created from 4-fold subdivision into \code{4^(d-1)} congruent
  sub-diamonds.

  A "diamond" in this package means the figure formed by two equilateral
  triangles connected across a common edge, or, two adjacent faces of
  an octahedron or icosahedron.

  Recursion in this package limits the useful size of problem.  A version
  with C routines is in development.
}
\value{
  A list with four components, each a pair of \code{(x, y)} coordinates
  describing a vertex of the diamond.
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
  lines (diamond.edges (base, 4))

  # kites for levels 2-4
  diamond.plot (diamond.edges (base, 1))
  lines (diamond.edges (base, 4), lwd=1, col="blue")
  lines (diamond.triedges (base, 4), lwd=1, col="blue")
  lines (diamond.dualedges (base, 4), lwd=1, col="blue")
  lines (diamond.edges (base, 3), lwd=2, col="green")
  lines (diamond.triedges (base, 3), lwd=2, col="green")
  lines (diamond.dualedges (base, 3), lwd=2, col="green")
  lines (diamond.edges (base, 2), lwd=3, col="red")
  lines (diamond.triedges (base, 2), lwd=3, col="red")
  lines (diamond.dualedges (base, 2), lwd=3, col="red")
  lines (diamond.edges (base, 1), lwd=2, col="black")
  title ("'Kites' for Levels 2-4", line=0)
}
\keyword{dplot}
