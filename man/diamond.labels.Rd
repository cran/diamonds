\name{diamond.labels}
\alias{diamond.labels}
\title{Labels for a diamond partition}
\description{
  Creates a vector of character labels for diamonds in a
  partition.
}
\usage{
  diamond.labels (d)
}
\arguments{
  \item{d}{depth of recursion of partition; if \code{d=1},
    then no partition}
}
\details{
  The labels are generated from a quadrant-recursive ordering,
  see reference.
}
\value{
  A vector of character labels for each diamond in the partition.
}
\references{
  White, D., (2000) Global grids from recursive diamond subdivisions
  of the surface of an octahedron or icosahedron, \emph{Environmental
  Monitoring and Assessment}, \bold{64}(1), 93-103.
}
\author{Denis White, \email{white.denis@epa.gov} }
\seealso{
  \code{\link{diamond.hierwalk}},
  \code{\link{diamond.randlabels}}
}
\examples{
  base <- diamond.base ()
  diamond.plot (diamond.edges (base, 1))
  lines (diamond.edges (base, 4), lwd=1, col="blue")
  lines (diamond.edges (base, 1), lwd=2, col="black")
  text (diamond.centers (base, 4), diamond.labels (4))
  title ("Morton/Peano Numbers at Level 4", line=0)
}
\keyword{dplot}
