\name{diamond.randlabels}
\alias{diamond.randlabels}
\title{Labels for a hierarchically randomized partition}
\description{
  Creates a vector of character labels for diamonds in a
  partition, where the labels are randomly permuted.
}
\usage{
  diamond.randlabels (d)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{d}{depth of recursion of partition; if \code{d=1},
    then no partition}
}
\details{
  The labels are generated from a quadrant-recursive ordering,
  then randomly permuted at each level in the recursive hierarchy,
  see reference.
}
\value{
  A vector of the hierarchically randomized character labels for
  each diamond in the partition.
}
\references{
  White, D., (2000) Global grids from recursive diamond subdivisions
  of the surface of an octahedron or icosahedron, \emph{Environmental
  Monitoring and Assessment}, \bold{64}(1), 93-103.
}
\author{Denis White, \email{white.denis@epa.gov} }
\seealso{
  \code{\link{diamond.randwalk}},
  \code{\link{diamond.labels}}
}
\examples{
  diamond.randlabels (2)
}
\keyword{dplot}
