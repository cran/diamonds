\name{diamond.edges}
\alias{diamond.edges}
\title{Coordinates of edges of diamonds}
\description{
  Creates \code{x} and \code{y} coordinates of edges of
  a diamond partition, suitable for \code{\link{lines}}.
}
\usage{
  diamond.edges (b, d)
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
  lines (diamond.edges (base, 4))
  lines (diamond.edges (base, 1), lwd=2)
  title ("Diamond Partition at Level 4", line=0)
}
\keyword{dplot}
