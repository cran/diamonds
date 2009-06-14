\name{diamond.centers}
\alias{diamond.centers}
\title{Coordinates of centers of diamonds}
\description{
  Creates \code{x} and \code{y} coordinates of center points of
  a diamond partition, suitable for \code{\link{points}}.
}
\usage{
  diamond.centers (b, d)
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
  points (diamond.centers (base, 4), pch=19, col="red")
  lines (diamond.edges (base, 1), lwd=2, col="black")
  title ("Diamond Centers at Level 4", line=0)
}
\keyword{dplot}
