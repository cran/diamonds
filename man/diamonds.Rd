\name{diamonds}
\alias{diamonds}
\alias{subdiamonds}
\alias{hexlines}
\alias{hiercenters}
\alias{subcenters}
\alias{hierlabels}
\alias{deep}
\alias{labeldeep}
\alias{centerdeep}
\alias{randlevel}
\alias{ranptdia}
\title{Internal functions for diamonds package}
\description{
  Various functions for internal use in the diamonds package.
}
\usage{
  diamonds(b, d)
  subdiamonds (b)
  hexlines (b, d)
  hiercenters (b, d)
  subcenters (b) 
  hierlabels (parent, level, d)
  deep (b, level, d)
  labeldeep (dt, parent, level, d)
  centerdeep (b, level, d)
  randlevel (tl, level, d)
  ranptdia (b)
}
\arguments{
  \item{b}{list of coordinates of vertices of a diamond}
  \item{d}{final depth of recursion of partition}
  \item{level}{current depth of recursion}
  \item{parent}{label of next higher level diamond}
  \item{dt}{list of labels being build recursively}
  \item{tl}{list of diamonds}
}
\details{See R codes.}
\author{Denis White, \email{white.denis@epa.gov} }
\seealso{\code{\link{diamond.base}}}
\keyword{dplot}
