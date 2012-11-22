\name{sampson}
\Rdversion{1.1}
\alias{sampson}
\docType{data}
\title{
  Sampson monk data
}
\description{
  Various relationships between several monks at a monastery collected
  over time.
}
\usage{data(sampson)}

\format{
  \code{sampson} is a list whose entries are 18x18 matrices representing
  the pairwise relationships between 18 monks.  The names of the monks
  are given as the row/column names of each matrix.
  
  Each matrix encodes a different relationship (there are a total of 10)
  described by the
  corresponding name field of the list.
}
\source{
  \cite{F. S. Sampson.  A novitiate in a period of change: An experimental and
  case study of social relationships.  PhD thesis, Cornell University. 1968.}
}
\seealso{
  \code{\link{mmsb.collapsed.gibbs.sampler}} is an example of a
  function which can model the structure of this data set.
}
\examples{
data(sampson)
}
\keyword{datasets}
