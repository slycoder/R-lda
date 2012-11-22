\name{poliblog}
\Rdversion{1.1}
\alias{poliblog}
\alias{poliblog.documents}
\alias{poliblog.vocab}
\alias{poliblog.ratings}
\docType{data}
\title{
  A collection of political blogs with ratings.
}
\description{
  A collection of 773 political blogs in LDA format with
  conservative/liberal ratings.
}
\usage{
data(poliblog.documents)
data(poliblog.vocab)
data(poliblog.ratings)
}
\format{
  \code{poliblog.documents} and \code{poliblog.vocab}
  comprise a corpus of 773 political blogs conforming to the LDA format.

  \code{poliblog.ratings} is a numeric vector of length 773 which gives
  a rating of liberal (-100) or conservative (100) to each document in
  the corpus.
}
\source{
  \cite{Blei, David M. and McAuliffe, John.  Supervised topic models.
    Advances in Neural Information Processing Systems, 2008.}
}
\seealso{
  \code{\link{lda.collapsed.gibbs.sampler}} for the format of the
  corpus.
}
\examples{
data(poliblog.documents)
data(poliblog.vocab)
data(poliblog.ratings)
}
\keyword{datasets}
