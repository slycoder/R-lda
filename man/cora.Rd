\name{cora}
\Rdversion{1.1}
\alias{cora}
\alias{cora.documents}
\alias{cora.vocab}
\alias{cora.cites}
\alias{cora.titles}
\docType{data}
\title{
  A subset of the Cora dataset of scientific documents.
}
\description{
  A collection of 2410 scientific documents in LDA format with links and titles
  from the Cora search engine.  
}
\usage{
data(cora.documents)
data(cora.vocab)
data(cora.cites)
data(cora.titles)
}

\format{
  \code{cora.documents} and \code{cora.vocab}
  comprise a corpus of 2410 documents conforming to the LDA format.

  \code{cora.titles} is a character vector of titles for each
  document (i.e., each entry of \code{cora.documents}).

  \code{cora.cites} is a list representing the citations between the
  documents in the collection (see related for format).
}
\source{
  \cite{Automating the construction of internet protals with machine
  learning.  McCallum et al.  Information Retrieval.  2000.}
}
\seealso{
  \code{\link{lda.collapsed.gibbs.sampler}} for the format of the
  corpus.

  \code{\link{rtm.collapsed.gibbs.sampler}} for the format of the
  citation links.
}
\examples{
data(cora.documents)
data(cora.vocab)
data(cora.links)
data(cora.titles)
}
\keyword{datasets}
