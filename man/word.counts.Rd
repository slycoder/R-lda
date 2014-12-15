\name{word.counts}
\Rdversion{1.1}
\alias{word.counts}
\alias{document.lengths}
\title{
  Compute Summary Statistics of a Corpus
}
\description{
  These functions compute summary statistics of a corpus.
  \code{word.counts} computes the word counts for a set of documents,
  while \code{documents.length} computes the length of the documents in
  a corpus.
}
\usage{
word.counts(docs, vocab = NULL)

document.lengths(docs)
}
\arguments{
  \item{docs}{
    A list of matrices specifying the corpus.  See
    \code{\link{lda.collapsed.gibbs.sampler}} for details on the
    format of this variable.
  }
  \item{vocab}{
    An optional character vector specifying the levels (i.e., labels) of
    the vocabulary words.  If unspecified (or \code{NULL}), the levels
    will be automatically inferred from the corpus.
  }
}
\value{
  \code{word.counts} returns an object of class \samp{table} which
  contains counts for the number of times each word appears in the input
  corpus.  If \var{vocab} is specified, then the levels of the table
  will be set to \var{vocab}.  Otherwise, the levels are automatically
  inferred from the corpus (typically integers \var{0:(V-1)}, where
  \var{V} indicates the number of unique words in the corpus).

  \code{documents.length} returns a integer vector of length
  \code{length(docs)}, each entry of which corresponds to the
  \dfn{length} (sum of the counts of all features) of each document in
  the corpus.
}
\author{
  Jonathan Chang (\email{slycoder@gmail.com})
}

\seealso{
  \code{\link{lda.collapsed.gibbs.sampler}} for the input format of
  these functions.

  \code{\link{read.documents}} and \code{\link{lexicalize}} for ways of
  generating the input to these functions.

  \code{\link{concatenate.documents}} for operations on a corpus.
}
\examples{
## Load the cora dataset.
data(cora.vocab)
data(cora.documents)

## Compute word counts using raw feature indices.
wc <- word.counts(cora.documents)
head(wc)
##   0   1   2   3   4   5 
## 136 876  14 111  19  29 

## Recompute them using the levels defined by the vocab file.
wc <- word.counts(cora.documents, cora.vocab)
head(wc)
##   computer  algorithms discovering    patterns      groups     protein 
##        136         876          14         111          19          29 

head(document.lengths(cora.documents))
## [1] 64 39 76 84 52 24
}
\keyword{ utilities }
