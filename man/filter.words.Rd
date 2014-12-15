\name{filter.words}
\Rdversion{1.1}
\alias{filter.words}
\alias{concatenate.documents}
\alias{shift.word.indices}

\title{
  Functions to manipulate text corpora in LDA format.
}
\description{
  \code{concatenate.documents} concatenates a set of documents.
  \code{filter.words} removes references to certain words
  from a collection of documents.
  \code{shift.word.indices} adjusts references to words by a fixed amount.
}
\usage{
concatenate.documents(...)
filter.words(documents, to.remove)
shift.word.indices(documents, amount)
}
\arguments{
  \item{\dots}{
    For \code{concatenate.documents}, the set of corpora to be merged.  All
  arguments to \code{\ldots} must be corpora of the same length.  The
  documents in the same position in each of the arguments will be
  concatenated, i.e., the new document 1 will be the concatenation of
  document 1 from argument 1, document 2 from argument 1, etc.
}
  \item{documents}{
    For \code{filter.words} and \code{shift.word.indices}, the corpus to
  be operated on.
 }
  \item{to.remove}{
    For \code{filter.words}, an integer vector of words to filter.
    The words in each document which also exist in \code{to.remove} will be removed. 
}
  \item{amount}{
    For \code{shift.word.indices}, an integer scalar by which to shift
  the vocabulary in the corpus.   \code{amount} will be added to each
  entry of the word field in the corpus.
}
}
\value{
  A corpus with the documents merged/words filtered/words shifted.  The format of the
  input and output corpora is described in \code{\link{lda.collapsed.gibbs.sampler}}.
}
\author{
  Jonathan Chang (\email{slycoder@gmail.com})    
}
\seealso{
  \code{\link{lda.collapsed.gibbs.sampler}} for the format of
  the return value.

  \code{\link{word.counts}} to compute statistics associated with a
  corpus.
}
\examples{
data(cora.documents)

## Just use a small subset for the example.
corpus <- cora.documents[1:6]
## Get the word counts.
wc <- word.counts(corpus)

## Only keep the words which occur more than 4 times.
filtered <- filter.words(corpus,
                         as.numeric(names(wc)[wc <= 4]))
## [[1]]
##      [,1] [,2] [,3] [,4] [,5]
## [1,]    1   23   34   37   44
## [2,]    4    1    3    4    1
##
## [[2]]
##      [,1] [,2]
## [1,]   34   94
## [2,]    1    1
## ... long output ommitted ...

## Shift the second half of the corpus.
shifted <- shift.word.indices(filtered[4:6], 100)
## [[1]]
##      [,1] [,2] [,3]
## [1,]  134  281  307
## [2,]    2    5    7
##
## [[2]]
##      [,1] [,2]
## [1,]  101  123
## [2,]    1    4
##
## [[3]]
##      [,1] [,2]
## [1,]  101  194
## [2,]    6    3

## Combine the unshifted documents and the shifted documents.
concatenate.documents(filtered[1:3], shifted)
## [[1]]
##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
## [1,]    1   23   34   37   44  134  281  307
## [2,]    4    1    3    4    1    2    5    7
##
## [[2]]
##      [,1] [,2] [,3] [,4]
## [1,]   34   94  101  123
## [2,]    1    1    1    4
##
## [[3]]
##      [,1] [,2] [,3] [,4] [,5] [,6]
## [1,]   34   37   44   94  101  194
## [2,]    4    1    7    1    6    3

}
\keyword{ utilities }

