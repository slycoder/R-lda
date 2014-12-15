\name{lexicalize}
\Rdversion{1.1}
\alias{lexicalize}
\title{
  Generate LDA Documents from Raw Text
}
\description{
  This function reads raw text in \dfn{doclines} format and returns a
  corpus and vocabulary suitable for the inference procedures
  defined in the \pkg{lda} package.
}
\usage{
lexicalize(doclines, sep = " ", lower = TRUE, count = 1L, vocab = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{doclines}{
    A character vector of document lines to be used to construct
    a corpus.  See details for a description of the format of these
    lines. 
}
  \item{sep}{
    Separator string which is used to tokenize the input strings
  (default \samp{ }).
}
  \item{lower}{
    Logical indicating whether or not to convert all tokens to
  lowercase (default \samp{TRUE}).
}
  \item{count}{
    An integer scaling factor to be applied to feature counts.  A single
    observation of a feature will be rendered as \var{count}
    observations in the return value (the default
    value, \samp{1}, is appropriate in most cases).
}
  \item{vocab}{
    If left unspecified (or \code{NULL}), the vocabulary for the corpus
    will be automatically inferred from the observed tokens.  Otherwise,
    this parameter should be a character vector specifying acceptable
    tokens.  Tokens not appearing in this list will be filtered from the
    documents. 
}
}
\details{
  This function first tokenizes a character vector by splitting each
  entry of the vector by \var{sep} (note that this is currently a fixed
  separator, not a regular expression).  If \var{lower} is \samp{TRUE},
  then the tokens are then all converted to lowercase.

  At this point, if \var{vocab} is \code{NULL}, then a vocabulary is
  constructed from the set of unique tokens appearing across all
  character vectors.  Otherwise, the tokens derived from the character
  vectors are filtered so that only those appearing in \var{vocab} are
  retained.

  Finally, token instances within each document (i.e., original
  character string) are tabulated in the format described in
  \code{\link{lda.collapsed.gibbs.sampler}}.
}
\value{
  If \var{vocab} is unspecified or \code{NULL}, a list with two components:
  \item{documents}{A list of document matrices in the format described in \code{\link{lda.collapsed.gibbs.sampler}}.}
  \item{vocab}{A character vector of unique tokens occurring in the corpus.}
}
\author{
  Jonathan Chang (\email{slycoder@gmail.com})  
}
\note{ Because of the limited tokenization and filtering capabilities of this
  function, it may not be useful in many cases.  This may be resolved
  in a future release.
}

\seealso{
  \code{\link{lda.collapsed.gibbs.sampler}} for the format of
  the return value.

  \code{\link{read.documents}} to generate the same output from a file
  encoded in LDA-C format.

  \code{\link{word.counts}} to compute statistics associated with a
  corpus.

  \code{\link{concatenate.documents}} for operations on a collection of documents.  
}
\examples{
## Generate an example.
example <- c("I am the very model of a modern major general",
             "I have a major headache")

corpus <- lexicalize(example, lower=TRUE)

## corpus$documents:
## $documents[[1]]
##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
## [1,]    0    1    2    3    4    5    6    7    8     9
## [2,]    1    1    1    1    1    1    1    1    1     1
## 
## $documents[[2]]
##      [,1] [,2] [,3] [,4] [,5]
## [1,]    0   10    6    8   11
## [2,]    1    1    1    1    1

## corpus$lexicon:
## $vocab
## [1] "i"        "am"       "the"      "very"     "model"    "of"      
## [7] "a"        "modern"   "major"    "general"  "have"     "headache"

## Only keep words that appear at least twice:
to.keep <- corpus$vocab[word.counts(corpus$documents, corpus$vocab) >= 2]

## Re-lexicalize, using this subsetted vocabulary
documents <- lexicalize(example, lower=TRUE, vocab=to.keep)

## documents:
## [[1]]
##      [,1] [,2] [,3]
## [1,]    0    1    2
## [2,]    1    1    1
## 
## [[2]]
##      [,1] [,2] [,3]
## [1,]    0    1    2
## [2,]    1    1    1
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ file }
