\name{newsgroups}
\Rdversion{1.1}
\alias{newsgroup}
\alias{newsgroup.train.documents}
\alias{newsgroup.test.documents}
\alias{newsgroup.train.labels}
\alias{newsgroup.test.labels}
\alias{newsgroup.vocab}
\alias{newsgroup.label.map}
\docType{data}
\title{
  A collection of newsgroup messages with classes.
}
\description{
  The 20 Newsgroups data set is a collection of approximately 20,000 newsgroup documents, 
  partitioned (nearly) evenly across 20 different newsgroups. 
}
\usage{
data(newsgroup.train.documents)
data(newsgroup.test.documents)
data(newsgroup.train.labels)
data(newsgroup.test.labels)
data(newsgroup.vocab)
data(newsgroup.label.map)
}
\format{
  \code{newsgroup.train.documents} and \code{newsgroup.test.documents}
  comprise a corpus of 20,000 newsgroup documents conforming to the LDA format,
  partitioned into 11269 training and 7505 training and test cases evenly distributed
  across 20 classes. 

  \code{newsgroup.train.labels} is a numeric vector of length 11269 which gives
  a class label from 1 to 20 for each training document in the corpus.
  
  \code{newsgroup.test.labels} is a numeric vector of length 7505 which gives
  a class label from 1 to 20 for each training document in the corpus.
  
  \code{newsgroup.vocab} is the vocabulary of the corpus.
  
  \code{newsgroup.label.map} maps the numeric class labels to actual class names.
  
}
\source{
  \cite{http://qwone.com/~jason/20Newsgroups/}
}
\seealso{
  \code{\link{lda.collapsed.gibbs.sampler}} for the format of the
  corpus.
}
\examples{
data(newsgroup.train.documents)
data(newsgroup.test.documents)
data(newsgroup.train.labels)
data(newsgroup.test.labels)
data(newsgroup.vocab)
data(newsgroup.label.map)
}
\keyword{datasets}
