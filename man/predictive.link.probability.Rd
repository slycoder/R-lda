\name{predictive.link.probability}
\Rdversion{1.1}
\alias{predictive.link.probability}
\title{
  Use the RTM to predict whether a link exists between two documents.
}
\description{
  This function takes a fitted LDA-type model (e.g., LDA or RTM) and
  makes predictions about the likelihood of a link existing between
  pairs of documents.
}
\usage{
predictive.link.probability(edgelist, document_sums, alpha, beta)
}
\arguments{
  \item{edgelist}{
    A two-column integer matrix where each row represents an edge on
    which to make a prediction.  An edge is expressed as a pair of
    integer indices (1-indexed) into the columns (i.e., documents) of
    \var{document_sums} (see below). 
  }
  \item{document_sums}{
    A \eqn{K \times D} matrix where each entry is a numeric proportional
    to the probability of seeing a topic (row) conditioned on document
    (column) (this entry is sometimes denoted \eqn{\theta_{d,k}} in the
    literature, see details).  The \var{document_sums} field or
    the \var{document_expects} field from the output of
    \code{\link{lda.collapsed.gibbs.sampler}} and
    \code{\link{rtm.collapsed.gibbs.sampler}} can be used.    
  }
  \item{alpha}{
    The value of the Dirichlet hyperparamter generating the distribution
  over \var{document_sums}.  This, in effect, smooths the similarity
  between documents. 
  }
  \item{beta}{
    A numeric vector of regression weights which is used to determine
    the similarity between two vectors (see details). Arguments will be
    recycled to create a vector of length \code{dim(document_sums)[1]}. 
  }
}
\details{
  Whether or not a link exists between two documents \eqn{i} and \eqn{j}
  is a function of the weighted inner product of the
  \code{document_sums[,i]} and \code{document_sums[,j]}.   After
  normalizing \code{document_sums} column-wise, this inner
  product is weighted by \var{beta}.

  This quantity is then passed to a
  link probability function.  Like
  \code{\link{rtm.collapsed.gibbs.sampler}} in this package, only the
  exponential link probability function is supported.  Note that
  quantities are automatically scaled to be between 0 and 1.
}
\value{
  A numeric vector of length \code{dim(edgelist)[1]}, representing the
  probability of a link existing between each pair of documents given in
  the edge list.  
}
\references{
  \cite{Chang, Jonathan and Blei, David M.
    Relational Topic Models for Document Networks.  Artificial
  intelligence and statistics. 2009.}
}
\author{
  Jonathan Chang (\email{slycoder@gmail.com})  
}
\seealso{
  \code{\link{rtm.collapsed.gibbs.sampler}} for the format of
  \var{document_sums}.  \code{\link{links.as.edgelist}} produces values
  for \var{edgelist}.  \code{\link{predictive.distribution}} makes
  predictions about document content instead.
}
\examples{
## See demo.

\dontrun{demo(rtm)}
}
\keyword{ utilities }
