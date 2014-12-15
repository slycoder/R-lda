\name{predictive.distribution}
\Rdversion{1.1}
\alias{predictive.distribution}
\title{
  Compute predictive distributions for fitted LDA-type models.
}
\description{
  This function takes a fitted LDA-type model and computes a predictive
  distribution for new words in a document.  This is useful for making
  predictions about held-out words.
}
\usage{
predictive.distribution(document_sums, topics, alpha, eta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{document_sums}{
    A \eqn{K \times D} matrix where each entry is a numeric proportional
    to the probability of seeing a topic (row) conditioned on document
    (column) (this entry is sometimes denoted \eqn{\theta_{d,k}} in the
    literature, see details).  Either the \var{document_sums} field or
    the \var{document_expects} field from the output of
    \code{\link{lda.collapsed.gibbs.sampler}} can be used.
  }
  \item{topics}{
    A \eqn{K \times V} matrix where each entry is a numeric proportional
    to the probability of seeing the word (column) conditioned on topic
    (row) (this entry is sometimes denoted \eqn{\beta_{w,k}} in the
    literature, see details).  The column names should correspond to the
    words in the vocabulary.  The \var{topics} field from the output of
    \code{\link{lda.collapsed.gibbs.sampler}} can be used.
  }
  \item{alpha}{
    The scalar value of the Dirichlet hyperparameter for
    topic proportions.  See references for details.
  }
  \item{eta}{
    The scalar value of the Dirichlet hyperparamater for topic
    multinomials.  See references for details.
  }
}
\details{
  The formula used to compute predictive probability is \eqn{p_d(w) =
    \sum_k (\theta_{d, k} + \alpha) (\beta_{w, k} + \eta)}.
}
\value{
  A \eqn{V \times D} matrix of the probability of seeing a word (row) in
  a document (column).   The row names of the matrix are set to the
  column names of \var{topics}.
}
\references{
  \cite{Blei, David M. and Ng, Andrew and Jordan, Michael. Latent Dirichlet allocation. Journal of Machine Learning Research, 2003.}
}
\author{
  Jonathan Chang (\email{slycoder@gmail.com})  
}

\seealso{
  \code{\link{lda.collapsed.gibbs.sampler}} for the format of
  \var{topics} and \var{document_sums} and details of the model.

  \code{\link{top.topic.words}} demonstrates another use for a fitted
  topic matrix.
}
\examples{
## Fit a model (from demo(lda)).
data(cora.documents)
data(cora.vocab)

K <- 10 ## Num clusters
result <- lda.collapsed.gibbs.sampler(cora.documents,
                                      K,  ## Num clusters
                                      cora.vocab,
                                      25,  ## Num iterations
                                      0.1,
                                      0.1) 

## Predict new words for the first two documents
predictions <-  predictive.distribution(result$document_sums[,1:2],
                                        result$topics,
                                        0.1, 0.1)

## Use top.topic.words to show the top 5 predictions in each document.
top.topic.words(t(predictions), 5)

##      [,1]         [,2]      
## [1,] "learning"   "learning"
## [2,] "algorithm"  "paper"   
## [3,] "model"      "problem" 
## [4,] "paper"      "results" 
## [5,] "algorithms" "system"  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }
