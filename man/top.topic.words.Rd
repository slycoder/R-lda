\name{top.topic.words}
\Rdversion{1.1}
\alias{top.topic.words}
\alias{top.topic.documents}
\title{
  Get the Top Words and Documents in Each Topic
}
\description{
  This function takes a model fitted using
  \code{\link{lda.collapsed.gibbs.sampler}} and returns a matrix of the
  top words in each topic.
}
\usage{
top.topic.words(topics, num.words = 20, by.score = FALSE)
top.topic.documents(document_sums, num.documents = 20, alpha = 0.1)
}
\arguments{
  \item{topics}{
    For \code{top.topic.words}, a \eqn{K \times V} matrix where each entry is a numeric proportional
    to the probability of seeing the word (column) conditioned on topic
    (row) (this entry is sometimes denoted \eqn{\beta_{w,k}} in the
    literature, see details).  The column names should correspond to the words in the
    vocabulary.  The \var{topics} field from the output of
    \code{\link{lda.collapsed.gibbs.sampler}} can be used.
  }
  \item{num.words}{
    For \code{top.topic.words}, the number of top words to return for each topic.
  }
  \item{document_sums}{
    For \code{top.topic.documents}, a \eqn{K \times D} matrix where each entry is a numeric proportional
    to the probability of seeing a topic (row) conditioned on the
    document (column) (this entry is sometimes denoted \eqn{\theta_{d,k}} in the
    literature, see details).  The \var{document_sums} field from the output of
    \code{\link{lda.collapsed.gibbs.sampler}} can be used.
  }
  \item{num.documents}{
    For \code{top.topic.documents}, the number of top documents to return for each topic.
  }
  \item{by.score}{
    If \var{by.score} is set to \code{FALSE} (default), then words in
    each topic will
    be ranked according to probability mass for each word \eqn{\beta_{w,
	k}}.  If \var{by.score} is \code{TRUE}, then words will be
    ranked according to a score defined by \eqn{\beta_{w, k} (\log
      \beta_{w,k} - 1 / K \sum_{k'} \log \beta_{w,k'})}.  
  }
  \item{alpha}{
  }
}
\value{
  For \code{top.topic.words}, a \eqn{num.words \times K} character matrix where each column contains
  the top words for that topic.  
  
  For \code{top.topic.documents}, a \eqn{num.documents \times K} integer matrix where each column contains
  the top documents for that topic.  The entries in the matrix are 
  column-indexed references into \code{document_sums}.
}
\references{
  \cite{Blei, David M. and Ng, Andrew and Jordan, Michael. Latent Dirichlet allocation. Journal of Machine Learning Research, 2003.}
}
\author{
  Jonathan Chang (\email{slycoder@gmail.com})  
}

\seealso{
  \code{\link{lda.collapsed.gibbs.sampler}} for the format of \var{topics}.

  \code{\link{predictive.distribution}} demonstrates another use for a fitted
  topic matrix.
}
\examples{
## From demo(lda).

data(cora.documents)
data(cora.vocab)

K <- 10 ## Num clusters
result <- lda.collapsed.gibbs.sampler(cora.documents,
                                      K,  ## Num clusters
                                      cora.vocab,
                                      25,  ## Num iterations
                                      0.1,
                                      0.1) 

## Get the top words in the cluster
top.words <- top.topic.words(result$topics, 5, by.score=TRUE)

## top.words:
##      [,1]             [,2]        [,3]       [,4]            [,5]      
## [1,] "decision"       "network"   "planning" "learning"      "design"  
## [2,] "learning"       "time"      "visual"   "networks"      "logic"   
## [3,] "tree"           "networks"  "model"    "neural"        "search"  
## [4,] "trees"          "algorithm" "memory"   "system"        "learning"
## [5,] "classification" "data"      "system"   "reinforcement" "systems" 
##      [,6]         [,7]       [,8]           [,9]           [,10]      
## [1,] "learning"   "models"   "belief"       "genetic"      "research" 
## [2,] "search"     "networks" "model"        "search"       "reasoning"
## [3,] "crossover"  "bayesian" "theory"       "optimization" "grant"    
## [4,] "algorithm"  "data"     "distribution" "evolutionary" "science"  
## [5,] "complexity" "hidden"   "markov"       "function"     "supported"
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }
