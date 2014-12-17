\name{rtm.collapsed.gibbs.sampler}
\Rdversion{1.1}
\alias{rtm.collapsed.gibbs.sampler}
\alias{rtm.em}

\title{
  Collapsed Gibbs Sampling for the Relational Topic Model (RTM).
}
\description{
  Fit a generative topic model which accounts for both the words which
  occur in a collection of documents as well as the links between the documents.
}
\usage{
rtm.collapsed.gibbs.sampler(documents, links, K, vocab, num.iterations,
  alpha, eta, beta, trace = 0L, test.start = length(documents) + 1L)
rtm.em(documents, links, K, vocab, num.e.iterations, num.m.iterations,
        alpha, eta,
        lambda = sum(sapply(links, length))/(length(links) * (length(links) -1)/2),
  initial.beta = rep(3, K), trace = 0L,
  test.start = length(documents) + 1L, tempering = 0.0)
}

\arguments{
  \item{documents}{
    A collection of documents in LDA format.  See
    \code{\link{lda.collapsed.gibbs.sampler}} for details.
  }
  \item{links}{
    A list representing the connections between the documents.  This
  list should be of the same length as the \var{documents}.  Each
  element, \code{links[[i]]}, is an integer vector expressing connections
  between document \var{i} and the 0-indexed documents pointed to by the
  elements of the vector.
}
  \item{K}{
  A scalar integer indicating the number of latent topics for the model.
}
  \item{vocab}{
    A character vector specifying the vocabulary words associated with
    the word indices used in \var{documents}.
}
  \item{num.iterations}{
    The number of sweeps of Gibbs sampling over the entire corpus to make.
}
  \item{num.e.iterations}{
    For \code{rtm.em}, the number of iterations in each Gibbs sampling E-step.
}
  \item{num.m.iterations}{
    For \code{rtm.em}, the number of M-step iterations.
}
  \item{alpha}{
    The scalar value of the Dirichlet hyperparameter for
    topic proportions.
}
  \item{eta}{
    The scalar value of the Dirichlet hyperparamater for topic
    multinomials.
}
\item{beta}{
    A length \code{K} numeric of regression coefficients expressing the
  relationship between each topic and the probability of link.
}
\item{lambda}{
  For \code{rtm.em}, the regularization parameter used when estimating
  beta.  \var{lambda} expresses the number of non-links to simulate
  among all
  possible connections between documents.
}
  \item{initial.beta}{
    For \code{rtm.em}, an initial value of \code{beta} at which to start
  the EM process.
  }
  \item{trace}{
    When \code{trace} is greater than zero, diagnostic messages will be
  output.  Larger values of \code{trace} imply more messages.
}
  \item{test.start}{
    Internal use only.
  }
  \item{tempering}{
    A numeric between 0 and 1 indicating how newly computed parameters should
    be averaged with the previous iterations parameters.  By default, the new
    values are used directly and the old value discarded.  When set to 1, the
    new values are ignored and the initial values retained indefinitely.
  }
}
\details{
  The Relational Topic Model uses LDA to model the content of
  documents but adds connections between documents as dependent on the
  similarity of the distribution of latent topic assignments.  (See
  reference for details).

  Only the exponential link probability function
  is implemented here.  Note that the collapsed Gibbs sampler is
  different than the variational inference procedure proposed in the
  paper and is extremely experimental.

  \code{rtm.em} provides an EM-wrapper around
  \code{rtm.collapsed.gibbs.sampler} which iteratively estimates the
  regression parameters \code{beta}.
}
\value{
  A fitted model as a list with the same components as returned by
  \code{\link{lda.collapsed.gibbs.sampler}}.
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
  See \code{\link{lda.collapsed.gibbs.sampler}} for a description of the
  input formats and similar models.

  \code{\link{nubbi.collapsed.gibbs.sampler}} is a different kind of
  model for document networks.

  \code{\link{predictive.link.probability}} makes predictions based on
  the output of this model.
}
\examples{
## See demo.

\dontrun{demo(rtm)}
}
\keyword{ models }

