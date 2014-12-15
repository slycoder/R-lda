\name{lda.collapsed.gibbs.sampler}
\Rdversion{1.1}
\alias{lda.collapsed.gibbs.sampler}
\alias{slda.em}
\alias{mmsb.collapsed.gibbs.sampler}
\alias{lda.cvb0}

\title{
  Functions to Fit LDA-type models
}
\description{
  These functions use a collapsed Gibbs sampler to fit three different
  models: latent Dirichlet allocation (LDA), the mixed-membership stochastic
  blockmodel (MMSB), and supervised LDA (sLDA).  These functions take
  sparsely represented input documents, perform inference, and return
  point estimates of the latent parameters using the state at the last
  iteration of Gibbs sampling. Multinomial logit for sLDA is supported 
  using the multinom function from nnet package .
}
\usage{
lda.collapsed.gibbs.sampler(documents, K, vocab, num.iterations, alpha,
eta, initial = NULL, burnin = NULL, compute.log.likelihood = FALSE,
  trace = 0L, freeze.topics = FALSE)

slda.em(documents, K, vocab, num.e.iterations, num.m.iterations, alpha,
eta, annotations, params, variance, logistic = FALSE, lambda = 10,
regularise = FALSE, method = "sLDA", trace = 0L, MaxNWts=3000)

mmsb.collapsed.gibbs.sampler(network, K, num.iterations, alpha,
beta.prior, initial = NULL, burnin = NULL, trace = 0L)

lda.cvb0(documents, K, vocab, num.iterations, alpha, eta, trace = 0L)
}

\arguments{
  \item{documents}{
    A list whose length is equal to the number of documents, D.   Each
  element of \var{documents} is an integer matrix with two rows.  Each
  column of \var{documents[[i]]} (i.e., document \eqn{i}) represents a
  word occurring in the document.

  \var{documents[[i]][1, j]} is a
  0-indexed word identifier for the jth word in document i.  That is,
  this should be an index - 1 into \var{vocab}.  \var{documents[[i]][2,
  j]} is an integer specifying the number of times that word appears in
  the document.
}
  \item{network}{
    For \code{mmsb.collapsed.gibbs.sampler}, a \eqn{D \times D} 
    matrix (coercible as logical) representing the adjacency matrix for
    the network.  Note that elements on the diagonal are ignored.
}
  \item{K}{
    An integer representing the number of topics in the model.
}
  \item{vocab}{
    A character vector specifying the vocabulary words associated with
  the word indices used in \var{documents}. 
}
  \item{num.iterations}{
    The number of sweeps of Gibbs sampling over the entire corpus to make.
}
  \item{num.e.iterations}{
    For \code{slda.em}, the number of Gibbs sampling sweeps to make over
    the entire corpus for each iteration of EM. 
}
  \item{num.m.iterations}{
    For \code{slda.em}, the number of EM iterations to make.
}
  \item{alpha}{
    The scalar value of the Dirichlet hyperparameter for
    topic proportions. 
  }
  \item{beta.prior}{
    For \code{mmsb.collapsed.gibbs.sampler}, the the beta hyperparameter
    for each entry of the block relations matrix.  This parameter should
    be a length-2 list whose entries are \eqn{K \times K} matrices.  The
    elements of the two matrices comprise the two parameters for each beta variable. 
  }
  \item{eta}{
    The scalar value of the Dirichlet hyperparamater for topic
    multinomials.
  }
  \item{initial}{
    A list of initial topic assignments for words.  It should be
  in the same format as the \var{assignments} field of the return
  value.  If this field is NULL, then the sampler will be initialized
  with random assignments.
  }
  \item{burnin}{
    A scalar integer indicating the number of Gibbs sweeps to consider
    as burn-in (i.e., throw away) for \code{lda.collapsed.gibbs.sampler}
    and \code{mmsb.collapsed.gibbs.sampler}.  If this parameter is non-NULL, it
    will also have the side-effect of enabling the
    \var{document_expects} field of the return value (see below for
    details).  Note that burnin iterations do NOT count towards \var{num.iterations}.
  }
  \item{compute.log.likelihood}{
    A scalar logical which when \code{TRUE} will cause the sampler to
  compute the log likelihood of the words (to within a constant
  factor) after each sweep over the variables.  The log likelihood for each
  iteration is stored in the \var{log.likelihood} field of the result.
  This is useful for assessing convergence, but slows things down a tiny
  bit.}
  \item{annotations}{
    A length D numeric vector of covariates associated with each
    document.  Only used by \code{slda.em} which models documents along
    with numeric annotations associated with each document. When using the 
    logistic option, annotations must be consecutive integers starting from 0.
}
  \item{params}{
    For \code{slda.em}, a length Kx(number of classes-1) numeric vector of 
    regression coefficients at which the EM algorithm should be initialized.
}
  \item{variance}{
    For \code{slda.em}, the variance associated with the Gaussian
    response modeling the annotations in \var{annotations}.
}
  \item{logistic}{
    For \code{slda.em}, a scalar logical which, when \code{TRUE}, causes
    the annotations to be modeled using a logistic response instead of a
    Gaussian (the covariates must be consecutive integers starting from 
    zero when used with sLDA).  
}
  \item{lambda}{
    When \var{regularise} is \code{TRUE}. This is a scalar that is the
    standard deviation of the Gaussian prior on the regression coefficients.
}
  \item{regularise}{
    When \code{TRUE}, a Gaussian prior is used for the regression
    coefficients. This requires the \code{penalized} package.
}
  \item{method}{
    For \code{slda.em}, a character indicating how to model the
  annotations.  Only \code{"sLDA"}, the stock model given in the
  references, is officially supported at the moment.
}
  \item{trace}{
    When \code{trace} is greater than zero, diagnostic messages will be
  output.  Larger values of \code{trace} imply more messages.
}
 \item{MaxNWts}{
    Input to the nnet's multinom function with a default value of 3000 maximum 
    weights. Increasing this value may be necessary when using logistic sLDA with
    a large number of topics at the necessary expense of longer run times. 
}
  \item{freeze.topics}{
    When \code{TRUE}, topic assignments will occur but the counts of
  words associated with topics will not change. \var{initial} should be
  set when this option is used. This is best use for sampling test
  documents.
}

}
\value{
A fitted model as a list with the following components:
  \item{assignments}{A list of length D.  Each element of the list, say
  \code{assignments[[i]]} is an integer vector of the same length as the
  number of columns in \code{documents[[i]]} indicating the topic
  assignment for each word.}  
  \item{topics}{A \eqn{K \times V} matrix where each entry indicates the
    number of times a word (column) was assigned to a topic (row).  The column
    names should correspond to the vocabulary words given in \var{vocab}.}
  \item{topic_sums}{A length K vector where each entry indicates the
    total number of times words were assigned to each topic.}
  \item{document_sums}{A \eqn{K \times D} matrix where each entry is an
    integer indicating the number of times words in each document
    (column) were assigned to each topic (column).}
  \item{log.likelihoods}{Only for \code{lda.collapsed.gibbs.sampler}.  A
    matrix with 2 rows and \code{num.iterations} columns of log likelihoods when the flag
    \code{compute.log.likelihood} is set to \code{TRUE}.  The first row
    contains the full log likelihood (including the prior), whereas the
    second row contains the log likelihood of the observations
    conditioned on the assignments.}
  \item{document_expects}{This field only exists if \var{burnin} is
    non-NULL. This field is like document_sums but instead of only
    aggregating counts for the last iteration, this field aggegates
    counts over all iterations after burnin.}  
  \item{net.assignments.left}{Only for
    \code{mmsb.collapsed.gibbs.sampler}.  A \eqn{D \times D} integer matrix of
    topic assignments for the source document corresponding to the link
    between one document (row) and another (column).}
  \item{net.assignments.right}{Only for
    \code{mmsb.collapsed.gibbs.sampler}.  A \eqn{D \times D} integer matrix of
    topic assignments for the destination document corresponding to the link
    between one document (row) and another (column).}
  \item{blocks.neg}{Only for
    \code{mmsb.collapsed.gibbs.sampler}.  A \eqn{K \times K} integer
    matrix indicating the number of times the source of a non-link was
    assigned to a topic (row) and the destination was assigned to
    another (column).}
  \item{blocks.pos}{Only for
    \code{mmsb.collapsed.gibbs.sampler}.  A \eqn{K \times K} integer
    matrix indicating the number of times the source of a link was
    assigned to a topic (row) and the destination was assigned to
    another (column).}
  \item{model}{For \code{slda.em}, a model of type \code{\link{lm}},
    the regression
    model fitted to the annotations.}
  \item{coefs}{For \code{slda.em}, a length Kx(number of classes-1) 
   numeric vector of coefficients for the regression model.}
}
\references{
  \cite{Blei, David M. and Ng, Andrew and Jordan, Michael. Latent
    Dirichlet allocation. Journal of Machine Learning Research, 2003.}

  \cite{Airoldi , Edoardo M.  and Blei, David M. and Fienberg, Stephen
    E. and Xing, Eric P.  Mixed Membership Stochastic
    Blockmodels. Journal of Machine Learning Research, 2008.}

  \cite{Blei, David M. and McAuliffe, John.  Supervised topic models.
    Advances in Neural Information Processing Systems, 2008.}

  \cite{Griffiths, Thomas L. and Steyvers, Mark.  Finding scientific
    topics.  Proceedings of the National Academy of Sciences, 2004.}

  \cite{Asuncion, A., Welling, M., Smyth, P., and Teh, Y. W.  On
  smoothing and inference for topic models.  Uncertainty in Artificial Intelligence,
  2009.}
}
\author{
  Jonathan Chang (\email{slycoder@gmail.com})
}
\note{
  WARNING: This function does not compute precisely the correct thing
    when the count associated with a word in a document is not 1 (this
    is for speed reasons currently).  A workaround when a word appears
    multiple times is to replicate the word across several columns of a
    document.  This will likely be fixed in a future version.
}

\seealso{
 \code{\link{read.documents}} and \code{\link{lexicalize}} can be used
    to generate the input data to these models.

    \code{\link{top.topic.words}},
    \code{\link{predictive.distribution}}, and \code{\link{slda.predict}} for operations on the fitted models. 
}
\examples{
## See demos for the three functions:

\dontrun{demo(lda)}

\dontrun{demo(slda)}

\dontrun{demo(mmsb)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }

