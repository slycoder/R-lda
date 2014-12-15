\name{slda.predict}
\alias{slda.predict}
\alias{slda.predict.docsums}
\title{
  Predict the response variable of documents using an sLDA model.
}
\description{
    These functions take a fitted sLDA model and predict the value of
    the response variable (or document-topic sums) for each given document.  
}
\usage{
slda.predict(documents, topics, model, alpha, eta,
num.iterations = 100, average.iterations = 50, trace = 0L)

slda.predict.docsums(documents, topics, alpha, eta,
num.iterations = 100, average.iterations = 50, trace = 0L)
}
\arguments{
  \item{documents}{
    A list of document matrices comprising a corpus, in the format
    described in \code{\link{lda.collapsed.gibbs.sampler}}.
}
\item{topics}{
    A \eqn{K \times V} matrix where each entry is an integer that is the
    number of times the word (column) has been allocated to the topic
    (row) (a normalised version of this is sometimes denoted
    \eqn{\beta_{w,k}} in the literature, see details).  The column names
    should correspond to the words in the vocabulary.  The \var{topics}
    field from the output of \code{\link{slda.em}} can be used.  
}
\item{model}{
  A fitted model relating a document's topic distribution to the
  response variable.  The \var{model} field from the output of
  \code{\link{slda.em}} can be used.
}
\item{alpha}{
  The scalar value of the Dirichlet hyperparameter for
  topic proportions.  See references for details.
}
\item{eta}{
  The scalar value of the Dirichlet hyperparamater for topic
  multinomials.
}
\item{num.iterations}{
  Number of iterations of inference to perform on the documents.
}
\item{average.iterations}{
  Number of samples to average over to produce the predictions.
}
\item{trace}{
  When \code{trace} is greater than zero, diagnostic messages will be
  output.  Larger values of \code{trace} imply more messages.
}
}
\details{
  Inference is first performed on the documents by using Gibbs sampling
  and holding the word-topic matrix \eqn{\beta_{w,k}} constant.  Typically
  for a well-fit model only a small number of iterations are required to
  obtain good fits for new documents.  These topic vectors are then
  piped through \code{model} to yield numeric predictions associated
  with each document.
}
\value{
  For \code{slda.predict}, a numeric vector of the same length as
  \code{documents} giving the predictions. For \code{slda.predict.docsums}, a
  \eqn{K \times N} matrix of document assignment counts.  
}
\references{
  \cite{Blei, David M. and McAuliffe, John.  Supervised topic models.
    Advances in Neural Information Processing Systems, 2008.}
}
\author{
  Jonathan Chang (\email{slycoder@gmail.com})
}

\seealso{
  See \code{\link{lda.collapsed.gibbs.sampler}} for a description of the
  format of the input data, as well as more details on the model.
  
  See \code{\link{predictive.distribution}} if you want to make
  predictions about the contents of the documents instead of the
  response variables.
}
\examples{
## The sLDA demo shows an example usage of this function.
\dontrun{demo(slda)}
}
\keyword{ utilities }
