\name{nubbi.collapsed.gibbs.sampler}
\Rdversion{1.1}
\alias{nubbi.collapsed.gibbs.sampler}
\title{
  Collapsed Gibbs Sampling for the Networks Uncovered By Bayesian
  Inference (NUBBI) Model.
}
\description{
  Fit a NUBBI model, which takes as input a collection of entities with
  corresponding textual descriptions as well as a set of descriptions
  for pairs of entities.  The NUBBI model the produces a latent space
  description of both the entities and the relationships between them.
}
\usage{
nubbi.collapsed.gibbs.sampler(contexts, pair.contexts, pairs, K.individual,
                              K.pair, vocab, num.iterations, alpha, eta, xi)
}
\arguments{
  \item{contexts}{
    The set of textual descriptions (i.e., documents) for individual
    entities in LDA format (see
    \code{\link{lda.collapsed.gibbs.sampler}} for details).
}
  \item{pair.contexts}{
    A set of textual descriptions for pairs of entities, also in LDA format.
  }
  \item{pairs}{
    Labelings as to which pair each element of \code{pair.contexts}
    refer to.  This parameter should be an integer matrix with two columns
    and the same number of rows as \code{pair.contexts}.  The two
    elements in each row of \code{pairs} are 0-indexed indices into
    \code{contexts} indicating which two entities that element of
    \code{pair.contexts} describes.
  }
\item{K.individual}{
  A scalar integer representing the number of topics for the individual entities.
}
\item{K.pair}{
  A scalar integer representing the number of topics for entity pairs.
}
  \item{vocab}{
    A character vector specifying the vocabulary words associated with
    the word indices used in \var{contexts} and \var{pair.contexts}.
}
  \item{num.iterations}{
    The number of sweeps of Gibbs sampling over the entire corpus to make.
}
  \item{alpha}{
    The scalar value of the Dirichlet hyperparameter for
    topic proportions.
}
  \item{eta}{
    The scalar value of the Dirichlet hyperparamater for topic
    multinomials.
}
  \item{xi}{
    The scalar value of the Dirichlet hyperparamater for source
    proportions.
}
}
\details{
  The NUBBI model is a switching model wherein the description of each
  entity-pair can be ascribed to either the first entity of the pair,
  the second entity of the pair, or their relationship.  The NUBBI model
  posits a latent space (i.e., topic model) over the individual entities, and a different
  latent space over entity relationships.

  The collapsed Gibbs sampler used in this model is different than the
  variational inference method proposed in the paper and is highly experimental.
}
\value{
  A fitted model as a list with the same components as returned by
  \code{\link{lda.collapsed.gibbs.sampler}} with the following additional components:
  \item{source_assignments}{A list of \code{length(pair.contexts)} whose
    elements \code{source_assignments[[i]]} are of the same length as
    \code{pair.contexts[[i]]} where each entry is either 0 if the
    sampler assigned the word to the first entity, 1 if the sampler
    assigned the word to the second entity, or 2 if the sampler assigned
    the word to the relationship between the two.}
  \item{document_source_sums}{A matrix with three columns and
    \code{length(pair.contexts)} rows where each row indicates how many
    words were assigned to the first entity of the pair, the second
    entity of the pair, and the relationship between the two,
    respectively.}
  \item{document_sums}{Semantically similar to the entry in
    \code{lda.collapsed.gibbs.sampler}, except that it is a list whose
    first \code{length(contexts)} correspond to the columns of the entry
    in \code{lda.collapsed.gibbs.sampler} for the individual contexts,
    and the remaining \code{length(pair.contexts)} entries correspond to
    the columns for the pair contexts.}
  \item{topics}{Like the entry in \code{lda.collapsed.gibbs.sampler},
    except that it contains the concatenation of the \code{K.individual}
    topics and the \code{K.pair} topics.}
}
\references{
  \cite{Chang, Jonathan and Boyd-Graber, Jordan and Blei, David M.
    Connections between the lines: Augmenting social networks with text.
    KDD, 2009.}
}
\author{
  Jonathan Chang (\email{slycoder@gmail.com})
}
\note{
  The underlying sampler is quite general and could potentially be used
  for other models such as the author-topic model (McCallum et al.) and the citation
  influence model (Dietz et al.).   Please examine the source code
  and/or contact the author(s) for further details.
}
\seealso{
  See \code{\link{lda.collapsed.gibbs.sampler}} for a description of the
  input formats and similar models.

  \code{\link{rtm.collapsed.gibbs.sampler}} is a different kind of
  model for document networks.
}
\examples{
## See demo.

\dontrun{demo(nubbi)}
}
\keyword{ models }
