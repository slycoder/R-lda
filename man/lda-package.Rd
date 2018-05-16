\name{lda-package}
\alias{lda-package}
\alias{lda}
\docType{package}
\title{
\packageTitle{lda}
}
\description{
\packageDescription{lda}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{lda}
\packageIndices{lda}
}
\author{
  \packageAuthor{lda}

  Maintainer: \packageMaintainer{lda}

  Special thanks to the following for their reports and comments:
  Edo Airoldi,
  Jordan Boyd-Graber,
  Christopher E. Cramer,
  Andrew Dai,
  James Danowski,
  Khalid El-Arini,
  Roger Levy,
  Solomon Messing,
  Joerg Reichardt,
  Dmitriy Selivanov
}
\references{
  \cite{Blei, David M. and Ng, Andrew and Jordan, Michael. Latent Dirichlet allocation. Journal of Machine Learning Research, 2003.}
}
\keyword{ package }
\seealso{
Functions to fit models:
\code{\link{lda.collapsed.gibbs.sampler}}
\code{\link{slda.em}}
\code{\link{mmsb.collapsed.gibbs.sampler}}
\code{\link{nubbi.collapsed.gibbs.sampler}}
\code{\link{rtm.collapsed.gibbs.sampler}}

Functions to read/create corpora:
\code{\link{lexicalize}}
\code{\link{read.documents}}
\code{\link{read.vocab}}

Functions to manipulate corpora:
\code{\link{concatenate.documents}}
\code{\link{filter.words}}
\code{\link{shift.word.indices}}
\code{\link{links.as.edgelist}}

Functions to compute summary statistics on corpora:
\code{\link{word.counts}}
\code{\link{document.lengths}}

Functions which use the output of fitted models:
\code{\link{predictive.distribution}}
\code{\link{top.topic.words}}
\code{\link{top.topic.documents}}
\code{\link{predictive.link.probability}}

Included data sets:
\code{\link{cora}}
\code{\link{poliblog}}
\code{\link{sampson}}
}
\examples{
## See demos for the following three common use cases:

\dontrun{demo(lda)}

\dontrun{demo(slda)}

\dontrun{demo(mmsb)}

\dontrun{demo(rtm)}
}
