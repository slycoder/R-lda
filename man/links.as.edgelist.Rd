\name{links.as.edgelist}
\Rdversion{1.1}
\alias{links.as.edgelist}
\title{
  Convert a set of links keyed on source to a single list of edges.
}
\description{
  This function takes as input a collection of links (as used/described
  by the model fitting functions in this package) and reproduces the
  links as a matrix.
}
\usage{
links.as.edgelist(links)
}
\arguments{
  \item{links}{
    A list of links; the format of this is described in \code{\link{rtm.collapsed.gibbs.sampler}}.
  }
}
\value{
  A two-column matrix where each row represents an edge.   Note that the
  indices in this matrix are 1-indexed rather than 0-indexed.
}
\author{
  Jonathan Chang (\email{slycoder@gmail.com})
}
\seealso{
  \code{\link{rtm.collapsed.gibbs.sampler}} for the input format.
  \code{\link{predictive.link.probability}} is a usage example of the output
  of this function.
}
\examples{
## Take the citations for the first few documents of Cora.
data(cora.cites)

links <- cora.cites[1:5]
links
## [[1]]
## [1] 484 389

## [[2]]
## integer(0)

## [[3]]
## integer(0)

## [[4]]
## [1] 177 416 533

## [[5]]
## [1] 153

links.as.edgelist(links)
##      [,1] [,2]
## [1,]    1  485
## [2,]    1  390
## [3,]    4  178
## [4,]    4  417
## [5,]    4  534
## [6,]    5  154
}
\keyword{ utilities }
