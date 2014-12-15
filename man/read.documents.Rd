\name{read.documents}
\Rdversion{1.1}
\alias{read.documents}
\alias{read.vocab}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read LDA-formatted Document and Vocabulary Files
}
\description{
  These functions read in the document and vocabulary files associated
  with a corpus.  The format of the files is the same as that used by
  LDA-C (see below for details).  The return value of these functions
  can be used by the inference procedures defined in the \pkg{lda} package.
}
\usage{
read.documents(filename = "mult.dat")

read.vocab(filename = "vocab.dat")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
    A length-1 character vector specifying the path to the
    document/vocabulary file.  These are set to \file{mult.dat} and
    \file{vocab.dat} by default.
}
}
\details{
  The details of the format are also described in the readme for LDA-C.

  The format of the documents file is appropriate for typical text data
  as it sparsely encodes observed features.  A single file encodes a
  \dfn{corpus} (a collection of documents).  Each line of the file
  encodes a single \dfn{document} (a feature vector).

  The line encoding a document begins with an integer followed by a
  number of \dfn{feature-count pairs}, all separated by spaces.  A
  feature-count pair consists of two integers separated by a colon.  The
  first integer indicates the feature (note that this is zero-indexed!)
  and the second integer indicates the count (i.e., value) of that
  feature.  The initial integer of a line indicates how many
  feature-count pairs are to be expected on that line.

  Note that we permit a feature to appear more than once on a line, in
  which case the value for that feature will be the sum of all instances
  (the behavior for such files is undefined for LDA-C).  For example, a
  line reading \samp{4 7:1 0:2 7:3 1:1}
  will yield a document with feature 0 occurring twice, feature 1
  occurring once, and feature 7 occurring four times, with all other
  features occurring zero times.

  The format of the vocabulary is a set of newline separated strings
  corresponding to features.  That is, the first line of the vocabulary
  file will correspond to the label for feature 0, the second for
  feature 1, etc.
}

\value{
  \code{read.documents} returns a list of matrices suitable as input for
  the inference routines in \pkg{lda}.  See
  \code{\link{lda.collapsed.gibbs.sampler}} for details.

  \code{read.vocab} returns a character vector of strings corresponding to
  features. 
}
\references{
  \cite{Blei, David M.  Latent Dirichlet Allocation in C. \url{http://www.cs.princeton.edu/~blei/lda-c/index.html}}
}
\author{
  Jonathan Chang (\email{slycoder@gmail.com})
}

\seealso{
  \code{\link{lda.collapsed.gibbs.sampler}} for the format of
  the return value of \code{read.documents}.

  \code{\link{lexicalize}} to generate the same output from raw text data.

  \code{\link{word.counts}} to compute statistics associated with a
  corpus.

  \code{\link{concatenate.documents}} for operations on a collection of documents.
}

\examples{
## Read files using default values.
\dontrun{setwd("corpus directory")}
\dontrun{documents <- read.documents()}
\dontrun{vocab <- read.vocab()}

## Read files from another location.
\dontrun{documents <- read.documents("corpus directory/features")}
\dontrun{vocab <- read.vocab("corpus directory/labels")}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ file }
