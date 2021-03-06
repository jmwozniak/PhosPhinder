\name{extractBackground}
\alias{extractBackground}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Background Motif Extraction
}
\description{
extractBackground() extracts a list of background sequences (of a particular width centered around a specific character) from an input file for subsequent use with tools such as motif-x.
}
\usage{
extractBackground(sequences, central, width)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sequences}{
  list of sequences from which background will be extracted.
}
  \item{central}{
  a character denoting the central residue for the background you want to extract (eg. S for Serine).
}
  \item{width}{
  an integer denoting the width of the motif you want to extract.
}
}
\details{
output will be limited to background sequences that are at least the size of the specified width. Background sequences that are shorter than this width (due to being near termini of proteins) will be removed.}
\value{
a list of strings
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Jacob M. Wozniak (jakewozniak@gmail.com)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
require(PhosPhinder)

examples.path = system.file("extdata", package = "PhosPhinder")
parsed_ref.path = paste(examples.path, "/Human_Uniprot_Parsed_Example.txt", sep="")

parsed_ref <- read.table(parsed_ref.path, header=FALSE, row.names=NULL, sep="\t")

s <- unlist(parsed_ref[,2])
c <- "S"
w <- 15

extractBack_Example <- extractBackground(s, c, w)

fileName2 <- paste(examples.path, "/Background_", c, "_" ,w, "_", Sys.Date(), ".txt", sep="")

write.table(extractBack_Example, fileName2, row.names=FALSE, col.names=FALSE, quote=FALSE)
}