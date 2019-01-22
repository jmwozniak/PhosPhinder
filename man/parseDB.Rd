\name{parseDB}
\alias{parseDB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parse a .fasta Database
}
\description{
parseDB() parses a .fasta database from Uniprot into a data table containing each protein accession ID in one column and the protein sequence in another column.
}
\usage{
parseDB(database)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{database}{
a .fasta database downloaded from Uniprot.
}
}
\details{
databases must be downloaded from Uniprot or they may not work. File extension of database should be changed from ".fasta" to ".txt" prior to use. This function was built to organize .fasta databases for easier manipulation in R (such as prior to input into the phindPhospho() function) or other analysis software.
}
\value{
a data table with two columns: protein accession ID and protein sequence
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
uniprot_ref.path = paste(examples.path, "/Human_Uniprot_Example.txt", sep="")

uniprot_ref = readLines(uniprot_ref.path)

parseDB_Example <- parseDB(uniprot_ref)

filename1 <- paste(examples.path, "/Human_Uniprot_Parsed_", Sys.Date(), ".txt", sep="")

write.table(parseDB_Example, filename1, row.names=FALSE, col.names=FALSE, sep="\t")
}

