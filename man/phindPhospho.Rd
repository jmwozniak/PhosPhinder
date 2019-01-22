\name{phindPhospho}
\alias{phindPhospho}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Localize Phospho-sites from Phospho-proteomic Data
}
\description{
phindPhospho() uses phospho-peptide information to localize the phospho-sites in the full-length protein and returns flanking sequences for subsequent use with tools such as motif-x.
}
\usage{
phindPhospho(data_table, reference_table)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_table}{
a data table containing your experimentally acquired phospho-peptide data. See "phindPhospho_Input_Example.csv" for an example.
}
  \item{reference_table}{
a data table containing full protein sequence information for the phospho-peptides in the data_table. This input can be created from any Uniprot database using the parseDB() function. See "Human_Uniprot_Parsed_Example.txt" for an example.
}
}
\details{
input files must contain the following 5 columns (exactly as below - see "phindPhospho_Input_Example.csv"):

  1) Protein_ID - Uniprot accession of protein

  2) Peptide_Seq - peptide sequence containing phospho-site(s)

  3) Total_Sites - total number of phospho-site(s) within phospho-peptide

  4) Phos_Locs - location of phospho-site(s) within phospho-peptide in the ptmRS convention. Multiple phospho-sites should be separated by a semi-colon  (eg. "S9(Phospho): 98.56; S18(Phospho): 99.84")

  5) Qualifier - unique qualifier of phospho-peptide

}
\value{
a data table with 7 columns:

  1) Protein ID - see above

  2) Qualifier - see above

  3) Pep_Loc - location(s) of phosphosite(s) on peptide

  4) Prot_Loc - location(s) of phosphosite(s) on protein

  5) Flank_Seq - flanking sequence(s) of phosphosite(s)

  6) Confidence - ambiguity of localized phosphosite(s)

  7) Prot_Seq - protein sequence

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
examples.path = system.file("extdata", package = "PhosPhinder")
phos_data.path = paste(examples.path, "/phindPhospho_Input_Example.csv", sep="")
parsed_ref.path = paste(examples.path, "/Human_Uniprot_Parsed_Example.txt", sep="")

phos_data <- read.csv(phos_data.path, header=TRUE, row.names = NULL, stringsAsFactors = FALSE)

parsed_ref <- read.table(parsed_ref.path, header=FALSE, row.names=NULL, sep="\t")

phindPhospho_Example <- phindPhospho(phos_data, parsed_ref)

fileName3 <- paste(examples.path, "/PhosPhound_", Sys.Date(), ".txt", sep="")

write.table(phindPhospho_Example, fileName3, row.names=FALSE, append = FALSE,  sep = "\t")
}

