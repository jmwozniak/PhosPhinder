phindPhospho <- function(data, reftab){
  library(stringr)

  #Format Phospho-data so that its easier to work with in R
  bestSites <- gsub(" ", "", data$Phos_Locs, fixed = TRUE)
  bestSites <- gsub("(", "<", bestSites, fixed = TRUE)
  noPeptides <- length(bestSites)

  #Extract phospho site locations from bestSites vector
  noPhosSites <- NULL
  confidence <- rep("",noPeptides)
  locList <- rep("",noPeptides)
  tempLoc<- NULL
  tempLoc2 <- NULL
  for(i in 1:noPeptides){
    noPhosSites <- str_count(bestSites[i],";")
    tempLoc <- substr(bestSites[i], 0, 3)

    if(noPhosSites > 0){
      for(y in 1:(noPhosSites)){
        temp <- gregexpr(";", bestSites[i])[[1]][y]
        tempLoc2 <- substr(bestSites[i], temp+1, temp+3)
        tempLoc <- c(tempLoc, tempLoc2)
      }
    }

    tempLoc <- paste(tempLoc, sep="", collapse = ".")
    locList[i] <- tempLoc
    tempLoc <- NULL
    if(noPhosSites + 1 == data$Total_Sites[i]) confidence[i] <- "Unambiguous" else confidence[i] <- "Ambiguous"
  }
  locList <- gsub("<", "", locList, fixed = TRUE)

  #Makes the Annotated sequences easier to work with
  seqs1 <- NULL
  seqs1 <- data$Peptide_Seq
  seqs1 <- toupper(seqs1)

  #Reformat reference table
  refID <- unlist(reftab[,1][reftab[,1] %in% data$Protein_ID])
  refSeq <- unlist(reftab[,2][reftab[,1] %in% data$Protein_ID])
  reftab2 <- c(as.character(refID),as.character(refSeq))
  reftab2 <- matrix(reftab2,ncol=2,byrow=FALSE)

  #Find seqs1 in the reference database and get new location of phosphosites in full protein
  fullProt <- NULL
  newLocList <- rep("",noPeptides)
  temp2 <- NULL
  temp3 <- NULL
  flanking <- NULL
  flankList <- rep("",noPeptides)
  rows <- nrow(reftab2)
  for(i in 1:noPeptides){
    for(y in 1:rows){
      if(data$Protein_ID[i]==reftab2[y,1])
      {
        fullProt <- c(fullProt, reftab2[y,2])
        tempLoc <- regexpr(seqs1[i], reftab2[y,2])[1]
        pepLocs <- strsplit(locList[i], ".", fixed = TRUE)

        for(x in 1:length(pepLocs[[1]])){
          temp <- as.numeric(substr(pepLocs[[1]][x], 2, nchar(pepLocs[[1]][x])))
          site <- substr(pepLocs[[1]][x], 1, 1)
          newLoc <- temp + tempLoc - 1
          flanking <- substr(fullProt[[i]], newLoc-7, newLoc+7)
          newSiteLoc <- paste(c(site, newLoc), sep="", collapse = "")
          temp2 <- c(temp2, newSiteLoc)
          temp3 <- c(temp3, flanking)
        }
        temp2 <- paste(temp2, sep="", collapse=".")
        newLocList[i] <- temp2
        temp2 <- NULL
        temp3 <- paste(temp3, sep="", collapse=".")
        flankList[i] <- temp3
        temp3 <- NULL
        break
      }
    }
  }

  #Create table with ProteinID, GeneName, PhosphoLoc on Peptide, Peptide sequence, phospholoc on Protein, and protein sequence
  writelist<- c(as.character(data$Protein_ID), as.character(data$Qualifier), as.character(locList), as.character(newLocList), as.character(flankList), as.character(confidence), as.character(fullProt))
  writetab <- matrix(writelist,ncol=7,byrow=FALSE)

  #Name the columns of the table
  colnames(writetab) <- c("Protein_ID", "Qualifier", "Pep_Loc",  "Prot_Loc", "Flank_Seq", "Confidence","Prot_Seq")

  return(writetab)
}



