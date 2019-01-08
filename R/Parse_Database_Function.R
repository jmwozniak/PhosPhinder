parseDB <- function(ref){
  ##Orders Uniprot fasta file into data frame with 2 columns: ID, sequence
  no_ofSeqs <- sum(substr(ref, 0, 1) == ">")
  reftab <- rep("",2*no_ofSeqs)
  temp <- NULL
  count <- 1
  for(i in 1:length(ref)){
    if(substr(ref[i], 0, 1) == ">"){
      reftab[count] <- ref[i]
      count <- count + 1
    }
    else{
      temp <- c(temp, ref[i])
      temp <- paste(temp, sep="", collapse="")
    }

    if((substr(ref[i+1], 0, 1) == ">")||(i == length(ref)))
    {
      reftab[count] <- temp
      count <- count + 1
      temp <- NULL
    }
  }
  reftab <- matrix(reftab,ncol=2,byrow=TRUE)
  temp <- reftab

  # Substitute "|" for "]"  because R is stupid
  tempIDs <- gsub("|", "]", reftab[,1], fixed = TRUE)

  #Get all accession IDs from column 1
  accID <- rep("",no_ofSeqs)
  for(i in 1:length(tempIDs)){
    left <- gregexpr(']',tempIDs[i])[[1]][1]
    right <- gregexpr(']',tempIDs[i])[[1]][2]
    accID[i] <- substr(tempIDs[i], left+1, right-1)
  }

  #Replace column 1 with the accessionIDs
  reftab[,1] <- accID

  return(reftab)
}

parseDB_SLOW <- function(ref){
  reftab <- NULL
  temp <- NULL
  for(i in 1:length(ref)){
    if(substr(ref[i], 0,1) == ">"){
      reftab <- c(reftab, ref[i])
    }else {
      temp <- c(temp, ref[i])
      temp <- paste(temp, sep="", collapse="")
    }

    if ((substr(ref[i+1], 0,1) == ">")||(i == length(ref)))
    {
      reftab <- c(reftab, temp)
      temp<-NULL
    }
  }
  reftab <- matrix(reftab,ncol=2,byrow=TRUE)

  # Substitute "|" for "]"  because R is stupid
  tempIDs <- gsub("|", "]", reftab[,1], fixed = TRUE)

  #Get all accession IDs from column 1
  accID<-NULL
  for(i in 1:length(tempIDs)){
    left <- gregexpr(']',tempIDs[i])[[1]][1]
    right <-gregexpr(']',tempIDs[i])[[1]][2]
    accID <- c(accID, substr(tempIDs[i], left+1, right-1))
  }

  #Replace column 1 with the accessionIDs
  reftab[,1] <- accID
  return(reftab)
}
