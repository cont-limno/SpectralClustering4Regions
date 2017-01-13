###############################################
#load spectral clustering functions
source("main.R")
###############################################
# Read data
dataTerr <- read.csv("terrData.csv", header=T)
dataFW <- read.csv("freshData.csv", header=T)
i <- which(colnames(dataFW)=="hu12_states")
dataTerrFW <- merge(dataTerr, dataFW[-i], by.x="zoneid", by.y="zoneid")
NB18876 <- read.csv("NB_18876.csv", header=T)
islands <- read.csv("islandIdx.csv", header=T)
latLong18876 <- read.csv("latLong18876.csv", header=T)

rm(i)
##################################################
# function definition
NBindex <- function(index){
  id <- which(index)
  NB <- data.frame()
  for( i in 1:nrow(NB18876) ){
    if( (NB18876[i,"row"] %in% id) & 
          (NB18876[i,"neighbor"] %in% id)){
     NB <- rbind(NB,NB18876[i,c("row","neighbor")])
    }  
  }
  hash <- 1:length(id)
  names(hash) <- id
  for(i in 1:nrow(NB)){
    NB[i,1] <- hash[as.character(NB[i,1])]
    NB[i,2] <- hash[as.character(NB[i,2])]
  }
  return( NB)
}

generateData <-function(type, islandsIn =F , states= vector(), conFactor=1 ){
  # Generate the data for clustering
  # 
  # 
  # Args:
  #     type: three options "dataTerr", "dataFW", and "dataTerrFW"
  #     islandIn: if T the islands will be included
  #     states: a vector of states names that have to be included           
  #     conFactor:    contiguity constraint factor 
  # Returns:
  #    out: a list with three elements: data, conMatrix, and latLong
  # 
  # Error handeling
  #type
  if(!identical(type, dataFW)& 
       !identical(type, dataTerrFW) & 
       !identical(type, dataTerr) ){
    stop("Wrong input for type variable: Make sure to choose one of
          dataFW, dataTerr, or dataTerrFW
         ")
  }
  #islandsIn
  if(!is.logical(islandsIn)){
    stop("islandsIn must be logical variable.")
  }
  #states
  allStates <- levels(type$hu12_states)
  allStates <- allStates[nchar(allStates)==2]
  if(!is.vector(states)){
    print(allStates)
    stop("The state variable must be a vectore containing a subset of
         the above state list")
  }
  if(sum(states %in% allStates)!= length(states)){
    stop("The state variable must be a vectore containing a subset of
         the above state list")
  }
  
  ####################################################  
  
  
  #############################################
  # finding the row index
  index <- rep(T, nrow(type))
  # make index of islands False if islands are not included
  if(islandsIn==F){
    index [c(islands)$x] <- F
  }
  
  if(length(states)>0){
    id <- rep(F, nrow(type))
    outStates <- allStates[!allStates %in% states]
    for(i in 1:length(outStates)){
      state <- outStates[i]
      id <- grepl(state,type$hu12_states)|id
    }
    index <- index & !id 
  }
  
  #############################################
  # generate the output data
  data <- type[index,-c(1,2)]
  n <- nrow(data)
  m <- ncol(data)
  data <- as.matrix(data)
  data <- as.numeric(data)
  data <- matrix(data,nrow=n,ncol=m)
  # Delete the constant columns
  colSum <- apply(data,2,var)
  constants <- which(colSum==0)
  if(length(constants)!=0){
    data <- data[,-constants] 
  }
  # 
  
  latLong <- latLong18876[index,]
  NB <- NBindex(index)
  conMatrix <- neighborMatrix(NB,conFactor = conFactor)
  out <- list(data=data, conMatrix = conMatrix, latLong = latLong)
  return(out)
}


# Second Load sepectralClustering file

save.image(file="HU12SpectralClustering.RData")
 