######################################################
#File description comment, including purpose of program, inputs
# and outputs

######################################################
neighborMatrix <- function(NB,conFactor=1){
  # Compute constraint Matrix 
  # 
  # Args:
  #     NB: The contiguity constriant data frame
  #     conFactor: contiguity constraint factor
  # Returns:
  #     conMatrix: Contiguity Matrix
  # Error handeling
  
  ####################################################
  conMatrix <- sparseMatrix(NB[,1],NB[,2],x=rep(1,nrow(NB)))
  if(conFactor >=2){
    n <- nrow(conMatrix)
    nb <- sparseMatrix(dims =c(nrow(conMatrix),ncol(conMatrix)),i={},j={})
    #saprse identity matrix
    NBt <- sparseMatrix(1:n,1:n, x= rep(1,n))
    for(i in 1:conFactor){
      NBt <- NBt%*%conMatrix
      nb <- nb + NBt
    }
    nb[nb!=0]<-1
    diag(nb)<-0
    conMatrix <- nb
  }
  return(conMatrix) 
}
outlierDetector <- function(data, outlier.Threshold = 0.2 ){
  # Compute the outlier of the data using Principal component
  # 
  # Args:
  #     data: a numeric data frame or matrix
  #     outlier.Threshold: The Threshold which makes a data outlier. 
  #
  # Returns:
  #    outId: A logical vecotor which specifies all the outliers.
  # 
  # Error handeling
  
  ####################################################
  #Principal Component
  pc <- prcomp(data,scale=TRUE, center=TRUE)
  var <-pc$sdev^2
  cvar <-var/sum(var)
  n <- 1
  s <- cvar[n]
  while( s < 0.85){
    n <- n+1
    s <- s +cvar[n]
  }
  dataNew <- pc$x[,1:n]
  rm("pc","data") 
 
  #Similarity calculation
  dist <- as.matrix(dist(dataNew) )
  sigma <- median(dist)
  dist <- exp(-dist^2/(2*sigma^2))
  
  #Outliers detection
  diag(dist)<-0
  i <- apply(dist,1,max)
  outId <-i <= outlier.Threshold
  return(outId)
}
prinComp <- function(data, outId, showPC = F){
  # Run the pricnipal componenet algroithm on the data
  # to reduce dimension 
  #
  # Args:
  #     data: a numeric data frame or matrix
  #     outId: A logical vecotor which specifies 
  #           all the outliers.
  #     showPC: A logical value indicating whether 
  #             principal compunent should be return
  #             or not.
  #
  # Returns:
  #     dataNew: After Principal component data
  # 
  # Error handeling
  ####################################################
  outSize <-sum(outId)
  if(outSize!=0){
    colmean <- apply(data,2,mean)
    data[outId,] <- matrix(colmean, nrow=outSize, ncol= length(colmean), byrow=T )
  }
  rm(colmean,outSize)
  pc <- prcomp(data,scale=TRUE, center=TRUE)
  var <-pc$sdev^2
  cvar <-var/sum(var)
  n <- 1
  s <- cvar[n]
  while( s < 0.85){
    n <- n+1
    s <- s +cvar[n]
  }
  dataNew <- pc$x[,1:n]
  if(showPC){
    return(pc)
  }
  return(dataNew)
}
