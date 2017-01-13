#############################################

library("Matrix")
library("geigen")
library("rARPACK")
library(maps)
library(WDI)
library(RColorBrewer)
library("maptools")

source("Preprocess.R")
source("SpectralClustering.R")
source("Postprocess.R")

######################################################
# speCluster()
speCluster <- function(data, conMatrix, cluster.number, 
                       iter.max=400, repetition= 400 ){
  # Perform Spectral Clustering on a data matrix 
  #
  # Args:
  #     data: A numeric data frame or matrix.
  #     conMatrix: Contiguity matrix.
  #     cluster.number: The number of clusters.
  #     iter.max: The maximum number of iterations allowed for 
  #               kmeans step.
  #     repetition: How  many  random  sets  should  be  chosen 
  #                 for  as  the  initial centers in kmeans step.
  #
  # Returns:
  #     A list contains two parts:
  #        clusters: A vector of integers(from 1:cluster.number)
  #                  indicating the cluster to which each point is 
  #                  allocated.
  #        SS: A list with two values SSW for Sum Squered Within and  
  #                  SSB for SumSquered Between
  # Error handeling
  ############################################
  #Preprocess
  outId <-outlierDetector(data)
  dataAfterPC <- prinComp(data=data,outId=outId)
  rm(data)
  ############################################
  # Spectral clustering Algorithm
  S <- similarity(data = dataAfterPC , neighbors=conMatrix)
  rm(outId, conMatrix)
  U <- produceU( similarity = S, ncol=cluster.number)
  rm(S)
  clusters <- kmeansU(data=U, cluster.number = cluster.number,iter.max=500)
  ############################################
  #postprocess
  SS <- sumSquares(data=dataAfterPC, clusters= clusters)
  ############################################
  out <- list(clusters= clusters,SS= SS)
  return(out)
}

stepOne <- function(data, conMatrix, ncol){
  # This function Computes the data after Principal component 
  #  
  #
  # Args:
  #     data: A numeric data frame or matrix.
  #     conMatrix: Contiguity matrix.
  #     ncol: number of columns of the output matrix U
  #               
  #
  # Returns:
  #     A list contains two parts:
  #     dataAfterPC: After Principal component data
  #     U: n by ncol numeric matrix that contains the ncol tops 
  #       eigenvectors of Laplacian matrix as column.
  #
  # Error handeling
  ############################################
  #Preprocess
  outId <-outlierDetector(data)
  dataAfterPC <- prinComp(data=data,outId=outId)
  rm(data)
  ############################################
  # Spectral clustering Algorithm
  S <- similarity(data = dataAfterPC, neighbors=conMatrix)
  rm(outId, conMatrix)
  U <- produceU( similarity = S, ncol=ncol)
  out <- list( dataAfterPC=dataAfterPC, U=U)
  return(out)
}
stepTwo <- function(data, U, cluster.number= cluster.number,
                    iter.max=400, repetition=400){
  # Perform Spectral Clustering on U matrix.
  #
  # Args:
  #     data: A numeric data frame or matrix.
  #     U: A numeric matrix 
  #     cluster.number: The number of clusters.
  #     iter.max: The maximum number of iterations allowed for 
  #               kmeans step.
  #     repetition: How  many  random  sets  should  be  chosen 
  #                 for  as  the  initial centers in kmeans step.
  #
  # Returns:
  #     A list contains two parts:
  #        clusters: A vector of integers(from 1:cluster.number)
  #                  indicating the cluster to which each point is 
  #                  allocated.
  #        SS: A list with two values SSW for Sum Squered Within and  
  #                  SSB for SumSquered Between
  # Error handeling
  ############################################
  
  clusters <- kmeansU(data=U, cluster.number = cluster.number,
                      iter.max=iter.max, repetition=repetition)
  SS <- sumSquares(data=data, clusters= clusters)
  ############################################
  out <- list(clusters= clusters,SS= SS)
  return(out)
}


