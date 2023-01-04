#' @noRd
#' @title Functions to prepare microbiome data for unit tests
#' @description 
#' Prepares the micobiome data with or without two outliers.
#' Two different subsets of the measurement components are available.
#' @description Cleaned. TM7, Cyanobacteria/Chloroplast, Actinobacteria, Proteobacteria, other

ppi_microbiomedata_cleaned_TCAP <- function(){
  data("microdata", package = "scorecompdir")
  microdata <- microdata[!microdata$IndividualID %in% c(2079, 2280), ] #remove two outlying measurements
  countdata=as.matrix(microdata[,12:31])

  #sample size
  n=92

  #dimension
  p=20

  #calculate totals
  tot=matrix(0,n,1)
  for (j in 1:p)
  {
   tot=tot+countdata[,j]
  }
  tot=as.vector(tot)

  #proportion data
  prop=countdata
  for (j in 1:n)
  {
    	prop[j,]=countdata[j,]/tot[j]
  }

  ####Reduce dimensions to p=5####


  #calculate 5D dataset
  comb=matrix(0,n,5)
  comb[,1]=prop[,"TM7"]
  comb[,2]=prop[,"Cyanobacteria/Chloroplast"]
  comb[,3]=prop[,"Actinobacteria"]
  comb[,4]=prop[,"Proteobacteria"]
  comb[,5]=abs(1-comb[,1]-comb[,2]-comb[,4]-comb[,3])
  propreal=comb
  colnames(propreal) <- c("TM7", "Cyanobacteria/Chloroplast", "Actinobacteria", "Proteobacteria", "pool")

  #dimension
  p=5

  #set beta (this is fixed here)
  beta0=matrix(0,p,1)
  beta0[1]=-0.8
  beta0[2]=-0.85
  beta0[3]=0
  beta0[4]=-0.2
  beta0[5]=0
  return(list(
    propreal = propreal,
    beta0 = beta0,
    p = p
  ))
}

#' @noRd
#' @description Not cleaned. TM7, Cyanobacteria/Chloroplast, Actinobacteria, Proteobacteria, other
ppi_microbiomedata_TCAP <- function(){
  data("microdata", package = "scorecompdir")
  countdata=as.matrix(microdata[,12:31])
  
  #sample size
  n=94
  
  #dimension
  p=20
  
  #calculate totals
  tot=matrix(0,n,1)
  for (j in 1:p)
  {
   tot=tot+countdata[,j]
  }
  tot=as.vector(tot)
  
  #proportion data
  prop=countdata
  for (j in 1:n)
  {
  	prop[j,]=countdata[j,]/tot[j]
  }
  
  ###Reduce dimensions to p=5
  
  
  #calculate 5D dataset
  comb=matrix(0,n,5)
  comb[,1]=prop[,"TM7"]
  comb[,2]=prop[,"Cyanobacteria/Chloroplast"]
  comb[,3]=prop[,"Actinobacteria"]
  comb[,4]=prop[,"Proteobacteria"]
  comb[,5]=abs(1-comb[,1]-comb[,2]-comb[,4]-comb[,3])
  propreal=comb
  colnames(propreal) <- c("TM7", "Cyanobacteria/Chloroplast", "Actinobacteria", "Proteobacteria", "pool")
  
  
  #dimension
  p=5
  
  return(list(
    propreal = propreal,
    p = p
  ))
}


#' @noRd
#' @description Not cleaned. Spirochates, Verrucomicrobia, Cyanobacteria/Chloroplast, TM7 and pooled
ppi_microbiomedata_SVCTP <- function(){
  data("microdata", package = "scorecompdir")
  countdata=as.matrix(microdata[,12:31])

  #sample size
  n=94

  #dimension
  p=20

  #calculate totals
  tot=matrix(0,n,1)
  for (j in 1:p)
  {
    tot=tot+countdata[,j]
  }
  tot=as.vector(tot)

  #proportion data
  prop=countdata
  for (j in 1:n)
  {
    prop[j,]=countdata[j,]/tot[j]
  }

  ##Reduce dimensions to p=5


  #dimension
  p=5

  #calculate 5D dataset
  comb=matrix(0,n,p)
  comb[,1]=prop[,"Spirochaetes"]
  comb[,2]=prop[,"Verrucomicrobia"]
  comb[,3]=prop[,"Cyanobacteria/Chloroplast"]
  comb[,4]=prop[,"TM7"]
  for (j in 1:sum(p,-1))
  {
    comb[,p]=comb[,p]+comb[,j]
  }
  comb[,p]=1-comb[,p]
  colnames(comb) <- c("Spirochaetes", "Verrucomicrobia", "Cyanobacteria/Chloroplast", "TM7", "pool")

  #save data
  propreal=comb

  return(list(
    propreal = propreal,
    p = ncol(propreal)
  ))
}
