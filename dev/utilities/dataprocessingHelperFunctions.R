#This file contains helper functions related to data processing

#' Calculating Z-scores for a matrix or data frame
#'
#' @description This function calculates a z-score (for each row by default).  
#' This normalises the rows so they have a mean of zero and a standard deviation of 1,
#' which allows for easier comparisons between 
#'
#' @param mat This is a numerical matrix or data frame.
#' @param byRow This specifies whether the z-scores will be calculated for every row or every column in the matrix or data frame
#'
#' @return a matrix in which every row (or column, if byRow=FALSE) has a mean of 0 and a standard deviation of 1.
#'
#' @examples mat=matrix(data=c(1, 3, 4, 5, 2, 4, 5, 2), byrow = T, nrow = 2)
#' zscore(mat)
#' zscore(mat, byRow=FALSE)
#' 
zscore<-function(mat, byRow=TRUE){
  side=1
  if(!byRow){side=2}
  
  temp=apply(mat, side, function(i){
    (i-mean(i))/sd(i)
  })
  
  if(dim(temp)[1]!=dim(mat)[1] & dim(temp)[2]!=dim(mat)[2]){
    t(temp)
  }else{temp}
  
}

#' Convert output of GENIE3 into an adjacency matrix
#'
#' @description This takes a matrix where rows are regulators and columns are target and values are scores
#' It converts the table to a three column layout where the score of each regulator-target pair is specified
#' 
#' @param network This is the network (usually the output of GENIE3)
#' @param threshold This is the threshold of score above which the edge will be included in the final table
#'
#' @return a three column table, containing the regulator, the target, and the score.  Only rows with scores exceeding the threshold will be provided.
#'
convertToAdjacency <- function(network, threshold){

  temp=unlist(sapply(rownames(network), function(i){
    paste(i, colnames(network)[which(network[i,]>threshold)], network[i,which(network[i,]>threshold)])
  }))
  
  col3=data.frame(t(sapply(temp, function(i){
    split=strsplit(i, ' ')
    c(split[[1]][1], split[[1]][2], split[[1]][3])
  })))
  
  col3[,3]=as.numeric(col3[,3])
  
  col3
}



#' Calculating pafway network
#'
#' @description This function finds associations between terms in the descriptions of nodes in a network
#'
#' @param GO a vector of descriptions of nodes in the network.  The vector needs names that match the nodes in the network, but can contain extra genes.
#' @param edges the network, as a matrix or data frame.  Column 1 is Source nodes and Column 2 is the targets.
#' @param GOtypes  a vector of terms that are of interest
#' @return a matrix in which each value represents the p-value of the term corresponding to the COLUMN is upstream of the term corresponding to the ROW.
#'
#' 
pafway <- function(GO, edges, GOtypes) {
  
  GOinNetwork = GO[unique(c(edges[, 1], edges[, 2]))]
  
  grepLen=sapply(GOtypes, function(i){
    length(grep(i, GOinNetwork))
  })
  names(grepLen)=GOtypes
  
  sapply(GOtypes, function(i) {
    if(grepLen[i]!=0){
      #find edges first:
      grepI=grep(i, GOinNetwork[edges[, 1]])
      
      sapply(GOtypes, function(j) {
        if(grepLen[j]!=0){
          
          a=length(grep(j, GOinNetwork[edges[grepI, 2]]))
          
          p_bot = grepLen[i]/length(GOinNetwork) * grepLen[j]/length(GOinNetwork)
          
          b = stats::binom.test(a, length(edges[, 1]), p = p_bot, alternative = c("greater"))
          b$p.value
        }else{1}
      })
    }else{rep(1, length(GOtypes))}
  })
}



#' Title
#'
#' @param i 
#'
#' @return
#' @export
#'
#' @examples
hello <-function(i){
  
  
  
  
}