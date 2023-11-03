#This file contains helper functions related to data processing

#' Calculating Z-scores for a matrix or data frame
#'
#' @param mat This is a numerical matrix or data frame.
#' @param byRow This specifies whether the z-scores will be calculated for every row or every column in the matrix or data frame
#'
#' @return a matrix in which every row (or column, if byRow=FALSE) has a mean of 0 and a standard deviation of 1.
#' @export
#'
#' @examples  
#' 
zscore<-function(mat, byRow=TRUE){
  side=1
  if(!byRow){side=2}
  
  apply(mat, side, function(i){
    (i-mean(i))/sd(i)
  })
  
}