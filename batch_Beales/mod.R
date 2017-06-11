#'@title mod
#'@description Function computes the modulus as computed in SAS.
#'@param x numeric vector on which to compute the modulus
#'@param m divisor
#'@return vector of moduli
#'@example
#'mod(c(-7,14,20,-81,23),4)


mod<-function(x,m){
  remain<-sapply(x,function(x) (x/m-trunc(x/m))*m)
  return(remain) 
}