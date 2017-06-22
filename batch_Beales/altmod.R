#'@title altmod
#'@description Function computes the alternative modulus that returns m when mod(x,m) = 0, 
#'for 0 <= x <= m, otherwise returns mod(x,m). Converted from FLUXMASTER in SAS
#'@param x numeric vector on which to compute the alternative modulus
#'@param m number by which to compute the alternative modulus (divisor)
#'@return vector of alternative moduli
#'@examples
#'altmod(c(-7,14,20,-81,23),4)


altmod<-function(x,m){
  return(m-mod(c(2*m-x),m))
}

