#'@title nsamp
#'@description Function computes the number of samples associated with stratbins. 
#'Requires previous definition of the stratbin-identifier vector,
#'called stratbin, and the indicies of the elements of stratbin corresponding
#'to observations that are sampled, called locsamp. 
#'@param stratbin stratbin-identifier vector 
#'@param ifhi vector of 0's and 1's, 1 indicates hiflow
#'@param ifsamp vector of 0's and 1's, 1 indicates sampled
#'@param stratbins vector indicating in stratbin
#'@param modulus 0 or 1 indicating hiflow or lowflow
#'@return numeric vector
#'@keywords internal

nsamp<-function(stratbin,ifhi,ifsamp,stratbins,modulus){
  instratbins<-data.frame(stratbin,ifhi,ifsamp)
  instratbins<-instratbins[which(mod(instratbins$stratbin,2)==modulus & instratbins$ifsamp==1),]
  instratbins<-aggregate(instratbins[c("stratbin")],by=list(instratbins$stratbin),length)
  
  result<-numeric(length(stratbins))
  result[which(stratbins %in% instratbins$Group.1)]<-instratbins$stratbin
  return(result)
}
