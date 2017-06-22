#'@title collapse_stratbins
#'@description Collapse the smallest stratbin given by stratbins, where the collapsed stratbin is merged 
#'with the smallest neighboring stratbin. Module modifies stratbins, stratbin, and strat_nsamp in parent.frame() 
#'@param stratbin stratbin-identifier vector 
#'@param stratbins vector indicating in stratbin
#'@param strat_nsamp number of samples in each stratbin given by `nsamp(stratbin,ifhi,ifsamp,stratbins,modulus)`
#'@keywords internal

collapse_stratbins<-function(stratbin,strat_nsamp,stratbins){
  
  if (length(stratbins)>1){
    j_collapse<-min(which(strat_nsamp==min(strat_nsamp))[1])
    j_merge<-altmod(c(j_collapse-1,j_collapse+1),length(strat_nsamp))
    j_merge<-j_merge[which(strat_nsamp[j_merge]==min(strat_nsamp[j_merge]))[1]]
    
    loc_collapse_stratbin<-which(stratbin==stratbins[j_collapse])
    
    stratbin[loc_collapse_stratbin]<- stratbins[j_merge]
    stratbins <- stratbins[c(seq(1:length(stratbins))!=j_collapse)]
    
    strat_nsamp[j_merge]<- sum(strat_nsamp[c(j_merge,j_collapse)])
    strat_nsamp<- strat_nsamp[c(seq(1:length(strat_nsamp))!=j_collapse)]
    
  }
  
  assign("stratbins",stratbins,parent.frame())
  assign("strat_nsamp",strat_nsamp,parent.frame())
  assign("stratbin",stratbin,parent.frame())
  
}

