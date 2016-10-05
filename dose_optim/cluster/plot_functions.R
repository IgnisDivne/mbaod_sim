
##Function to get power
get_power <- function(sdlogX){
  f<-function(s,n,sigma){
    2*((n-1)/2/sigma**2)**((n-1)/2)/gamma(0.5*(n-1))*s**(n-2)*exp(-(n-1)*s**2/2/sigma**2)}
  nmin<-4
  nmax<-25
  nsub<-nmax-nmin+1
  result<-rep(0,nsub)
  sd<-sdlogX #standard deviation of logX from adult data
  for (i in 1:nsub) {
    n<-nmin+i-1
    tv<-qt(0.975, n-1)
    sup<-log(1.4)*sqrt(n)/tv
    g<-function(x){f(x,n,sd)}
    poweri<-integrate(g, 0, sup, subdivisions=100)
    result[i]<-poweri[[1]]
  }
  power<-data.frame(nsub=c(nmin:nmax), power=result)
  return(power)
}


###Function to get the total number of children and children per cohort
get_nID <- function(res_path,n_iter,full_path=T){
  rep_nr <- c()
  n_ID <- c()
  doses_tot <- c()
  nID_tot <- 1:length(n_iter)
  
  if(full_path==F){
    results_all <- list()
    file_dir <- res_path
    for(i in 1:n_iter){
      
      #data frame with the cumulative number of individuals per group per iteration
      load(paste(file_dir,"/results_rep_",i,".Rdata",sep=""))
      results_all[[paste("iteration_",i,sep="")]] <- aod_res
    }
  }else{    
    load(res_path)
  }
  for(i in 1:n_iter){
    doses <- results_all[[paste("iteration_",i,sep="")]]$final_design$a
    nID_tot[i] <- sum(results_all[[paste("iteration_",i,sep="")]]$final_design$groupsize)
    nID <- results_all[[paste("iteration_",i,sep="")]]$final_design$groupsize
    
    rep_nr <- c(rep_nr,(rep(i,length(doses))))
    n_ID <- c(n_ID,nID)
    doses_tot <- c(doses_tot,doses)

  }
  df_nID <- data.frame(Iteration = rep_nr, n_ID,doses_tot)
  return(list(df_nID,nID_tot))
}