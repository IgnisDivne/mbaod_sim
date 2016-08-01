param_sim_mat <- function(CL=1.0,V=3,omCL=0.05,omV=0.05,TM50=75,PMA,WT,p=F){
  if(p==F){
    LCL <-rnorm(length(PMA),CL,sqrt(omCL))
    LV <- rnorm(length(PMA),V,sqrt(omV))
  }else{
    LCL <-qnorm(p,CL,sqrt(omCL))
    LV <- qnorm(p,V,sqrt(omV))
    
  }
  LCLi=((PMA)/(PMA + TM50))
  LVi=LV+log(WT/70)
  return(list(LCLi,LVi))
  
}
#Function to simulate LCL and LV parameters
param_sim <- function(CL=1.0,V=3,sdCL=0.05,sdV=0.05,TM50=3.651,HILL=1.5261,PMA,WT,p=F){
  if(p==F){
    LCL <-rnorm(length(PMA),CL,sqrt(sdCL))
    LV <- rnorm(length(PMA),V,sqrt(sdV))
  }else{
    LCL <-qnorm(p,CL,sqrt(sdCL))
    LV <- qnorm(p,V,sqrt(sdV))
    
  }
  HILL=exp(HILL)
  TM50=exp(TM50)
  LCLi=LCL+0.75*log(WT/70)+log((PMA^HILL)/(PMA^HILL + TM50^HILL))
  LVi=LV+log(WT/70)
  return(list(LCLi,LVi))
  
}
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
  n_children <- c()
  age_group <- c()
  nID_group <-data.frame(Group_1=rep(0,100),Group_2=rep(0,100),Group_3=rep(0,100),
                         Group_4=rep(0,100),Group_5=rep(0,100),Group_6=rep(0,100))
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
    groups <- results_all[[paste("iteration_",i,sep="")]]$final_design$x[2:length(results_all[[paste("iteration_",i,sep="")]]$final_design$x)]
    nID_tot[i] <- sum(results_all[[paste("iteration_",i,sep="")]]$final_design$groupsize)-100
    nID <- results_all[[paste("iteration_",i,sep="")]]$final_design$groupsize[2:length(results_all[[paste("iteration_",i,sep="")]]$final_design$groupsize)]
    
    rep_nr <- c(rep_nr,(rep(i,length(groups))))
    n_children <- c(n_children,cumsum(nID))
    age_group <- c(age_group,groups)

  }
  df_nID <- data.frame(Iteration = rep_nr, n_children,age_group)
  return(list(df_nID,nID_tot))
}