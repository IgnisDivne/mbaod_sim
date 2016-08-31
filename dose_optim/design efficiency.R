#merge results
source_dir <- "C:/Users/erist836/Desktop/MBAOD2/PKPD_ED_2in1g_misspec_odsc_cov_run_dir_2"
target_dir <- "C:/Users/erist836/Desktop/MBAOD2/PKPD_ED_2in1g_misspec_odsc_cov2_run_dir_1"
n_source<-58
n_target<-1 + 43

for(i in 1:n_source){
  load(paste(source_dir,"/results_rep_",i,".Rdata",sep=""))
  
  
  save(aod_res, file=file.path(target_dir,paste("results_rep_", n_target,".Rdata",sep="")))
  n_target <- n_target +1
}





target_dir <- "C:/Users/erist836/Desktop/MBAOD2/PKPD_ED_misspec_run_dir_1"
niter <- 50
effs <- data.frame("Iteration" = rep( -1,50*niter),"Cohort" =rep( -1,50*niter),"Efficiency" =rep( -1,50*niter))
effs2 <- effs
counter <-1
for(j in 1:niter){
  load(paste(target_dir,"/results_rep_",j,".Rdata",sep=""))
  for(i in 1:length(aod_res$stop_res)){
    
    effs[counter,]<- c(j,i,aod_res$stop_res[[paste("cohort_",i,sep="")]][[4]])
    effs2[counter,]<- c(j,i,aod_res$stop_res[[paste("cohort_",i,sep="")]][[3]])
    counter <- counter +1
  }
}
effs  <- subset(effs,effs$Iteration !=-1)
effs2 <- subset(effs2,effs2$Iteration !=-1)

cohort_eff_df <- data.frame("Cohort"=1:max(effs$Cohort),"lower" =rep(0,max(effs$Cohort)),
                            "median"=max(effs$Cohort), "upper"=rep(0,max(effs$Cohort)))


for(i in 1:max(effs$Cohort)){
  tmp <- subset(effs,effs$Cohort==i)
  cohort_eff_df[i,2:4] <- quantile(tmp[,3],probs=c(0.025,0.5,0.975))#c(min(tmp[,3]),median(tmp[,3]),max(tmp[,3]))  
}

  
library(ggplot2)
eff_plot <- ggplot(data=cohort_eff_df,aes(x=Cohort,y=median)) + geom_line() + geom_errorbar(aes(x=Cohort,ymin=lower,ymax=upper))
eff_plot
