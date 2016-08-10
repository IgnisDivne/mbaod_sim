setwd("C:/Users/erist836/Documents/GitHub/mbaod_sim/non-linear/propofol_run_run_dir_5")



reps <- 40

df_res <- data.frame(rep=1:reps, nCohorts = rep(0,reps),
                     COVs = rep(0,reps), passed=rep(0,reps), emp_dCrit = rep(0,reps))
sc_fail <- list()

df_bias <- data.frame(rep = 1:reps, TM50=rep(0,reps), hill=rep(0,reps),allCovs=rep(0,reps))


for(j in 1:reps){ 
  ncovs <- 0
  load(paste("results_rep_",j,".Rdata",sep=""))
  for(i in 1:length(aod_res$stop_res)){
    ncovs <- ncovs + !is.null(aod_res[[paste("cohort_",i,sep="")]]$est_result$cov_mat)
  }
  df_res$nCohorts[j] <- length(aod_res$stop_res)
  df_res$COVs[j] <- ncovs#/length(aod_res$stop_res)
  df_res$passed[j] <- aod_res$stop_res[[length(aod_res$stop_res)]][[1]]
  df_res$emp_dCrit[j]<- (1/det(as.matrix(aod_res[[paste("cohort_",length(aod_res)-3,sep="")]]$est_result$cov_mat[1:4,1:4]))^(1/4))
  
  
  
  df_bias[j,2:3] <- 100*(exp(c(aod_res$est_summary$THETA_3[length(aod_res$est_summary$THETA_3)],aod_res$est_summary$THETA_4[length(aod_res$est_summary$THETA_4)]))-exp(c(3.651,1.5261)))/exp(c(3.651,1.5261))
  if(length(aod_res$stop_res)==ncovs){
    df_bias$allCovs[j] <-1
  }
  
  
  if(df_res$passed[j]==F){
    sc_fail[[paste("iteration_",j,sep="")]] <- aod_res
  }
  
} 
df_bias$passed <- df_res$passed


library(reshape2)
library(ggplot2)

#plot for the number of cohorts/covs
df_res_m <- melt(df_res[,1:4], id.vars = c("rep","passed"))
df_passed <- data.frame(rep=((1:reps)*df_res$passed))
df_passed <-subset(df_passed,df_passed$rep>0)
df_passed$value <- df_res$nCohorts[df_passed$rep] +2 

df_res_passed <- subset(df_res,df_res$passed==1)

res_plot <- ggplot(data = df_res_m, aes(x=rep,y=value,fill=variable))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_grey() +theme_classic()+ geom_hline(aes(yintercept=median(df_res_passed$nCohorts), color="Median nCohorts"))+
  geom_hline(aes(yintercept = median(subset(df_res_passed,df_res_passed$nCohorts==df_res_passed$COVs)$nCohorts),color="Median nCOVS=nCohorts"))+
  geom_hline(aes(yintercept = median(subset(df_res_passed,df_res_passed$nCohorts!=df_res_passed$COVs)$nCohorts),color="Median nCOVS!=nCohorts"))+
  xlab("Simulation Nr") + ylab("Number of Cohorts/COVs")+
  geom_text(data=df_passed,label="*",aes(fill=rep("*Passed SC",length(df_passed$rep))))
res_plot




#plot for the empirical D-criterion
df_res_crit <- (df_res[,c(1,5)])
df_res_crit[df_res_crit==Inf] <-0
df_inf <-subset(df_res_crit,df_res_crit$emp_dCrit==0)
df_inf$emp_dCrit <- 2 
df_res_crit[df_res_crit==0] <-1


res_plot_crit <- ggplot(data = df_res_crit, aes(x=rep,y=emp_dCrit,fill="D-Criterion"))+
  geom_bar(stat="identity")+
  scale_fill_grey() +theme_classic()+ 
  xlab("Simulation Nr") + ylab("empirical D-criterion")+
  geom_text(data=df_inf,label="*",aes(fill=rep("* inf",length(df_inf$rep))))+scale_y_log10()
  
res_plot_crit
  

#Plot for the parameter bias

df_bias_m <- melt(df_bias[,1:4], id.vars=c("rep","allCovs"))

df_bias_p <- subset(df_bias_m,df_bias_m$allCovs==0)
bias_plot_f <- ggplot(data=df_bias_p,aes(x=variable,y=value)) +
  geom_boxplot() + ggtitle("nCovs != nCohorts")+
  xlab("Parameter") + ylab("Relative Bias [%]") +
  scale_y_continuous(breaks=seq(-100, 650, 100),limits=c(-100,650))
print("Total Relative bias (nCovs != nCohorts):")
sum(abs(df_bias_p$value))/length(df_bias_p$value)


df_bias_p <- subset(df_bias_m,df_bias_m$allCovs==1)
bias_plot_s <- ggplot(data=df_bias_p,aes(x=variable,y=value)) +
  geom_boxplot() + ggtitle("nCovs == nCohorts")+
  xlab("Parameter") + ylab("Relative Bias [%]")+
  scale_y_continuous(breaks=seq(-100, 650, 100),limits=c(-100,650))
print("Total Relative bias (nCovs == nCohorts):")
sum(abs(df_bias_p$value))/length(df_bias_p$value)


df_bias_p <- df_bias_m
bias_plot_all <- ggplot(data=df_bias_p,aes(x=variable,y=value)) +
  geom_boxplot() + ggtitle("All")+
  xlab("Parameter") + ylab("Relative Bias [%]")+
  scale_y_continuous(breaks=seq(-100, 650, 100),limits=c(-100,650))
print("Total Relative bias (All):")
sum(abs(df_bias_p$value))/length(df_bias_p$value)


multiplot(bias_plot_all,bias_plot_f,bias_plot_s,cols=3)





df_bias_m <- melt(df_bias[-df_inf$rep,1:4], id.vars=c("rep","allCovs"))

df_bias_p <- subset(df_bias_m,df_bias_m$allCovs==0)
bias_plot_f <- ggplot(data=df_bias_p,aes(x=variable,y=value)) +
  geom_boxplot() + ggtitle("nCovs != nCohorts")+
  xlab("Parameter") + ylab("Relative Bias [%]") +
  scale_y_continuous(breaks=seq(-100, 650, 100),limits=c(-100,650))
print("Total Relative bias (nCovs != nCohorts):")
sum(abs(df_bias_p$value))/length(df_bias_p$value)


df_bias_p <- subset(df_bias_m,df_bias_m$allCovs==1)
bias_plot_s <- ggplot(data=df_bias_p,aes(x=variable,y=value)) +
  geom_boxplot() + ggtitle("nCovs == nCohorts")+
  xlab("Parameter") + ylab("Relative Bias [%]")+
  scale_y_continuous(breaks=seq(-100, 650, 100),limits=c(-100,650))
print("Total Relative bias (nCovs == nCohorts):")
sum(abs(df_bias_p$value))/length(df_bias_p$value)


df_bias_p <- df_bias_m
bias_plot_all <- ggplot(data=df_bias_p,aes(x=variable,y=value)) +
  geom_boxplot() + ggtitle("All")+
  xlab("Parameter") + ylab("Relative Bias [%]")+
  scale_y_continuous(breaks=seq(-100, 650, 100),limits=c(-100,650))
print("Total Relative bias (All):")
sum(abs(df_bias_p$value))/length(df_bias_p$value)


multiplot(bias_plot_all,bias_plot_f,bias_plot_s,cols=3)


#Get the empirical d-criterion for the failed runs
for(i in 1:length(sc_fail)){
  print(1/det(as.matrix(sc_fail[[i]][[paste("cohort_",length(sc_fail[[i]])-3,sep="")]]$est_result$cov_mat[1:4,1:4]))^(1/4))
  
}

sc_fail$iteration_35$cohort_50$opt_result$opt_output$poped.db$design$x