setwd("C:/Users/erist836/Documents/GitHub/mbaod_sim/non-linear/propofol_run_run_dir_1")

df_res <- data.frame(rep=1:16, nCohorts = rep(0,16),
                     COVs = rep(0,16), passed=rep(0,16))
for(j in 1:16){
  ncovs <- 0
  load(paste("results_rep_",j,".Rdata",sep=""))
  for(i in 1:length(aod_res$stop_res)){
    ncovs <- ncovs + !is.null(aod_res[[paste("cohort_",i,sep="")]]$est_result$cov_mat)
  }
  df_res$nCohorts[j] <- length(aod_res$stop_res)
  df_res$COVs[j] <- ncovs#/length(aod_res$stop_res)
  df_res$passed[j] <- aod_res$stop_res[[length(aod_res$stop_res)]][[1]]
  
}

library(reshape2)
library(ggplot2)


df_res_m <- melt(df_res, id.vars = c("rep","passed"))
df_passed <- data.frame(rep=((1:16)*df_res$passed))
df_passed <-subset(df_passed,df_passed$rep>0)
df_passed$value <- df_res$nCohorts[df_passed$rep] +2 

res_plot <- ggplot(data = df_res_m, aes(x=rep,y=value,fill=variable))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_grey() +theme_classic()+ 
  xlab("Simulation Nr") + ylab("Number of Cohorts/COVs")+
  geom_text(data=df_passed,label="*",aes(fill=rep("*Passed SC",length(df_passed$rep))))
res_plot

  