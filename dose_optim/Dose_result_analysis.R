
#set path to the directory where this file is located
setwd("C:/Users/erist836/Desktop/MBAOD2")

# remove things from the global environment
rm(list=ls())

library(PopED)
library(mvtnorm)
library(xpose4)

# load the MBAOD package change to your relevant load call
devtools::load_all("C:/Users/erist836/Documents/GitHub/MBAOD/R")
# load the weight estimator file
source("stop_crit_dose2.R")
source("plot_functions.R")



#mbaod d
res_dir <- ("C:/Users/erist836/Desktop/MBAOD2/PKPD_D_2in1g_misspec_run_run_dir_1")
df_d <- unlist(get_nID(res_dir,n_iter=100,full_path = F)[2])


#mbaod ed
res_dir <- "C:/Users/erist836/Desktop/MBAOD2/PKPD_ED_misspec_run_dir_1"
df_ed <- unlist(get_nID(res_dir,n_iter=100,full_path = F)[2])


#####Restricted OD#########

#ODsc D
file_dir <- c("C:/Users/erist836/Desktop/MBAOD2/PKPD_D_2in1g_misspec_odsc_cov_run_dir_3")
df_OD_d <- unlist(get_nID(res_dir,n_iter=100,full_path = F)[2])



#ODsc ED
file_dir <- "C:/Users/erist836/Desktop/MBAOD2/PKPD_ED_2in1g_misspec_odsc_cov2_run_dir_1"
df_OD_ed <- unlist(get_nID(file_dir,n_iter=100,full_path=F)[2])


##Create the data frames
df_box_d <- data.frame("D MBAOD"=df_d,"ED MBAOD"=df_ed,"D ODsc"=df_OD_d,"ED ODsc"=df_OD_ed)
                           


library(reshape2)
df_box_d <- melt(df_box_d)

library(ggplot2)

##Plot the total number of IDS
p_box_true <- ggplot(data=df_box_d,aes(x=variable,y=value)) + geom_boxplot() +
  xlab(" ") + ylab("Number of Individuals")+ ggtitle("")+
  scale_x_discrete(labels = c("MBAOD lnD","MBAOD ELD","ODsc lnD","ODsc ELD"))+scale_y_continuous(breaks=seq(0, 45, 5),limits=c(0,45))+
  theme_classic()
p_box_true



#############################################
#########MBAOD and OD Group Selection Plots#####
############################################


eval_design <- function(result_dir=NULL,niterations=100,old=T){
  df_params <- data.frame(iteration=1:100, base=rep(0,100),EMAX=rep(0,100),EC50=rep(0,100),HILL=rep(0,100))
  power_res <- 0
  not_there <- 0
  if(is.null(result_dir)) {
    print("Please give a results directory")
  }else{
    for(i in 1:niterations){
      load(paste(result_dir,"/results_rep_",i,".Rdata",sep=""))
      
      
      df_params$base[i] <- aod_res$est_summary$THETA_4[length(aod_res$est_summary$THETA_4)]
      df_params$EMAX[i] <- aod_res$est_summary$THETA_5[length(aod_res$est_summary$THETA_5)]
      df_params$EC50[i] <- aod_res$est_summary$THETA_6[length(aod_res$est_summary$THETA_6)]
      df_params$HILL[i] <- aod_res$est_summary$THETA_7[length(aod_res$est_summary$THETA_7)]
      
      
      if(length(aod_res$stop_res)==0){
        not_there <- not_there +1
        
      }else{
        if(old==T){
          power_res <- power_res + sum(aod_res$stop_res$cohort_1[[1]])
        }else{
          
          stopres<- aod_res$stop_res[[length(aod_res$stop_res)]]
          power_res <- power_res + sum(stopres[[1]]) 
        }
        
      }
      
      
    }
  }
  power_res <- 100*power_res/(100-not_there)
  print(paste(i-not_there, "of",i,"stop_res results found."))
  
  return(list(df_params,power_res))
} 

for(i in 1:100){
  load(paste(result_dir,"/results_rep_",i,".Rdata",sep=""))
  print(length(aod_res$stop_res))
}


df_LLT <- data.frame(Age_group_MBAOD = rep(1:7),
                     nID_MBAOD=rep(0,7),
                     nID_ODSC=rep(0,7),
                     nID_Scaling=rep(0,7),
                     nID_Adult_Prior=rep(0,7),
                     nID_OD=rep(0,7))
df_LLM <- df_LLT

#No misspecification values


result_dir <- ("C:/Users/erist836/Desktop/MBAOD2/PKPD_D_2in1g_misspec_run_run_dir_1")
MBAOD_d <- eval_design(result_dir = result_dir,old = F) 



result_dir <- "C:/Users/erist836/Desktop/MBAOD2/PKPD_ED_misspec_run_dir_1"
MBAOD_ed <- eval_design(result_dir = result_dir,old = F) 

#ODsc 

result_dir <- c("C:/Users/erist836/Desktop/MBAOD2/PKPD_D_2in1g_misspec_odsc_cov_run_dir_3")
OD_d <- eval_design(result_dir = result_dir,old = F) 


result_dir <- "C:/Users/erist836/Desktop/MBAOD2/PKPD_ED_2in1g_misspec_odsc_cov2_run_dir_1"
OD_ed <- eval_design(result_dir = result_dir,old = F) 


bias_df_true <- data.frame(DA = rep(c("MBAOD lnD", "MBAOD ELD", "ODsc lnD","ODsc ELD"), 1, each=100),
                           base = rep(0,400),
                           EMAX  = rep(0,400),
                           EC50 = rep(0,400),
                           Hill = rep(0,400)) 

bias_df_true$base <- 100*((c(MBAOD_d[[1]]$base,MBAOD_ed[[1]]$base, OD_d[[1]]$base , OD_ed[[1]]$base ))-1)/1

bias_df_true$EMAX <- 100*((c(MBAOD_d[[1]]$EMAX,MBAOD_ed[[1]]$EMAX, OD_d[[1]]$EMAX,OD_ed[[1]]$EMAX))-100)/100

bias_df_true$EC50 <- 100*((c(MBAOD_d[[1]]$EC50,MBAOD_ed[[1]]$EC50, OD_d[[1]]$EC50,OD_ed[[1]]$EC50 ))-7)/7

bias_df_true$Hill <- 100*((c(MBAOD_d[[1]]$HILL,MBAOD_ed[[1]]$HILL, OD_d[[1]]$HILL,OD_ed[[1]]$HILL )) -2)/2
library(reshape) 
library(ggplot2)


i=4
titles<- c("MBAOD lnD", "MBAOD ELD", "ODsc lnD","ODsc ELD")  #400X400PX
df_box_tmp <- data.frame(Base = subset(bias_df_true,bias_df_true$DA==titles[i])$base,
                         EMAX = subset(bias_df_true,bias_df_true$DA==titles[i])$EMAX,
                         EC50 = subset(bias_df_true,bias_df_true$DA==titles[i])$EC50,
                         Gamma = subset(bias_df_true,bias_df_true$DA==titles[i])$Hill) 

df_box_m <- melt(df_box_tmp)

p_biasT <- ggplot(data=df_box_m,aes(x=variable,y=value)) + geom_hline(data=data.frame(rep(0,4)))+ 
  xlab("Parameter") + ylab("% BIAS")+ ggtitle(titles[i])+
  scale_y_continuous(breaks=seq(-50, 50, 10),limits=c(-50,50)) + geom_boxplot() +
  theme_classic() 
p_biasT


library(scales)

res_dir <- ("C:/Users/erist836/Desktop/MBAOD2/PKPD_D_2in1g_misspec_run_run_dir_1")
d_mbaod_nid <- get_nID(res_dir,n_iter=100,full_path = F)[[1]]
d_mbaod_p <- ggplot(data=d_mbaod_nid,aes(x=doses_tot,y=..count../sum(..count..))) + geom_histogram(binwidth=5) +theme_classic() + 
  ggtitle("MBAOD lnD") + xlab("Dose") + ylab("Relative Number of Cohorts [%]") +
  scale_y_continuous(limits = c(0,0.26),breaks = seq(0,0.25,0.05), labels = seq(0,25,5))
d_mbaod_p 

res_dir <- "C:/Users/erist836/Desktop/MBAOD2/PKPD_ED_misspec_run_dir_1"
ed_mbaod_nid <- get_nID(res_dir,n_iter=100,full_path = F)[[1]]
ed_mbaod_p <- ggplot(data=ed_mbaod_nid,aes(x=doses_tot,y=..count../sum(..count..))) + geom_histogram(binwidth=5) +theme_classic() + 
  ggtitle("MBAOD ELD") + xlab("Dose") + ylab("Relative Number of Cohorts [%]") +
  scale_y_continuous(limits = c(0,0.26),breaks = seq(0,0.25,0.05), labels = seq(0,25,5))#+  scale_y_continuous(labels = percent_format())
ed_mbaod_p 

res_dir <- c("C:/Users/erist836/Desktop/MBAOD2/PKPD_D_2in1g_misspec_odsc_cov_run_dir_3")
d_odsc_nid <- get_nID(res_dir,n_iter=100,full_path = F)[[1]]
d_odsc_p <- ggplot(data=d_odsc_nid,aes(x=doses_tot,y=..count../sum(..count..))) + geom_histogram(binwidth=5) +theme_classic() + 
  ggtitle("ODsc lnD") + xlab("Dose") + ylab("Relative Number of Cohorts [%]") +
  scale_y_continuous(limits = c(0,0.26),breaks = seq(0,0.25,0.05), labels = seq(0,25,5)) + #+  scale_y_continuous(labels = percent_format())
  scale_x_continuous(limits=c(0,500))
d_odsc_p 

res_dir <- c("C:/Users/erist836/Desktop/MBAOD2/PKPD_ED_2in1g_misspec_odsc_cov2_run_dir_1")
ed_odsc_nid <- get_nID(res_dir,n_iter=100,full_path = F)[[1]]
ed_odsc_p <- ggplot(data=ed_odsc_nid,aes(x=doses_tot,y=..count../sum(..count..))) + geom_histogram(binwidth=5) +theme_classic() + 
  ggtitle("ODsc ELD") + xlab("Dose") + ylab("Relative Number of Cohorts [%]") +
  scale_y_continuous(limits = c(0,0.26),breaks = seq(0,0.25,0.05), labels = seq(0,25,5)) + #+  scale_y_continuous(labels = percent_format())
  scale_x_continuous(limits=c(0,500))
ed_odsc_p 






target_dir <- c("C:/Users/erist836/Desktop/MBAOD2/PKPD_D_2in1g_misspec_run_run_dir_1",
                "C:/Users/erist836/Desktop/MBAOD2/PKPD_ED_misspec_run_dir_1",
              "C:/Users/erist836/Desktop/MBAOD2/PKPD_D_2in1g_misspec_odsc_cov_run_dir_3",
              "C:/Users/erist836/Desktop/MBAOD2/PKPD_ED_2in1g_misspec_odsc_cov2_run_dir_1")


titles<- c("MBAOD lnD", "MBAOD ELD", "ODsc lnD", "ODsc ELD")  



t=4
niter <- 50
effs <- data.frame("Iteration" = rep( -1,50*niter),"Cohort" =rep( -1,50*niter),"Efficiency" =rep( -1,50*niter))
effs2 <- effs
counter <-1
for(j in 1:niter){
  load(paste(target_dir[t],"/results_rep_",j,".Rdata",sep=""))
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
eff_plot <- ggplot(data=cohort_eff_df,aes(x=Cohort,y=median)) + geom_line()+ 
geom_errorbar(aes(x=Cohort,ymin=lower,ymax=upper), color = "grey", width=0.3) + theme_classic() +
  ylab("Relative efficiency") + ggtitle(titles[t]) + scale_x_continuous(limits=c(0,18))+
  scale_y_continuous(limits=c(0.5,1.1))
eff_plot




























##Coverage tests of CI in stopping criteria##

##The code will for one misspecification of TM50,
#extract the final values of CL and V and 


##Set path to directory of this file
setwd("C:/Users/erist836/Desktop/MBAOD_project")

# remove things from the global environment
rm(list=ls())

library(PopED)
library(mvtnorm)
library(xpose4)
library(reshape2)

# load the MBAOD package (CHange to your path)
devtools::load_all("C:/Users/erist836/Desktop/WarfarinPKPD/R/MBAOD")

# load the PopED model file
source("PopED_files/poped.mod.PK.1.comp.maturation_Xlog.R")
# load the weight estimator file
source("get_weight.R")
source("stop_critX.R")
source("plot_functions.R")
corr<-0




#MBAOD Propofol True
res_dir <-c("C:/Users/erist836/Desktop/MBAOD_project/merge/")

#MBAOD Propofol GFR
res_dir<-("C:/Users/erist836/Desktop/MBAOD_project/merge_gfr/")





##For the propofol examples
size_mat_scaling_p <- function(params,age,wt){
  
  base<-exp(params[,1])
  TM50<- exp(params[,3])
  HILL<-exp(params[,4])
  V <- exp(params[,2])
  
  
  
  CL <-  base * (wt/70)^0.75 * (age^HILL)/(age^HILL+TM50^HILL)
  V  <-  V *(wt/70)
  
  df <- data.frame(CL,V)
  return(df)
}

group_distribution <-data.frame(age_group=1:6,
                                nID=rep(0,6))
# alpha <- 0.05
# results_i1c3 <- results_all$iteration_1$cohort_3$est_result$thetas
# results_i4c3 <- results_all$iteration_4$cohort_3$est_result$thetas
# CL_vec <-1:100
# V_vec <- 1:100
# tm50_vec <- 1:100
# HILL_vec <- 1:100


WT <- c(6.906,8.4,9.29,11.11,13.72,15.96,18.14,21.15,23.96,26.74,31.59,36.03,40.53,47.06,51.908,58.02,62.78,66.96,69.11,71.48,73.49)
AGE <- get_age_from_median_weight(WT)
AGE_WT <- data.frame(AGE,WT)




if(max(as.numeric(run_summary$cohort))>1) run_summary <- subset(run_summary, run_summary$cohort!=1) 

final_estimates <- data.frame(iter=1:100, CL=rep(0,100),V=rep(0,100),TM50=rep(0,100),HILL=rep(0,100),cohort=rep(0,100))


d9_CI_CL <- data.frame(iter=rep(1:100,rep(21,100)),pma=rep(AGE$PMA,100),d9_lower=rep(0,2100),d9_higher=rep(0,2100),
                       np=rep(0,2100))
d9_CI_V <-d9_CI_CL 


for( i in 1:100){
  load(paste(res_dir,"results_rep_",i,".Rdata",sep=""))
  tmp <- data.frame(aod_res$est_summary)
  
  
  final_estimates[i,2:5] <-tmp[length(tmp[,1]),2:5] 
  final_estimates[i,6] <-max(as.numeric(tmp$cohort))
  
  iteration <- aod_res
  
  x  <- iteration$final_design$x[-1]
  gs <- iteration$final_design$groupsize[-1]
  
  for(j in 1:length(x)){
    group_distribution[x[j],2] <- group_distribution[x[j],2] + gs[j]
    
  }
  
  if(length(iteration$stop_res)>0){
    stop_res<-iteration$stop_res[[length(iteration$stop_res)]]
    
    final_cohort_crit_res_CL <- stop_res[[3]]
    final_cohort_crit_res_V <- stop_res[[4]]
  }else{
    final_cohort_crit_res_CL <- matrix(nrow=21,ncol=4)
    final_cohort_crit_res_V <- matrix(nrow=21,ncol=4)
  }
  for(k in 1:21){
    #     
    age = AGE_WT$PMA[k]
    wt = AGE_WT$WT[k]
    
    d9_CI_CL[k+21*(i-1),5] <- size_mat_scaling_p(final_estimates[i,2:5],age,wt)[1]
    d9_CI_V[k+21*(i-1),5]  <- size_mat_scaling_p(final_estimates[i,2:5],age,wt)[2]
    #     
    #  
    
    d9_CI_CL[k+21*(i-1),3:4]= final_cohort_crit_res_CL[k,3:4] #c(exp(-qt((1-alpha/2),df)*SElcl),exp(qt((1-alpha/2),df)*SElcl))
    d9_CI_V[k +21*(i-1),3:4]= final_cohort_crit_res_V[k,3:4] #c(exp(-qt((1-alpha/2),df)*SE_V),exp(qt((1-alpha/2),df)*SE_V))
    
  }
  print(i)
}


df_CL_coverage <- data.frame(pma=AGE_WT$PMA,d9_ll=rep(0,21),d9_lh=rep(0,21),
                             d9_hl=rep(0,21),d9_hh=rep(0,21),
                             np_lower=rep(0,21),np_higher=rep(0,21))

df_V_coverage <- df_CL_coverage

for(i in 1:21){
  tmp2 <- subset(d9_CI_CL,pma==AGE_WT$PMA[i])
  df_CL_coverage[i,2:3] <- quantile((tmp2$d9_lower),probs=c(0.025,0.975),na.rm=T)#c(min(tmp2$d9_lower,na.rm=T),max(tmp2$d9_lower,na.rm=T))
  df_CL_coverage[i,4:5] <- quantile((tmp2$d9_higher),probs=c(0.025,0.975),na.rm=T)#c(min(tmp2$d9_higher,na.rm=T),max(tmp2$d9_higher,na.rm=T))
  df_CL_coverage[i,6:7] <- quantile((tmp2$np),probs=c(0.025,0.975),na.rm=T)/gm_mean((tmp2$np))
  
  
  tmp3 <- subset(d9_CI_V,pma==AGE_WT$PMA[i])
  df_V_coverage[i,2:3] <- quantile((tmp3$d9_lower),probs=c(0.025,0.975),na.rm=T)
  df_V_coverage[i,4:5] <- quantile((tmp3$d9_higher),probs=c(0.025,0.975),na.rm=T)
  df_V_coverage[i,6:7] <- quantile((tmp3$np),probs=c(0.025,0.975),na.rm=T)/gm_mean((tmp3$np))
}

df_CL_coverage_m<-data.frame(pma=rep(AGE_WT$PMA,3),Legend=c(rep("Precicted CI Limits",42),rep("CI of Estimates",21)),
                             lower=c(df_CL_coverage[,2],df_CL_coverage[,4],df_CL_coverage[,6]),
                             higher=c(df_CL_coverage[,3],df_CL_coverage[,5],df_CL_coverage[,7]))

df_V_coverage_m<-data.frame(pma=rep(AGE_WT$PMA,3),Legend=c(rep("Predicted CI Limits",42),rep("CI of Estimates",21)),
                            lower=c(df_V_coverage[,2],df_V_coverage[,4],df_V_coverage[,6]),
                            higher=c(df_V_coverage[,3],df_V_coverage[,5],df_V_coverage[,7]))


x_step=c(rep(seq(0.8,20.8,1),2),seq(1.2,21.2,1))

#950x400
coverage_CI_CL <- ggplot(data=(df_CL_coverage_m),aes(x=x_step,y=lower)) +
  geom_errorbar(aes(ymin=lower*100, ymax=higher*100,colour=Legend)) +
  xlab("PMA") + ylab("% of geometric mean")+ggtitle("Propofol GFR: CI coverage of CLp")+
  theme(axis.text = element_text(size=20,colour="black"),axis.title = element_text(size=23,face="bold"),
        axis.text.x = element_text(angle = 60,vjust = 0.65),plot.title=element_text(size=23,face="bold"))+
  scale_x_discrete(breaks=1:21,labels=(AGE_WT$PMA))+scale_y_continuous(breaks=seq(50,150,10),limits=c(50,150))
coverage_CI_CL

coverage_CI_V <- ggplot(data=(df_V_coverage_m),aes(x=x_step,y=lower)) +
  geom_errorbar(aes(ymin=lower*100, ymax=higher*100,colour=Legend)) +
  xlab("PMA") + ylab("% of geometric mean")+ggtitle("Propofol GFR: CI coverage of Vp")+
  theme(axis.text = element_text(size=20,colour="black"),axis.title = element_text(size=23,face="bold"),
        axis.text.x = element_text(angle = 60,vjust = 0.65),plot.title=element_text(size=23,face="bold"))+
  scale_x_discrete(breaks=1:21,labels=(AGE_WT$PMA))+scale_y_continuous(breaks=seq(90,110,2),limits=c(90,110))
coverage_CI_V










rm(list=ls())

#MBAOD LLT
res_dir <-("Page_results/bridging_maturation_model_xoptim_restricted_group_true_run_dir_1/")

#MBAOD LLM
res_dir <-("Page_results/bridging_maturation_model_xoptim_restricted_group_large2_run_dir_1/")   

res_dir <- "C:/Users/erist836/Documents/GitHub/mbaod_sim/bridging_maturation_model_xoptim_restricted_group_large_run_dir_1/"
library(reshape)
library(ggplot2)
source("get_weight.R")
source("stop_critX.R")


##For the LL examples
size_mat_scaling_l <- function(params,age,wt){
  
  base<-exp(params[,1])
  TM50<- (params[,3])
  HILL<-1
  V <- exp(params[,2])
  
  
  
  CL <-  base * (wt/70)^0.75 * (age^HILL)/(age^HILL+TM50^HILL)
  V  <-  V *(wt/70)
  
  df <- data.frame(CL,V)
  return(df)
}

size_mat_scaling_log <- function(params,age,wt){ 
  
  base<-params[,1]
  TM50<- params[,3]
  HILL<-1
  V <- params[,2]
  
  
  
  CL <-  base + log((wt/70))*0.75 + log((age^HILL)/(age^HILL+TM50^HILL))
  V  <-  V +(wt/70) 
  
  df <- data.frame(CL,V) 
  return(df)
}
group_distribution <-data.frame(age_group=1:6,
                                nID=rep(0,6))
# alpha <- 0.05
# results_i1c3 <- results_all$iteration_1$cohort_3$est_result$thetas
# results_i4c3 <- results_all$iteration_4$cohort_3$est_result$thetas
# CL_vec <-1:100
# V_vec <- 1:100
# tm50_vec <- 1:100
# HILL_vec <- 1:100


WT <- c(6.906,8.4,9.29,11.11,13.72,15.96,18.14,21.15,23.96,26.74,31.59,36.03,40.53,47.06,51.908,58.02,62.78,66.96,69.11,71.48,73.49)
AGE <- get_age_from_median_weight(WT)
AGE_WT <- data.frame(AGE,WT)




if(max(as.numeric(run_summary$cohort))>1) run_summary <- subset(run_summary, run_summary$cohort!=1) 

final_estimates <- data.frame(iter=1:100, CL=rep(0,100),V=rep(0,100),TM50=rep(0,100),HILL=rep(1,100),cohort=rep(0,100))


d9_CI_CL <- data.frame(iter=rep(1:100,rep(21,100)),pma=rep(AGE$PMA,100),d9_lower=rep(0,2100),d9_higher=rep(0,2100),
                       np=rep(0,2100))
d9_CI_V <-d9_CI_CL 


for( i in 1:100){
  load(paste(res_dir,"results_rep_",i,".Rdata",sep=""))
  tmp <- data.frame(aod_res$est_summary)
  
  
  final_estimates[i,2:4] <-tmp[length(tmp[,1]),2:4] 
  final_estimates[i,6] <-max(as.numeric(tmp$cohort))
  
  iteration <- aod_res
  
  x  <- iteration$final_design$x[-1]
  gs <- iteration$final_design$groupsize[-1]
  
  for(j in 1:length(x)){
    group_distribution[x[j],2] <- group_distribution[x[j],2] + gs[j]
    
  }
  
  if(length(iteration$stop_res)>0){
    stop_res<-iteration$stop_res[[length(iteration$stop_res)]]
    
    final_cohort_crit_res_CL <- stop_res[[3]]
    final_cohort_crit_res_V <- stop_res[[4]]
  }else{
    final_cohort_crit_res_CL <- matrix(nrow=21,ncol=4)
    final_cohort_crit_res_V <- matrix(nrow=21,ncol=4)
  }
  for(k in 1:21){
    #     
    age = AGE_WT$PMA[k]
    wt = AGE_WT$WT[k]
    
    d9_CI_CL[k+21*(i-1),5] <- size_mat_scaling_log(final_estimates[i,2:4],age,wt)[1]
    d9_CI_V[k+21*(i-1),5]  <- size_mat_scaling_log(final_estimates[i,2:4],age,wt)[2]
    #     
    #  
    
    d9_CI_CL[k+21*(i-1),3:4]= final_cohort_crit_res_CL[k,3:4] #c(exp(-qt((1-alpha/2),df)*SElcl),exp(qt((1-alpha/2),df)*SElcl))
    d9_CI_V[k +21*(i-1),3:4]= final_cohort_crit_res_V[k,3:4] #c(exp(-qt((1-alpha/2),df)*SE_V),exp(qt((1-alpha/2),df)*SE_V))
    
  }
  print(i)
}


df_CL_coverage <- data.frame(pma=AGE_WT$PMA,d9_ll=rep(0,21),d9_lh=rep(0,21),
                             d9_hl=rep(0,21),d9_hh=rep(0,21),
                             np_lower=rep(0,21),np_higher=rep(0,21))

df_V_coverage <- df_CL_coverage

for(i in 1:21){
  tmp2 <- subset(d9_CI_CL,pma==AGE_WT$PMA[i])
  df_CL_coverage[i,2:3] <- c(min(tmp2$d9_lower,na.rm=T),max(tmp2$d9_lower,na.rm=T))#quantile((tmp2$d9_lower),probs=c(0.025,0.975),na.rm=T)
  df_CL_coverage[i,4:5] <- c(min(tmp2$d9_higher,na.rm=T),max(tmp2$d9_higher,na.rm=T))#quantile((tmp2$d9_higher),probs=c(0.025,0.975),na.rm=T)
  df_CL_coverage[i,6:7] <- quantile((exp(tmp2$np)),probs=c(0.025,0.975),na.rm=T)/gm_mean(exp(tmp2$np))
  
  
  tmp3 <- subset(d9_CI_V,pma==AGE_WT$PMA[i])
  df_V_coverage[i,2:3] <- c(min(tmp3$d9_lower,na.rm=T),max(tmp3$d9_lower,na.rm=T))#quantile((tmp3$d9_lower),probs=c(0.025,0.975),na.rm=T)
  df_V_coverage[i,4:5] <- c(min(tmp3$d9_higher,na.rm=T),max(tmp3$d9_higher,na.rm=T))#quantile((tmp3$d9_higher),probs=c(0.025,0.975),na.rm=T)
  df_V_coverage[i,6:7] <- quantile((tmp3$np),probs=c(0.025,0.975),na.rm=T)/mean((tmp3$np))
}

df_CL_coverage_m<-data.frame(pma=rep(AGE_WT$PMA,3),Legend=c(rep("black",42),rep("grey",21)),
                             lower=c(df_CL_coverage[,2],df_CL_coverage[,4],df_CL_coverage[,6]),
                             higher=c(df_CL_coverage[,3],df_CL_coverage[,5],df_CL_coverage[,7]))

df_V_coverage_m<-data.frame(pma=rep(AGE_WT$PMA,3),Legend=c(rep("Predicted CI Limits",42),rep("CI of Estimates",21)),
                            lower=c(df_V_coverage[,2],df_V_coverage[,4],df_V_coverage[,6]),
                            higher=c(df_V_coverage[,3],df_V_coverage[,5],df_V_coverage[,7]))


x_step=c(rep(seq(0.8,20.8,1),2),seq(1.2,21.2,1))

#950x400
coverage_CI_CL <- ggplot(data=(df_CL_coverage_m),aes(x=x_step,y=lower)) +
  geom_errorbar(aes(ymin=lower*100, ymax=higher*100,colour=Legend)) +
  xlab("Pediatric Sub Group") + ylab("% of geometric mean")+ggtitle("LLM: CI coverage of CLp")+
  theme_bw()+ theme(axis.line = element_line(colour = "black"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank())+
  
  scale_x_discrete(breaks=1:21,labels=(AGE_WT$PMA))+scale_y_continuous(breaks=seq(50,150,10),limits=c(50,150))
coverage_CI_CL

coverage_CI_V <- ggplot(data=(df_V_coverage_m),aes(x=x_step,y=lower)) +
  geom_errorbar(aes(ymin=lower*100, ymax=higher*100,colour=Legend)) +
  xlab("Pediatric Sub Group") + ylab("% of geometric mean")+ggtitle("LLM: CI coverage of CLp")+
  theme_bw()+ theme(axis.line = element_line(colour = "black"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank())+
  
  scale_x_discrete(breaks=1:21,labels=(AGE_WT$PMA))+scale_y_continuous(breaks=seq(50,150,10),limits=c(50,150))
coverage_CI_V




coverage_CI_CL <- ggplot(data=(df_CL_coverage_m),aes(x=x_step,y=lower)) +
  geom_errorbar(aes(ymin=lower*100, ymax=higher*100,colour=Legend)) +
  xlab("PMA") + ylab("% of geometric mean")+ggtitle("LLM: CI coverage of CLp")+
  theme(axis.text = element_text(size=20,colour="black"),axis.title = element_text(size=23,face="bold"),
        axis.text.x = element_text(angle = 60,vjust = 0.65),plot.title=element_text(size=23,face="bold"))+
  scale_x_discrete(breaks=1:21,labels=(AGE_WT$PMA))+scale_y_continuous(breaks=seq(50,150,10),limits=c(50,150))
coverage_CI_CL

