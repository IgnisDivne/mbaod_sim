

setwd("C:/Users/Eric Stromberg/Desktop/MBAOD_project")

# remove things from the global environment
rm(list=ls())

library(PopED)
library(mvtnorm)
library(xpose4)

# load the MBAOD package
devtools::load_all("C:/Users/Eric Stromberg/Desktop/WarfarinPKPD/R/MBAOD")

# load the PopED model file
source("PopED_files/poped.mod.PK.1.comp.maturation_Xlog.R")
# load the weight estimator file
source("get_weight.R")
source("stop_critX.R")



#Model based prior
step_1=list(
  design = list(
    groupsize = 1000,
    x = c(age_group=c(7)),#6=adults ##t(rbind(PMA=c(1346)),WT=c(70))),  #AGE in PMA (weeks) assuming normal 39 week pregnancy
    xt = c(0.1, 2, 6, 12, 24)
  ),
  optimize=NULL,
  simulate=list(target="NONMEM", model="C:/Users/Eric Stromberg/Desktop/MBAOD_project/NONMEM_files/sim_log_NCA.mod",
                data=list(dosing = list(list(AMT=1000,Time=0)),
                          manipulation = list(expression(PMA <- age_group_2_PMA(ID,age_group)),
                                              expression(WT <- get_weight(ID,PMA,probs = c(0.1,0.9))),
                                              expression(AMT <- AMT*WT/70)
                          )
                )
                
  ),
  estimate=list(target="NONMEM", model="C:/Users/Eric Stromberg/Desktop/MBAOD_project/NONMEM_files/est_red_log_NCA.mod")
)

results_nca_sim <- mbaod_simulate(cohorts=list(step_1), # anything after step_3 is the same as step_3
                                      ncohorts=1, # number of steps or cohorts in one AOD
                                      rep=1, #number of times to repeat the MBAOD simulation 
                                      name="nca_sim_", 
                                      zip_directories=F,
                                      description="Simulated previous adult study",
                                      seednr=321)

NCA_PARAMS <- read.table("C:/Users/Eric Stromberg/Desktop/MBAOD_project/nca_sim__run_dir_NA/rep_1/cohort_1/mc_sim_1.tab", quote="\"")

sdLCL <- sd(NCA_PARAMS$V8)
power_adult_est <- get_power(sdLCL)
cvLCL <- sd(NCA_PARAMS$V7)/exp(1.01)

trans_sdLCL <- sqrt(log(cvLCL^2 +1))




######################################################
#Estimated number of individuals based on adult parameter estimates

#and a guess of scaling (same initials as for the MBAOD and OD)
######################################################


age_group <- 1:6
df <- data.frame(ID=1:60000)
df$AGE_GROUP <- sort(rep(age_group,10000))
df$PMA <- age_group_2_PMA(df$ID,df$AGE_GROUP)
df$WT <- get_weight(df$ID,df$PMA,probs = c(0.1,0.9))

for(j in 1:3){
  scaling_power_CL_tmp <- data.frame(sd = rep(0,6),nsub = rep(0,6), power=rep(0,6))
  scaling_power_V_tmp <-  data.frame(sd = rep(0,6),nsub = rep(0,6), power=rep(0,6))
  
  df_tmp <-df
  df_tmp$CL=unlist(param_sim(TM50=100-25*(j-1),PMA=df$PMA,WT=df$WT)[1])
  df_tmp$V=unlist(param_sim(TM50=100-25*(j-1),PMA=df$PMA,WT=df$WT)[2])
  
  for(i in 1:length(unique(df$AGE_GROUP))){
   
   sd_tmp_CL <- sd(subset(df_tmp,df_tmp$AGE_GROUP==i)$CL)
   power_tmp_CL <- get_power(sd_tmp_CL) 
   scaling_power_CL_tmp[i,2:3]  <- subset(power_tmp_CL,power_tmp_CL$power==min(subset(power_tmp_CL, power_tmp_CL$power >= 0.8)))
   scaling_power_CL_tmp[i,1] <- sd_tmp_CL 
   
   sd_tmp_V <- sd(subset(df_tmp,df_tmp$AGE_GROUP==i)$V)
   power_tmp_V <- get_power(sd_tmp_V) 
   scaling_power_V_tmp[i,2:3] <- subset(power_tmp_V,power_tmp_V$power==min(subset(power_tmp_V, power_tmp_V$power >= 0.8)))
   scaling_power_V_tmp[i,1] <- sd_tmp_V
  
  }
  if(j==1){
    scaling_power_CL <- scaling_power_CL_tmp
    scaling_power_V <- scaling_power_V_tmp
  }else{
    scaling_power_CL<- rbind(scaling_power_CL,scaling_power_CL_tmp)
    scaling_power_V<- rbind(scaling_power_V,scaling_power_V_tmp)
  }
}

scaling_power_CL$Misspecification <- c(rep("None",6),rep("Small",6),rep("Large",6))
scaling_power_CL$Age_Group <- rep(1:6,3)

scaling_power_V$Misspecification <-scaling_power_CL$Misspecification
scaling_power_V$Age_Group <-scaling_power_CL$Age_Group

param_sim <- function(CL=1.0,V=3,omCL=0.05,omV=0.05,TM50=75,PMA,WT,p=F){
  if(p==F){
  LCL <-rnorm(length(PMA),CL,sqrt(omCL))
  LV <- rnorm(length(PMA),V,sqrt(omV))
  }else{
    LCL <-qnorm(p,CL,sqrt(sdCL))
    LV <- qnorm(p,V,sqrt(sdV))
    
  }
  LCLi=LCL+0.75*log(WT/70)+log((PMA)/(PMA + TM50))
  LVi=LV+log(WT/70)
  return(list(LCLi,LVi))
  
}

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
    for(j in 1:length(groups)){
      nID_group[i,groups[j]] <- nID_group[i,groups[j]] + nID[j] 
    }
  }
  df_nID <- data.frame(Iteration = rep_nr, n_children,age_group)
  return(list(df_nID,nID_tot))
}


#############################################
#########Total Number of children plots######
############################################
#mbaod true
res_dir <- ("C:/Users/Eric Stromberg/Desktop/MBAOD_project/bridging_maturation_model_xoptim_restricted_group_true_run_dir_1/results_all.Rdata")
df_true <- unlist(get_nID(res_dir,n_iter=100)[2])

#mbaod small
res_dir <- c("C:/Users/Eric Stromberg/Desktop/MBAOD_project/bridging_maturation_model_xoptim_restricted_group_small_fixed_run_dir_3/results_all.Rdata")
df_small <- unlist(get_nID(res_dir,n_iter=100)[2])

#mbaod large
res_dir <-c("C:/Users/Eric Stromberg/Desktop/MBAOD_project/bridging_maturation_model_xoptim_restricted_group_large_fixed_run_dir_2/results_all.Rdata")  
df_large <- unlist(get_nID(res_dir,n_iter=100)[2])

###############################SIM CRIT#######################
res_dir <- ("C:/Users/Eric Stromberg/Desktop/MBAOD_project/bridging_maturation_model_xoptim_restricted_group_true_run_dir_1/results_all.Rdata")
df_true <- unlist(get_nID(res_dir,n_iter=100)[2])

#mbaod small
res_dir <- c("C:/Users/Eric Stromberg/Desktop/MBAOD_project/bridging_maturation_model_xoptim_restricted_group_small_fixed_run_dir_3/results_all.Rdata")
df_small <- unlist(get_nID(res_dir,n_iter=100)[2])

#mbaod large
res_dir <-c("C:/Users/Eric Stromberg/Desktop/MBAOD_project/bridging_maturation_model_xoptim_restricted_group_large_fixed_run_dir_2/results_all.Rdata")  
df_large <- unlist(get_nID(res_dir,n_iter=100)[2])



#####Restricted OD#########
#OD_small
file_dir <- "C:/Users/Eric Stromberg/Desktop/MBAOD_project/bridging_maturation_model_Xlog_small_OD_run_dir_3"
df_OD_small <- unlist(get_nID(file_dir,n_iter=20,full_path=F)[2])

#OD Large
file_dir <- "C:/Users/Eric Stromberg/Desktop/MBAOD_project/bridging_maturation_model_Xlog_OD_large_run_dir_3"
df_OD_large <- unlist(get_nID(file_dir,n_iter=10,full_path=F)[2])

#OD true
file_dir <- "C:/Users/Eric Stromberg/Desktop/MBAOD_project/bridging_maturation_model_Xlog_OD_true_run_dir_1"
df_OD_true <- unlist(get_nID(file_dir,n_iter=100,full_path=F)[2])
df_od_nid <- get_nID(file_dir,n_iter=100,full_path=F)[[1]]


##From adult parameters and guess of scaling
nca_true <- max(sum(subset(scaling_power_CL,scaling_power_CL$Misspecification=="None")$nsub),
                 sum(subset(scaling_power_CL,scaling_power_V$Misspecification=="None")$nsub))
nca_small <- max(sum(subset(scaling_power_CL,scaling_power_CL$Misspecification=="Small")$nsub),
                 sum(subset(scaling_power_CL,scaling_power_V$Misspecification=="Small")$nsub))
nca_large <- max(sum(subset(scaling_power_CL,scaling_power_CL$Misspecification=="Large")$nsub),
                 sum(subset(scaling_power_CL,scaling_power_CL$Misspecification=="Large")$nsub))

##From adult sd(logCLi)
nca_adult_est <- 36


df_box_true <- data.frame("True MBAOD"=df_true,"True OD"=df_OD_true,"True Scaling"= nca_true,
                           "Adult Prior"=nca_adult_est)
df_box_small <- data.frame("Small MBAOD"=df_small,"Small OD"=df_OD_small,"Small Scaling"= nca_small,
"Adult Prior"=nca_adult_est)
df_box_large <- data.frame("Large MBAOD"=df_large,"Large OD"=df_OD_large,"Large Scaling"= nca_large,
                           "Adult Prior"=nca_adult_est)



library(reshape2)
df_box_true <- melt(df_box_true)
df_box_true <- subset(df_box_true,df_box_true$value != -100)

df_box_small <- melt(df_box_small)
df_box_small <- subset(df_box_small,df_box_small$value != -100)

df_box_large <- melt(df_box_large)
df_box_large <- subset(df_box_large,df_box_large$value != -100)
library(ggplot2)
p_box_true <- ggplot(data=df_box_true,aes(x=variable,y=value)) + geom_boxplot() +
  xlab(" ") + ylab("N Children in total")+ ggtitle("No Misspecification")+
  scale_x_discrete(labels = c("MBAOD","OD","Scaling","Adult Prior"))+scale_y_continuous(breaks=seq(0, 45, 5),limits=c(0,45))+
  theme(axis.text = element_text(size=20,colour="black"),axis.title = element_text(size=25,face="bold"),
        axis.text.x = element_text(angle = 45, vjust = 0.65),plot.title=element_text(size=25,face="bold"))
p_box_true


p_box_small <- ggplot(data=df_box_small,aes(x=variable,y=value)) + geom_boxplot() +
  xlab("Design Approach") + ylab("N Children in total")+ ggtitle("Small Misspecification")+
  scale_x_discrete(labels = c("MBAOD","OD","Scaling","Adult Prior"))+scale_y_continuous(breaks=seq(0, 45, 5),limits=c(0,45))+
  theme(axis.text = element_text(size=20,colour="black"),axis.title = element_text(size=25,face="bold"),
        axis.text.x = element_text(angle = 45, vjust = 0.65),plot.title=element_text(size=25,face="bold"))
p_box_small

p_box_large <- ggplot(data=df_box_large,aes(x=variable,y=value)) + geom_boxplot() +
  xlab(" ") + ylab("N Children in total")+ ggtitle("Large Misspecification")+
  scale_x_discrete(labels = c("MBAOD","OD","Scaling","Adult Prior"))+scale_y_continuous(breaks=seq(0, 45, 5),limits=c(0,45))+
theme(axis.text = element_text(size=20,colour="black"),axis.title = element_text(size=25,face="bold"),
      axis.text.x = element_text(angle = 45, vjust = 0.65),plot.title=element_text(size=25,face="bold"))
p_box_large



#############################################
#########Number of Children per Cohort######
############################################
#mbaod true
res_dir <- ("C:/Users/Eric Stromberg/Desktop/MBAOD_project/bridging_maturation_model_xoptim_restricted_group_true_run_dir_1/results_all.Rdata")



#mbaod small
res_dir<- ("C:/Users/Eric Stromberg/Desktop/MBAOD_project/bridging_maturation_model_xoptim_restricted_group_small_run_dir_2/results_all.Rdata")


#mbaod large
res_dir<- ("C:/Users/Eric Stromberg/Desktop/MBAOD_project/bridging_maturation_model_xoptim_restricted_group_large_run_dir_2/results_all.Rdata")  

#####Restricted OD#########
#OD_small

file_dir <- "C:/Users/Eric Stromberg/Desktop/MBAOD_project/bridging_maturation_model_Xlog_OD_run_dir_NA"



#OD Large

file_dir <- "C:/Users/Eric Stromberg/Desktop/MBAOD_project/bridging_maturation_model_Xlog_OD_large_run_dir_1"


for(i in 1:6){
  df_box <- data.frame("Small MBAOD"=df_mbaod_small[i],"Large MBAOD" = df_mbaod_large[i])
  df_box <- melt(df_box)
  df_box <- subset(df_box,df_box$value != -100)
  p_box <- ggplot(data=df_box,aes(x=variable,y=value)) + geom_boxplot() +
    xlab("Misspecification on TM50") + ylab("N Children in total")
  
  p[i] <- p_box
}


#df_mbaod_small <- data.frame(Iteration, n_children,age_group)
#min iterations
df1 <- subset(df_mbaod_small,df_mbaod_small$Iteration==c(25))

#Max iterations
df2 <- subset(df_mbaod_small,df_mbaod_small$Iteration==c(49))
#Median
df3 <- subset(df_mbaod_small,df_mbaod_small$Iteration==50)#as.integer(rownames(subset(df_small,df_small$nID == 31))))
#Combined
df <- rbind(df1,df2,df3)
df
##Full Data
df <- df_mbaod_small
age <- c("3 to <6 months", "6 to <12 months", "1 to <2 years","2 to <6 years","6 to <12 years","12 to <18 years" )
p4 <- ggplot(data=df,aes(x=n_children,y=age_group,group=factor(Iteration),color=factor(Iteration))) +
  geom_point() + geom_line() + scale_color_manual(values=c("black","red","blue"),guide=FALSE)+
  xlab("Total Number of Included Children") + ylab("Age Group")+
  theme(axis.text = element_text(size=12,colour="black"),axis.title = element_text(size=15,face="bold"))+
  scale_y_continuous(breaks=1:6,labels=age) + scale_x_continuous(breaks=seq(10,70,5))
p4



######Plot of how CL varies with median weight and median age##########
median_weights <- c(6.906,8.4,9.29,11.11,13.72,15.96,18.14,21.15,23.96,26.74,31.59,36.03,40.53,47.06,51.908,58.02,62.78,66.96,69.11,71.48,73.49,78.367)

AGE <- get_age_from_median_weight(median_weights)$PMA

df_full <- data.frame(PMA = AGE,Weight = median_weights,CL2.5= 1:length(AGE),CL_median = 1:length(AGE),CL97.5= 1:length(AGE),
                      V2.5= 1:length(AGE),V_median = 1:length(AGE),V97.5= 1:length(AGE))

for(i in 1:length(median_weights)){
 df_full[i,c(4,7)]<- exp(unlist(param_sim(PMA=AGE[i],WT=median_weights[i],p=0.5)))
 df_full[i,c(3,6)]<- exp(unlist(param_sim(PMA=AGE[i],WT=median_weights[i],p=0.025)))
 df_full[i,c(5,8)]<- exp(unlist(param_sim(PMA=AGE[i],WT=median_weights[i],p=0.975)))

}


#####################
library(reshape2)

#dataframe with only CL values
df <- df_full[,1:5]

p5 <- ggplot(data=df,aes(x=PMA,y=CL_median)) + ggtitle("Scaling of Clearance")+
   xlab("PMA") + ylab("Clearance")+geom_ribbon(x=df$PMA,ymin=df$CL2.5,ymax=df$CL97.5,fill = "red", alpha = 0.4) +geom_line()+
  scale_y_continuous(breaks=seq(0, 5, 0.5),limits=c(0,5)) + scale_x_continuous(breaks=seq(0, 1400, 250),limits=c(0,1300))+
  theme(axis.text = element_text(size=30,colour="black"),axis.title = element_text(size=30,face="bold"),plot.title=element_text(size=30,face="bold"))  
p5

#dataframe with only V values
df <- df_full[,c(1:2,6:8)]

p6 <- ggplot(data=df,aes(x=PMA,y=V_median)) + ggtitle("Scaling of Volume of Distribution")+
  xlab("PMA") + ylab("Volume of Distribution")+geom_ribbon(x=df$PMA,ymin=df$V2.5,ymax=df$V97.5,fill = "red", alpha = 0.4) +geom_line()+
  scale_y_continuous(breaks=seq(0, 35, 5),limits=c(0,35)) + scale_x_discrete(breaks=seq(0, 1400, 250),limits=c(0,1300))+
  theme(axis.text = element_text(size=30,colour="black"),axis.title = element_text(size=30,face="bold"),plot.title=element_text(size=30,face="bold"))


p6

