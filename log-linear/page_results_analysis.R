
##Set path to directory of this file
setwd("C:/Users/erist836/Desktop/To Andy/MBAOD_project")

# remove things from the global environment
rm(list=ls())

library(PopED)
library(mvtnorm)
library(xpose4)

# load the MBAOD package (CHange to your path)
devtools::load_all("C:/Users/erist836/Desktop/WarfarinPKPD/R/MBAOD")

# load the PopED model file
source("PopED_files/poped.mod.PK.1.comp.maturation_Xlog.R")
# load the weight estimator file
source("get_weight.R")
source("stop_critX.R")
source("plot_functions.R")


##########################################################################
###Simulate the Adult CL and V to get estimates of SD in all age groups###
##########################################################################
step_1=list(
  design = list(
    groupsize = 1000,
    x = c(age_group=c(7)),#6=adults ##t(rbind(PMA=c(1346)),WT=c(70))),  #AGE in PMA (weeks) assuming normal 39 week pregnancy
    xt = c(0.1, 2, 6, 12, 24)
  ),
  optimize=NULL,
  simulate=list(target="NONMEM", model="NONMEM_files/sim_log_NCA.mod",
                data=list(dosing = list(list(AMT=1000,Time=0)),
                          manipulation = list(expression(PMA <- age_group_2_PMA(ID,age_group)),
                                              expression(WT <- get_weight(ID,PMA,probs = c(0.1,0.9))),
                                              expression(AMT <- AMT*WT/70)
                          )
                )
                
  ),
  estimate=list(target="NONMEM", model="NONMEM_files/est_red_log_NCA.mod")
)

results_nca_sim <- mbaod_simulate(cohorts=list(step_1), # anything after step_3 is the same as step_3
                                      ncohorts=1, # number of steps or cohorts in one AOD
                                      rep=1, #number of times to repeat the MBAOD simulation 
                                      name="nca_sim_test", 
                                      zip_directories=F,
                                      description="Simulated previous adult study",
                                      seednr=321)




###############################################################################################################
####TO GENERATE PLOTS FOR THE RESULTS PRESENTED AT PAGE (in folder Page_results), RUN THE CODE FROM HERE#######
###############################################################################################################
##Read the table with simulated adult parameters
NCA_PARAMS <- read.table("nca_sim_test_run_dir_NA/rep_1/cohort_1/mc_sim_1.tab", quote="\"")
##Calculate the Standard deviation of CL
sdLCL <- sd(NCA_PARAMS$V8)
##Get the power for different nsub
power_adult_sd_est <- get_power(sdLCL)
##get the min nsub => power >=0.8
nsub_adult_sd_est <- subset(power_adult_sd_est,power_adult_sd_est$power==min(subset(power_adult_sd_est$power,power_adult_sd_est$power>=0.8)))$nsub
##Calulate the log transformed based on CV instead of SD according to the paper of Wang et al.
cvLCL <- sd(NCA_PARAMS$V7)/exp(1.01)
L_trans_sdLCL <- sqrt(log(cvLCL^2 +1))
power_adult_L_trans_est <- get_power(L_trans_sdLCL)
nsub_adult_L_trans_est <- subset(power_adult_L_trans_est,power_adult_L_trans_est$power==min(subset(power_adult_L_trans_est$power,power_adult_L_trans_est$power>=0.8)))$nsub

##Nsub for using SD of CL and SD transformed according to wang et al is the same...6 per group


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
  if(j==1) tm50 <- 100
  if(j==2) tm50 <- 75
  if(j==3) tm50 <- 150
  
  
  df_tmp <-df
  df_tmp$CL=unlist(param_sim(TM50=tm50,PMA=df$PMA,WT=df$WT)[1])
  df_tmp$V=unlist(param_sim(TM50=tm50,PMA=df$PMA,WT=df$WT)[2])
  
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



##Read the total number of children from MBAOD and OD
#############################################
#########Total Number of children plots######
############################################
#mbaod true
res_dir <- ("Page_results/bridging_maturation_model_xoptim_restricted_group_true_run_dir_1/results_all.Rdata")
df_true <- unlist(get_nID(res_dir,n_iter=100)[2])

#mbaod small
res_dir <- c("Page_results/bridging_maturation_model_xoptim_restricted_group_small_new_run_dir_1/results_all.Rdata")
df_small <- unlist(get_nID(res_dir,n_iter=100)[2])

#mbaod large
res_dir <-c("Page_results/bridging_maturation_model_xoptim_restricted_group_large2_run_dir_1/results_all.Rdata")  
df_large <- unlist(get_nID(res_dir,n_iter=100)[2])



#####Restricted OD#########
#OD_small
file_dir <- "Page_results/bridging_maturation_model_Xlog_small_OD_run_dir_3"
df_OD_small <- unlist(get_nID(file_dir,n_iter=20,full_path=F)[2])

#OD Large
file_dir <- "Page_results/bridging_maturation_model_Xlog_OD_large2_run_dir_1"
df_OD_large <- unlist(get_nID(file_dir,n_iter=25,full_path=F)[2])

#OD true
file_dir <- "Page_results/bridging_maturation_model_Xlog_OD_true_run_dir_1"
df_OD_true <- data.frame(unlist(get_nID(file_dir,n_iter=100,full_path=F)[2]))
df_od_nid <- get_nID(file_dir,n_iter=100,full_path=F)[[1]]



 

##From adult parameters and guess of scaling
nca_true <- max(sum(subset(scaling_power_CL,scaling_power_CL$Misspecification=="None")$nsub),
                 sum(subset(scaling_power_CL,scaling_power_V$Misspecification=="None")$nsub))
nca_small <- max(sum(subset(scaling_power_CL,scaling_power_CL$Misspecification=="Small")$nsub),
                 sum(subset(scaling_power_CL,scaling_power_V$Misspecification=="Small")$nsub))
nca_large <- max(sum(subset(scaling_power_CL,scaling_power_CL$Misspecification=="Large")$nsub),
                 sum(subset(scaling_power_CL,scaling_power_CL$Misspecification=="Large")$nsub))

##From adult sd(logCLi)
nca_adult_est <- nsub_adult_sd_est

##Create the data frames
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

##Plot the total number of IDS
p_box_true <- ggplot(data=df_box_true,aes(x=variable,y=value)) + geom_boxplot() +
  xlab(" ") + ylab("N Children in total")+ ggtitle("No Misspecification")+
  scale_x_discrete(labels = c("MBAOD","OD","Scaling","Adult Prior"))+scale_y_continuous(breaks=seq(0, 45, 5),limits=c(0,45))+
  theme(axis.text = element_text(size=20,colour="black"),axis.title = element_text(size=23,face="bold"),
        axis.text.x = element_text(angle = 45, vjust = 0.65),plot.title=element_text(size=23,face="bold"))
p_box_true


p_box_small <- ggplot(data=df_box_small,aes(x=variable,y=value)) + geom_boxplot() +
  xlab("Design Approach") + ylab("N Children in total")+ ggtitle("Small Misspecification")+
  scale_x_discrete(labels = c("MBAOD","OD","Scaling","Adult Prior"))+scale_y_continuous(breaks=seq(0, 45, 5),limits=c(0,45))+
  theme(axis.text = element_text(size=20,colour="black"),axis.title = element_text(size=23,face="bold"),
        axis.text.x = element_text(angle = 45, vjust = 0.65),plot.title=element_text(size=23,face="bold"))
p_box_small

p_box_large <- ggplot(data=df_box_large,aes(x=variable,y=value)) + geom_boxplot() +
  xlab(" ") + ylab("N Children in total")+ ggtitle("Large Misspecification")+
  scale_x_discrete(labels = c("MBAOD","OD","Scaling","Adult Prior"))+scale_y_continuous(breaks=seq(0, 45, 5),limits=c(0,45))+
theme(axis.text = element_text(size=20,colour="black"),axis.title = element_text(size=23,face="bold"),
      axis.text.x = element_text(angle = 45, vjust = 0.65),plot.title=element_text(size=23,face="bold"))
p_box_large
##Plots done


#############################################
#########MBAOD and OD Group Selection Plots#####
############################################


##### Plot for group selection of Unrestricted OD#########

#####Read Data and create data.frames
#OD_small
file_dir <- "Page_results/bridging_maturation_model_Xlog_small_OD_unrestricted_small_run_dir_1"
df_OD_small <- data.frame(get_nID(file_dir,n_iter=9,full_path=F)[1], Miss = "Small")
#OD Large
file_dir <- "Page_results/bridging_maturation_model_Xlog_small_OD_unrestricted_large2_run_dir_1/results_all.Rdata"
df_OD_large <- data.frame(get_nID(file_dir,n_iter=100,full_path=T)[1],Miss = "Large")
#OD true
file_dir <- "Page_results/bridging_maturation_model_Xlog_small_OD_unrestricted_true_run_dir_1/results_all.Rdata"
df_OD_true <- data.frame(get_nID(file_dir,n_iter=100,full_path=T)[[1]],Miss = "True")


df_adult_true <- data.frame(Iteration = 1:max(df_OD_true$Iteration), n_children = rep(0,max(df_OD_true$Iteration)), 
                            age_group=rep(7,max(df_OD_true$Iteration)), Miss = rep("True", max(df_OD_true$Iteration)))
df_adult_large <- data.frame(Iteration = 1:max(df_OD_large$Iteration), n_children = rep(0,max(df_OD_large$Iteration)), 
                            age_group=rep(7.1,max(df_OD_large$Iteration)), Miss = rep("Large", max(df_OD_large$Iteration)))
df_adult_small <- data.frame(Iteration = 1:max(df_OD_small$Iteration), n_children = rep(0,max(df_OD_small$Iteration)), 
                            age_group=rep(6.9,max(df_OD_small$Iteration)), Miss = rep("Small", max(df_OD_small$Iteration)))
df_adult <- rbind(df_adult_true,df_adult_small,df_adult_large)
df_OD_small$age_group <- df_OD_small$age_group -0.1 
df_OD_large$age_group <- df_OD_large$age_group +0.1 

df <- rbind(df_OD_true, df_OD_small, df_OD_large)
df <- rbind(df,df_adult)
df$ID <- as.integer(df$Iteration) + as.integer(df$Miss)*10000


age <- c("3 to <6 months", "6 to <12 months", "1 to <2 years","2 to <6 years","6 to <12 years","12 to <18 years","Adults" )
#####End data.frame creations#######

####Plot the data####
unres_OD_plot <- ggplot(data=df,aes(x=n_children,y=age_group,group=ID,color=as.factor(Miss),guide=T)) +
  geom_point(aes(size=2,guide=F)) + geom_line(size=0.75) + scale_color_manual(values=rep(c("Black","Blue","Red"),50),name  ="TM50",
                                                                      breaks=c("True", "Small","Large"),
                                                                      labels=c("100 (T)", "75","150"))+
  xlab("N Children in Total") + ylab("Age Group")+ ggtitle("Group Selection of Unrestricted OD")+
  theme(axis.text = element_text(size=15,colour="black"),axis.title = element_text(size=20,face="bold"),plot.title=element_text(size=20,face="bold"),
        legend.title = element_text(colour="black", size=16, face="bold"),
        legend.text = element_text(colour="black", size=15))+
  scale_y_continuous(limits=c(0.5,7.5),breaks=1:7,labels=age) + scale_x_continuous(breaks=seq(9,15,2))
unres_OD_plot
###End Plot
###End Unres. OD plot#####

#mbaod true
res_dir <- ("Page_results/bridging_maturation_model_xoptim_restricted_group_true_run_dir_1/results_all.Rdata")
df_true_MBAOD <- data.frame(get_nID(res_dir,n_iter=100)[1])
df_true_mbaod_tot<- data.frame(get_nID(res_dir,n_iter=100)[2])
df_true_mbaod_tot$iter <- 1:100

mbaod_max <- df_true_mbaod_tot[df_true_mbaod_tot[,1]==max(df_true_mbaod_tot[,1]),] #3, 43
mbaod_min <- df_true_mbaod_tot[df_true_mbaod_tot[,1]==min(df_true_mbaod_tot[,1]),] #1,12,64,66
mbaod_median <- df_true_mbaod_tot[df_true_mbaod_tot[,1]==median(df_true_mbaod_tot[,1]),] #2,7,14,17,18,26,29,31,35,37,38,40....

df_mbaod_min <- df_true_MBAOD[df_true_MBAOD$Iteration == c(66),]
df_mbaod_max <- df_true_MBAOD[df_true_MBAOD$Iteration == c(3),]
df_mbaod_med <- df_true_MBAOD[df_true_MBAOD$Iteration == c(17),]
df_mbaod_min$age_group <-df_mbaod_min$age_group -0.1 
df_mbaod_max$age_group <-df_mbaod_max$age_group +0.1
df_mbaod_min$Miss <- "Small"
df_mbaod_med$Miss <- "Med"
df_mbaod_max$Miss <- "Large"

df_adult_min <- data.frame(Iteration = 66, n_children = 0, age_group=6.9, Miss = "Small")
df_adult_med <- data.frame(Iteration = 17, n_children = 0, age_group=7, Miss = "Med")
df_adult_max <- data.frame(Iteration = 3, n_children = 0, age_group=7.1, Miss = "Large")



df_mbaod <- rbind(df_adult_min,df_mbaod_min,df_adult_med,df_mbaod_med,df_adult_max,df_mbaod_max)
df_mbaod$ID <- as.integer(df_mbaod$Iteration) + as.integer(df_mbaod$Miss)*10000

age <- c("3 to <6 months", "6 to <12 months", "1 to <2 years","2 to <6 years","6 to <12 years","12 to <18 years","Adults" )
res_mbaod_plot <- ggplot(data=df_mbaod,aes(x=n_children,y=age_group,group=ID, color=as.factor(Miss))) +
  geom_point(aes(size=2)) + geom_line(size=0.75) + scale_color_manual(values=rep(c("Black","Blue","Red"),50),name  ="Simulation Nr",
                                                                      breaks=c("Small", "Med","Large"),
                                                                      labels=c("66", "17","3"))+
  xlab("N Children in Total") + ylab("Age Group")+ ggtitle("Group Selection of Restricted MBAOD")+
  theme(axis.text = element_text(size=15,colour="black"),axis.title = element_text(size=20,face="bold"),plot.title=element_text(size=20,face="bold"),
        legend.title = element_text(colour="black", size=16, face="bold"),
        legend.text = element_text(colour="black", size=15))+
  scale_y_continuous(limits=c(0.5,7.5),breaks=1:7,labels=age) + scale_x_continuous(breaks=seq(9,30,2))
res_mbaod_plot


#OD true
file_dir <- "Page_results/bridging_maturation_model_Xlog_OD_true_run_dir_1"
df_OD_true <- data.frame(unlist(get_nID(file_dir,n_iter=100,full_path=F)[2]))
df_od_nid <- get_nID(file_dir,n_iter=100,full_path=F)[[1]]
df_od_nid <- df_od_nid[df_od_nid$Iteration == c(1),]
df_od_nid$age_group <-df_od_nid$age_group -0.2
df_od_nid$Miss <- "OD"
df_od_nid_adult <- data.frame(Iteration = 1, n_children = 0, age_group=6.8, Miss = "OD")
df_mbaod <- rbind(df_adult_min,df_mbaod_min,df_adult_med,df_mbaod_med,df_adult_max,df_mbaod_max,df_od_nid_adult,df_od_nid)
df_mbaod$ID <- as.integer(df_mbaod$Iteration) + as.integer(df_mbaod$Miss)*10000

res_mbaod_OD_plot <- ggplot(data=df_mbaod,aes(x=n_children,y=age_group,group=ID, color=as.factor(Miss))) +
  geom_point(aes(size=2)) + geom_line(size=0.75) + scale_color_manual(values=rep(c("Black","Blue","Red", "dark green"),50),name  ="Simulation Nr",
                                                                      breaks=c("Small", "Med","Large","OD"),
                                                                      labels=c("66", "17","3","OD"))+
  xlab("N Children in Total") + ylab("Age Group")+ ggtitle("Group Selection of Restricted MBAOD and OD")+
  theme(axis.text = element_text(size=15,colour="black"),axis.title = element_text(size=20,face="bold"),plot.title=element_text(size=20,face="bold"),
        legend.title = element_text(colour="black", size=16, face="bold"),
        legend.text = element_text(colour="black", size=15))+
  scale_y_continuous(limits=c(0.5,7.5),breaks=1:7,labels=age) + scale_x_continuous(breaks=seq(9,30,2))
res_mbaod_OD_plot
##########################################
####End Plots for group Selection#########
##########################################

######Plot of how CL varies with median weight and median age##########
median_weights <- c(5.1,6.906,8.4,9.29,11.11,13.72,15.96,18.14,21.15,23.96,26.74,31.59,36.03,40.53,47.06,51.908,58.02,62.78,66.96,69.11,71.48,73.49,78.367)

AGE <- get_age_from_median_weight(median_weights)$PMA
AGE <- c(0,AGE)
median_weights<-c(0,median_weights)

df_full <- data.frame(PMA = AGE,Weight = median_weights,CL_true= 1:length(AGE),CL_small = 1:length(AGE),CL_large= 1:length(AGE),
                      V_true= 1:length(AGE),V_small = 1:length(AGE),V_large= 1:length(AGE))

for(i in 1:length(median_weights)){
 df_full[i,c(4,7)]<- exp(unlist(param_sim(PMA=AGE[i],WT=median_weights[i],TM50=75,sdCL=0,sdV=0)))
 df_full[i,c(3,6)]<- exp(unlist(param_sim(PMA=AGE[i],WT=median_weights[i],TM50=100,sdCL=0,sdV=0)))
 df_full[i,c(5,8)]<- exp(unlist(param_sim(PMA=AGE[i],WT=median_weights[i],TM50=150,sdCL=0,sdV=0)))

}
df_full_mat <-  data.frame(PMA = AGE,Weight = median_weights,CL_true= 1:length(AGE),CL_small = 1:length(AGE),CL_large= 1:length(AGE),
                           V_true= 1:length(AGE),V_small = 1:length(AGE),V_large= 1:length(AGE))
for(i in 1:length(median_weights)){
  df_full_mat[i,c(4,7)]<- (unlist(param_sim_mat(PMA=AGE[i],WT=median_weights[i],TM50=75,sdCL=0,omV=0)))
  df_full_mat[i,c(3,6)]<- (unlist(param_sim_mat(PMA=AGE[i],WT=median_weights[i],TM50=100,sdCL=0,omV=0)))
  df_full_mat[i,c(5,8)]<- (unlist(param_sim_mat(PMA=AGE[i],WT=median_weights[i],TM50=150,sdCL=0,omV=0)))
  
}

#####################
library(reshape2)

#dataframe with only CL values
p <- NULL
colors <- c("Blue","black","red")
for (i in 1:3){
df <- df_full[,1:(2+i)]
df <- melt(df,id.vars=c("PMA","Weight"))

p <- list(p, ggplot(data=df,aes(x=PMA,y=value,group=as.factor(variable))) +geom_line(size=1.0,colour=colors[as.integer(df$variable)],guide=T) + ggtitle("Scaling of Clearance")+
             xlab("PMA") + ylab("Median Clearance")+ 
             scale_y_continuous(breaks=seq(0, 3, 0.5),limits=c(0,3)) + scale_x_continuous(breaks=seq(0, 1400, 250),limits=c(0,1300))+
             theme(axis.text = element_text(size=20,colour="black"),axis.title = element_text(size=23,face="bold"),plot.title=element_text(size=25,face="bold"),
                   legend.title = element_text(colour="black", size=16, face="bold"),
                   legend.text = element_text(colour="black", size=15))  
)
}

#dataframe with only V values
df <- df_full[,c(1:2,6:8)]
df <- melt(df,id.vars=c("PMA","Weight"))
p6 <- ggplot(data=df,aes(x=Weight,y=V_median)) + ggtitle("Scaling of Volume of Distribution")+
  xlab("PMA") + ylab("Volume of Distribution")+geom_ribbon(x=df$Weight,ymin=df$V2.5,ymax=df$V97.5,fill = "red", alpha = 0.4) +geom_line()+
  scale_y_continuous(breaks=seq(0, 35, 5),limits=c(0,35)) + scale_x_discrete(breaks=seq(0, 1400, 250),limits=c(0,1300))+
  theme(axis.text = element_text(size=30,colour="black"),axis.title = element_text(size=30,face="bold"),plot.title=element_text(size=30,face="bold"))


p6

library(reshape2)

#dataframe with only CL values
df <- df_full_mat[,1:5]
colors <-c("blue","black","red")
df <- melt(df,id.vars=c("PMA","Weight"))
p7 <- ggplot(data=df,aes(x=PMA,y=value*100,group=variable)) + ggtitle("Maturation of Clearance")+
  xlab("PMA") + ylab("% Maturation") +geom_line(size=1.5,colour=colors[as.integer(df$variable)])+
  scale_y_continuous(breaks=seq(0, 100, 25),limits=c(0,101)) + scale_x_continuous(breaks=seq(0, 1400, 250),limits=c(0,1300))+
  theme(axis.text = element_text(size=15,colour="black"),axis.title = element_text(size=20,face="bold"),plot.title=element_text(size=20,face="bold"))  
p7

##########################
###Power Evals results####
##########################

#No Misspecification results
load("Page_results/bridging_maturation_model_eval_OD_true_run_dir_3/results_all.Rdata")
stop_crit_res <- 1:100
for(i in 1:100){
  stop_crit_res[i] <- results_all[[paste("iteration_",i,sep="")]]$stop_res$cohort_1[[1]]
  
}
power_true_od <- sum(stop_crit_res)   #72



#Large Misspecification results
load("Page_results/bridging_maturation_model_eval_OD_large2_run_dir_1/results_all.Rdata")
stop_crit_res <- 1:100
for(i in 1:100){
  stop_crit_res[i] <- results_all[[paste("iteration_",i,sep="")]]$stop_res$cohort_1[[1]]
  
}
power_large_od <- sum(stop_crit_res) #70

