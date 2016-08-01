
#set path to the directory where this file is located
setwd("C:/Users/erist836/Desktop/MBAOD_project")

# remove things from the global environment
rm(list=ls())

library(PopED)
library(mvtnorm)
library(xpose4)

# load the MBAOD package change to your relevant load call
devtools::load_all("C:/Users/erist836/Documents/GitHub/MBAOD/R")

# load the PopED model file
source("PopED_files/poped.mod.PK.1.comp.maturation_real.R")
# load the weight estimator file
source("get_weight.R")
source("stop_critX_2.R")
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
  simulate=list(target="NONMEM", model="NONMEM_files/sim_propofol.mod",
                data=list(dosing = list(list(AMT=1000,Time=0)),
                          manipulation = list(expression(PMA <- age_group_2_PMA(ID,age_group)),
                                              expression(WT <- get_weight(ID,PMA,probs = c(0.1,0.9))),
                                              expression(AMT <- AMT*WT/70)
                          )
                )
                
  ),
  estimate=list(target="NONMEM", model="NONMEM_files/est_red_propofol_NCA.mod")
)

results_nca_sim <- mbaod_simulate(cohorts=list(step_1), # anything after step_3 is the same as step_3
                                      ncohorts=1, # number of steps or cohorts in one AOD
                                      rep=1, #number of times to repeat the MBAOD simulation 
                                      name="nca_sim_test_propofol", 
                                      zip_directories=F,
                                      description="Simulated previous adult study",
                                      seednr=321)




###############################################################################################################
####TO GENERATE PLOTS FOR THE RESULTS PRESENTED AT PAGE (in folder Page_results), RUN THE CODE FROM HERE#######
###############################################################################################################
##Read the table with simulated adult parameters
NCA_PARAMS <- read.table("nca_sim_test_propofol_run_dir_1/rep_1/cohort_1/NCA_PARAMS.tab", quote="\"")
##Calculate the Standard deviation of CL
sdLCL <- sd(NCA_PARAMS$V6)
##Get the power for different nsub
power_adult_sd_est <- get_power(sdLCL)
##get the min nsub => power >=0.8
nsub_adult_sd_est <- subset(power_adult_sd_est,power_adult_sd_est$power==min(subset(power_adult_sd_est$power,power_adult_sd_est$power>=0.8)))$nsub

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

for(j in 1:2){
  scaling_power_CL_tmp <- data.frame(sd = rep(0,6),nsub = rep(0,6), power=rep(0,6))
  scaling_power_V_tmp <-  data.frame(sd = rep(0,6),nsub = rep(0,6), power=rep(0,6))
  if(j==1) {
    tm50 <- 3.65065
    hill <- 1.526056
    }
  if(j==2) {
    tm50 <- 3.86073
    hill <- 1.223775
  } 
  
  
  df_tmp <-df
  df_tmp$CL=unlist(param_sim(TM50=tm50,HILL=hill,PMA=df$PMA,WT=df$WT)[1])
  df_tmp$V=unlist(param_sim(TM50=tm50,HILL=hill,PMA=df$PMA,WT=df$WT)[2])
  
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

scaling_power_CL$Misspecification <- c(rep("None",6),rep("Large",6))
scaling_power_CL$Age_Group <- rep(1:6,2)

scaling_power_V$Misspecification <-scaling_power_CL$Misspecification
scaling_power_V$Age_Group <-scaling_power_CL$Age_Group



##Read the total number of children from MBAOD and OD
#############################################
#########Total Number of children plots######
############################################
#mbaod true
res_dir <- ("C:/Users/erist836/Desktop/MBAOD_project/merge")
df_true <- unlist(get_nID(res_dir,n_iter=100,full_path = F)[2])


#mbaod large
res_dir <- ("C:/Users/erist836/Desktop/MBAOD_project/merge_gfr")
df_large <- unlist(get_nID(res_dir,n_iter=100,full_path = F)[2])



#####Restricted OD#########

#OD true
file_dir <- c("C:/Users/erist836/Desktop/MBAOD_project/propofol_odsc")
df_OD_true <- data.frame(unlist(get_nID(file_dir,n_iter=100,full_path=F)[2]))
df_od_nid <- get_nID(file_dir,n_iter=100,full_path=F)[[1]]


#OD Large
file_dir <- c("C:/Users/erist836/Desktop/MBAOD_project/gfr_odsc")
df_OD_large <- unlist(get_nID(file_dir,n_iter=25,full_path=F)[2])



 

##From adult parameters and guess of scaling
nca_true <- max(sum(subset(scaling_power_CL,scaling_power_CL$Misspecification=="None")$nsub),
                 sum(subset(scaling_power_CL,scaling_power_V$Misspecification=="None")$nsub))
nca_large <- max(sum(subset(scaling_power_CL,scaling_power_CL$Misspecification=="Large")$nsub),
                 sum(subset(scaling_power_CL,scaling_power_CL$Misspecification=="Large")$nsub))

##From adult sd(logCLi)
nca_adult_est <- 6*nsub_adult_sd_est


##Regular Optimal Design results
OD_LLT <-data.frame(Age_groups=c(6,3,2), nID=c(12,12,12))
OD_LLM <- OD_LLT   ###The designs did not depend on misspecification here


##Create the data frames
df_box_true <- data.frame("True MBAOD"=df_true,"True OD"=df_OD_true,"True Scaling"= nca_true,
                           "Adult Prior"=nca_adult_est, "Regular OD True" = sum(subset(OD_LLT,OD_LLT$Age_groups <7)$nID))


df_box_large <- data.frame("Large MBAOD"=df_large,"Large OD"=df_OD_large,"Large Scaling"= nca_large,
                           "Adult Prior"=nca_adult_est,"Regular OD Missp." = sum(subset(OD_LLM,OD_LLM$Age_groups <7)$nID))



library(reshape2)
df_box_true <- melt(df_box_true)
df_box_true <- subset(df_box_true,df_box_true$value != -100)


df_box_large <- melt(df_box_large)
df_box_large <- subset(df_box_large,df_box_large$value != -100)
library(ggplot2)

##Plot the total number of IDS
p_box_true <- ggplot(data=df_box_true,aes(x=variable,y=value)) + geom_boxplot() +
  xlab(" ") + ylab("N Children in total")+ ggtitle("Propofol True")+
  scale_x_discrete(labels = c("MBAOD","ODSC","Scaling","Adult Prior","OD"))+scale_y_continuous(breaks=seq(0, 45, 5),limits=c(0,45))+
  theme(axis.text = element_text(size=20,colour="black"),axis.title = element_text(size=23,face="bold"),
        axis.text.x = element_text(angle = 45, vjust = 0.65),plot.title=element_text(size=23,face="bold"))
p_box_true


p_box_large <- ggplot(data=df_box_large,aes(x=variable,y=value)) + geom_boxplot() +
  xlab(" ") + ylab("N Children in total")+ ggtitle("Propofol GFR")+
  scale_x_discrete(labels = c("MBAOD","ODSC","Scaling","Adult Prior","OD"))+scale_y_continuous(breaks=seq(0, 45, 5),limits=c(0,45))+
theme(axis.text = element_text(size=20,colour="black"),axis.title = element_text(size=23,face="bold"),
      axis.text.x = element_text(angle = 45, vjust = 0.65),plot.title=element_text(size=23,face="bold"))
p_box_large
##Plots done


#############################################
#########MBAOD and OD Group Selection Plots#####
############################################


eval_design <- function(result_dir=NULL,niterations=100,old=T){
  df_ids    <- data.frame(age_group=1:7, nID=rep(0,7))
  df_params <- data.frame(iteration=1:100, CL=rep(0,100),V=rep(0,100),TM50=rep(0,100),HILL=rep(0,100))
  power_res <- 0
  not_there <- 0
  if(is.null(result_dir)) {
    print("Please give a results directory")
  }else{
    for(i in 1:niterations){
      load(paste(result_dir,"/results_rep_",i,".Rdata",sep=""))
      for(j in 1:length(aod_res$final_design$x))
        
        
        if(aod_res$final_design$x[j] ==7 & aod_res$final_design$groupsize[j]==100){
          df_ids$nID[aod_res$final_design$x[j]] + df_ids$nID[aod_res$final_design$x[j]] +0  #due to it then being an adult prior group
        }else{
          df_ids$nID[aod_res$final_design$x[j]] <- df_ids$nID[aod_res$final_design$x[j]] + aod_res$final_design$groupsize[j]
        }
      
      df_params$CL[i] <- aod_res$est_summary$THETA_1[length(aod_res$est_summary$THETA_1)]
      df_params$V[i] <- aod_res$est_summary$THETA_2[length(aod_res$est_summary$THETA_2)]
      df_params$TM50[i] <- aod_res$est_summary$THETA_3[length(aod_res$est_summary$THETA_3)]
      df_params$HILL[i] <- aod_res$est_summary$THETA_4[length(aod_res$est_summary$THETA_4)]
      
      
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
  
  return(list(df_ids,df_params,power_res))
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


result_dir <- c("C:/Users/erist836/Desktop/MBAOD_project/merge")
MBAOD_true <- eval_design(result_dir = result_dir,old = F) 
df_LLT$nID_MBAOD <- MBAOD_true[[1]]$nID


result_dir=c("C:/Users/erist836/Desktop/MBAOD_project/PropofolT_power_eval_odsc_run_dir_3")
ODSC_true <- eval_design(result_dir = result_dir,old = F)   #39% power
df_LLT$nID_ODSC <-  ODSC_true[[1]]$nID    
ODSC_true[[3]]


result_dir="C:/Users/erist836/Desktop/MBAOD_project/propofolT_power_eval_scaling_run_dir_1"
scaling_true <- eval_design(result_dir = result_dir,old = F) #40% power
scaling_true[[3]]


result_dir="C:/Users/erist836/Desktop/MBAOD_project/propofolT_power_eval_adult_prior_run_dir_2"
adult_prior_true <- eval_design(result_dir = result_dir,old=F) #29% power
adult_prior_true[[3]]

result_dir="C:/Users/erist836/Desktop/MBAOD_project/Propofol_power_eval_OD_run_dir_1"
OD_true <- eval_design(result_dir = result_dir,old = F)#100% power
OD_true[[3]]

df_LLT$nID_Scaling[1:6] <- subset(scaling_power_CL, scaling_power_CL$Misspecification=="None")$nsub 
df_LLT$nID_Adult_Prior[1:6] <- nsub_adult_sd_est 
df_LLT$nID_OD <- c(0,12,0,0,12,12,0)#c(18,0,0,0,0,0,18) 


bias_df_true <- data.frame(DA = rep(c("MBAOD", "ODSC", "Scaling", "Adult Prior","OD" ), 1, each=100),
                           CL = rep(0,500),
                           V  = rep(0,500),
                           TM50 = rep(0,500),
                           Hill = rep(0,500)) 

bias_df_true$CL <- 100*(exp(c(MBAOD_true[[2]]$CL, ODSC_true[[2]]$CL, scaling_true[[2]]$CL, adult_prior_true[[2]]$CL, OD_true[[2]]$CL))-exp(1))/exp(1)

bias_df_true$V <- 100*(exp(c(MBAOD_true[[2]]$V, ODSC_true[[2]]$V, scaling_true[[2]]$V, adult_prior_true[[2]]$V, OD_true[[2]]$V))-exp(3))/exp(3) 

bias_df_true$TM50 <- 100*(exp(c(MBAOD_true[[2]]$TM50, ODSC_true[[2]]$TM50, scaling_true[[2]]$TM50, adult_prior_true[[2]]$TM50, OD_true[[2]]$TM50))-exp(3.6521))/exp(3.6521)

bias_df_true$Hill <- 100*(exp(c(MBAOD_true[[2]]$HILL, ODSC_true[[2]]$HILL, scaling_true[[2]]$HILL, adult_prior_true[[2]]$HILL, OD_true[[2]]$HILL)) -exp(1.5261))/exp(1.5261)
library(reshape) 
library(ggplot2)
i=4
titles<- c("MBAOD", "ODSC","Scaling","Adult Prior","OD")
df_box_tmp <- data.frame(CL = subset(bias_df_true,bias_df_true$DA==titles[i])$CL,
                         V = subset(bias_df_true,bias_df_true$DA==titles[i])$V,
                         TM50 = subset(bias_df_true,bias_df_true$DA==titles[i])$TM50,
                         Hill = subset(bias_df_true,bias_df_true$DA==titles[i])$Hill) 

df_box_m <- melt(df_box_tmp)

p_biasT <- ggplot(data=df_box_m,aes(x=variable,y=value)) + geom_boxplot() +
  xlab("Parameter") + ylab("% BIAS")+ ggtitle(paste("Propofol True - ", titles[i],sep=""))+
  scale_y_continuous(breaks=seq(-100, 150, 50),limits=c(-100,150))
#      theme(axis.text = element_text(size=20,colour="black"),axis.title = element_text(size=23,face="bold"),
#            axis.text.x = element_text(angle = 45, vjust = 0.65),plot.title=element_text(size=23,face="bold"))

p_biasT

#Misspecification values
result_dir <-c("C:/Users/erist836/Desktop/MBAOD_project/merge_gfr")
MBAOD_M <- eval_design(result_dir = result_dir,old = F)
df_LLM$nID_MBAOD <- MBAOD_M[[1]]$nID

result_dir=c("C:/Users/erist836/Desktop/MBAOD_project/gfrT_power_eval_odsc_run_dir_2")
ODSC_M <- eval_design(result_dir = result_dir, old = F)
df_LLM$nID_ODSC <- ODSC_M[[1]]$nID
ODSC_M[[3]]  #7% power

result_dir="C:/Users/erist836/Desktop/MBAOD_project/gfrT_power_eval_scaling_run_dir_NA"
scaling_M <- eval_design(result_dir = result_dir,old = F)
scaling_M[[3]] #33% power

result_dir="C:/Users/erist836/Desktop/MBAOD_project/gfrT_power_eval_adult_prior_run_dir_NA"
adult_prior_M <- eval_design(result_dir = result_dir,old=F)
adult_prior_M[[3]] #8% power

result_dir="C:/Users/erist836/Desktop/MBAOD_project/gfrT_power_eval_OD_run_dir_1"
OD_M <- eval_design(result_dir = result_dir,old=F)

OD_M[[3]]



df_LLM$nID_Scaling[1:6] <- subset(scaling_power_CL, scaling_power_CL$Misspecification=="Large")$nsub
df_LLM$nID_Adult_Prior[1:6] <- nsub_adult_sd_est
df_LLM$nID_OD <- c(0,12,0,0,12,12,0)#c(18,0,0,0,0,0,18) 



bias_df_m <- data.frame(DA = rep(c("MBAOD", "ODSC", "Scaling", "Adult Prior","OD" ), 1, each=100),
                        CL = rep(0,500),
                        V  = rep(0,500),
                        TM50 = rep(0,500),
                        Hill = rep(0,500))

bias_df_m$CL <- 100*(exp(c(MBAOD_M[[2]]$CL, ODSC_M[[2]]$CL, scaling_M[[2]]$CL, adult_prior_M[[2]]$CL, OD_M[[2]]$CL))-exp(1))/exp(1)

bias_df_m$V <- 100*(exp(c(MBAOD_M[[2]]$V, ODSC_M[[2]]$V, scaling_M[[2]]$V, adult_prior_M[[2]]$V, OD_M[[2]]$V))-exp(3))/exp(3) 

bias_df_m$TM50 <- 100*(exp(c(MBAOD_M[[2]]$TM50, ODSC_M[[2]]$TM50, scaling_M[[2]]$TM50, adult_prior_M[[2]]$TM50, OD_M[[2]]$TM50))-exp(3.6521))/exp(3.6521)

bias_df_m$Hill <- 100*(exp(c(MBAOD_M[[2]]$HILL, ODSC_M[[2]]$HILL, scaling_M[[2]]$HILL, adult_prior_M[[2]]$HILL, OD_M[[2]]$HILL)) -exp(1.5261))/exp(1.5261)
library(reshape)
i=4
titles<- c("MBAOD", "ODSC","Scaling","Adult Prior","OD")
df_box_tmp <- data.frame(CL = subset(bias_df_m,bias_df_m$DA==titles[i])$CL,
                         V = subset(bias_df_m,bias_df_m$DA==titles[i])$V,
                         TM50 = subset(bias_df_m,bias_df_m$DA==titles[i])$TM50,
                         Hill = subset(bias_df_true,bias_df_true$DA==titles[i])$Hill)

df_box_m <- melt(df_box_tmp)

p_biasm <- ggplot(data=df_box_m,aes(x=variable,y=value)) + geom_boxplot() +
  xlab("Parameter") + ylab("% BIAS")+ ggtitle(paste("Propofol GFR - ", titles[i],sep=""))+
  scale_y_continuous(breaks=seq(-100, 150, 50),limits=c(-100,150))
#      theme(axis.text = element_text(size=20,colour="black"),axis.title = element_text(size=23,face="bold"),
#            axis.text.x = element_text(angle = 45, vjust = 0.65),plot.title=element_text(size=23,face="bold"))

p_biasm


#Relative number of individuals
df_LLT_rel <- df_LLT
df_LLT_rel$nID_MBAOD       <- df_LLT$nID_MBAOD/sum(df_LLT$nID_MBAOD)
df_LLT_rel$nID_ODSC        <- df_LLT$nID_ODSC/sum(df_LLT$nID_ODSC)
df_LLT_rel$nID_Scaling     <- df_LLT$nID_Scaling/sum(df_LLT$nID_Scaling)
df_LLT_rel$nID_Adult_Prior <- df_LLT$nID_Adult_Prior/sum(df_LLT$nID_Adult_Prior)
df_LLT_rel$nID_OD          <- df_LLT$nID_OD/sum(df_LLT$nID_OD)

df_LLM_rel <- df_LLM
df_LLM_rel$nID_MBAOD       <- df_LLM$nID_MBAOD/sum(df_LLM$nID_MBAOD)
df_LLM_rel$nID_ODSC        <- df_LLM$nID_ODSC/sum(df_LLM$nID_ODSC)
df_LLM_rel$nID_Scaling     <- df_LLM$nID_Scaling/sum(df_LLM$nID_Scaling)
df_LLM_rel$nID_Adult_Prior <- df_LLM$nID_Adult_Prior/sum(df_LLM$nID_Adult_Prior)
df_LLM_rel$nID_OD          <- df_LLM$nID_OD/sum(df_LLM$nID_OD)

df_LLT_rel[,2:6] <- df_LLT_rel[,2:6] *100 
df_LLM_rel[,2:6] <- df_LLM_rel[,2:6] *100 



i=5
df_box_tmp <- data.frame(ag= 1:7, rel_nid = df_LLT_rel[,i+1])

titles<- c("MBAOD", "ODSC","Scaling","Adult Prior","OD")
p_LLT <- ggplot(data=df_box_tmp,aes(x=ag,y=rel_nid)) + geom_bar(stat="identity") +
  xlab("Age Group") + ylab("% of individuals")+ ggtitle(paste("Propofol True - ", titles[i],sep=""))+
  scale_x_continuous(breaks=seq(1:7)) + scale_y_continuous(breaks=seq(0, 70, 5),limits=c(0,75))
#      theme(axis.text = element_text(size=20,colour="black"),axis.title = element_text(size=23,face="bold"),
#            axis.text.x = element_text(angle = 45, vjust = 0.65),plot.title=element_text(size=23,face="bold"))

p_LLT

df_box_tmp <- data.frame(ag= 1:7, rel_nid = df_LLM_rel[,i+1])
p_LLM <- ggplot(data=df_box_tmp,aes(x=ag,y=rel_nid)) + geom_bar(stat="identity") +
  xlab("Age Group") + ylab("% of individuals")+ ggtitle(paste("Propofol GFR - ", titles[i],sep=""))+
  scale_x_continuous(breaks=seq(1:7)) + scale_y_continuous(breaks=seq(0, 70, 5),limits=c(0,75))
#      theme(axis.text = element_text(size=20,colour="black"),axis.title = element_text(size=23,face="bold"),
#            axis.text.x = element_text(angle = 45, vjust = 0.65),plot.title=element_text(size=23,face="bold"))

p_LLM


   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
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


library(reshape)
library(ggplot2)
source("get_weight.R")
source("stop_critX_2.R")


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
    
    d9_CI_CL[k+21*(i-1),5] <- size_mat_scaling_l(final_estimates[i,2:4],age,wt)[1]
    d9_CI_V[k+21*(i-1),5]  <- size_mat_scaling_l(final_estimates[i,2:4],age,wt)[2]
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
  df_CL_coverage[i,6:7] <- quantile((tmp2$np),probs=c(0.025,0.975),na.rm=T)/gm_mean((tmp2$np))
  
  
  tmp3 <- subset(d9_CI_V,pma==AGE_WT$PMA[i])
  df_V_coverage[i,2:3] <- c(min(tmp3$d9_lower,na.rm=T),max(tmp3$d9_lower,na.rm=T))#quantile((tmp3$d9_lower),probs=c(0.025,0.975),na.rm=T)
  df_V_coverage[i,4:5] <- c(min(tmp3$d9_higher,na.rm=T),max(tmp3$d9_higher,na.rm=T))#quantile((tmp3$d9_higher),probs=c(0.025,0.975),na.rm=T)
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
  xlab("PMA") + ylab("% of geometric mean")+ggtitle("LLM: CI coverage of CLp")+
  theme(axis.text = element_text(size=20,colour="black"),axis.title = element_text(size=23,face="bold"),
        axis.text.x = element_text(angle = 60,vjust = 0.65),plot.title=element_text(size=23,face="bold"))+
  scale_x_discrete(breaks=1:21,labels=(AGE_WT$PMA))+scale_y_continuous(breaks=seq(50,150,10),limits=c(50,150))
coverage_CI_CL

coverage_CI_V <- ggplot(data=(df_V_coverage_m),aes(x=x_step,y=lower)) +
  geom_errorbar(aes(ymin=lower*100, ymax=higher*100,colour=Legend)) +
  xlab("PMA") + ylab("% of geometric mean")+ggtitle("LLM: CI coverage of Vp")+
  theme(axis.text = element_text(size=20,colour="black"),axis.title = element_text(size=23,face="bold"),
        axis.text.x = element_text(angle = 60,vjust = 0.65),plot.title=element_text(size=23,face="bold"))+
  scale_x_discrete(breaks=1:21,labels=(AGE_WT$PMA))+scale_y_continuous(breaks=seq(90,110,2),limits=c(90,110))
coverage_CI_V

   
   
   
   
   
