#Power evaluation and bias calc for alternate design approaches.
#Partial results presented at Population Group Approach Europe 2015
#(LINK HERE)


#set path to the directory where this file is located
setwd("C:/Users/erist836/Documents/GitHub/mbaod_sim/non-linear")

# load the MBAOD package. Change to your relevant load call
devtools::load_all("C:/Users/erist836/Documents/GitHub/MBAOD/R")


# remove things from the global environment
rm(list=ls())

library(PopED)
library(mvtnorm)
library(xpose4)

# load the PopED model file
source("PopED_files/poped.mod.PK.1.comp.maturation_Xlog.R")
# load the weight estimator file
source("get_weight.R")
source("stop_critX_power.R")
source("mbaod_simulate.R") ##Change to your dir
###########LLT ODSC###############

#PAGE 2015 results gave 72%

###########LLT Scaling###########
#PAGE 2015 results gave 100%

step_1_LLT_scaling=list(
  design = list(
    groupsize = c(100,6,8,7,6,6,6),
    m=7,
    x = t(rbind(age_group=c(7,6,5,4,3,2,1))),
    xt = c(0.1, 2, 6, 12, 24)
  ),
  optimize=NULL,
  simulate=list(target="NONMEM", model="NONMEM_files/sim_logX.mod",
                data=list(dosing = list(list(AMT=1000,Time=0)),
                          manipulation = list(expression(PMA <- age_group_2_PMA(ID,age_group)),
                                              expression(WT <- get_weight(ID,PMA,probs = c(0.1,0.9))),
                                              expression(AMT <- AMT*(WT/70))
                          )
                )
                
  ),
  estimate=list(target="NONMEM", model="NONMEM_files/est_full_log_true2.mod")
)


results_LLM_scaling <- mbaod_simulate(cohorts=list(step_1_LLT_scaling), # anything after step_3 is the same as step_3
                                      ncohorts=1, # number of steps or cohorts in one AOD
                                      rep=100, #number of times to repeat the MBAOD simulation 
                                      name="LLT_power_eval_scaling", 
                                      description="25 steps, 1st step one group, steps 2-10 have 1 added group per step",
                                      seednr=123,stop_crit_fun =stop_critX_power)

###########LLT Adult Prior#######
#PAGE 2015 results gave 100%

step_1_LLT_adult_prior=list(
  design = list(
    groupsize = c(100,6,6,6,6,6,6),
    m=7,
    x = t(rbind(age_group=c(7,6,5,4,3,2,1))),
    xt = c(0.1, 2, 6, 12, 24)
  ),
  optimize=NULL,
  simulate=list(target="NONMEM", model="NONMEM_files/sim_logX.mod",
                data=list(dosing = list(list(AMT=1000,Time=0)),
                          manipulation = list(expression(PMA <- age_group_2_PMA(ID,age_group)),
                                              expression(WT <- get_weight(ID,PMA,probs = c(0.1,0.9))),
                                              expression(AMT <- AMT*(WT/70))
                          )
                )
                
  ),
  estimate=list(target="NONMEM", model="NONMEM_files/est_full_log_true2.mod")
)


results_LLT_adult_prior <- mbaod_simulate(cohorts=list(step_1_LLT_adult_prior), # anything after step_3 is the same as step_3
                                          ncohorts=1, # number of steps or cohorts in one AOD
                                          rep=100, #number of times to repeat the MBAOD simulation 
                                          name="LLT_power_eval_adult_prior", 
                                          description="25 steps, 1st step one group, steps 2-10 have 1 added group per step",
                                          seednr=123,stop_crit_fun =stop_critX_power)





###########LLM ODSC#############
#PAGE 2015 results gave 70%

###########LLM Scaling#########
#PAGE 2015 results gave 100%
#6 6 6 8 8 6 0
#> df_LLT$nID_Scaling
#[1] 6 6 6 7 8 6 0

step_1_LLM_scaling=list(
  design = list(
    groupsize = c(100,6,8,8,6,6,6),
    m=7,
    x = t(rbind(age_group=c(7,6,5,4,3,2,1))),
    xt = c(0.1, 2, 6, 12, 24)
  ),
  optimize=NULL,
  simulate=list(target="NONMEM", model="NONMEM_files/sim_logX.mod",
                data=list(dosing = list(list(AMT=1000,Time=0)),
                          manipulation = list(expression(PMA <- age_group_2_PMA(ID,age_group)),
                                              expression(WT <- get_weight(ID,PMA,probs = c(0.1,0.9))),
                                              expression(AMT <- AMT*(WT/70))
                          )
                )
                
  ),
  estimate=list(target="NONMEM", model="NONMEM_files/est_full_log_large2.mod")
)


results_LLM_scaling <- mbaod_simulate(cohorts=list(step_1_LLM_scaling), # anything after step_3 is the same as step_3
                              ncohorts=1, # number of steps or cohorts in one AOD
                              rep=100, #number of times to repeat the MBAOD simulation 
                              name="LLM_power_eval_scaling", 
                              description="25 steps, 1st step one group, steps 2-10 have 1 added group per step",
                              seednr=123,stop_crit_fun =stop_critX_power)



##########LLM Adult Prior#####
#PAGE 2015 results gave 100%
step_1_LLM_adult_prior=list(
  design = list(
    groupsize = c(100,6,6,6,6,6,6),
    m=7,
    x = t(rbind(age_group=c(7,6,5,4,3,2,1))),
    xt = c(0.1, 2, 6, 12, 24)
  ),
  optimize=NULL,
  simulate=list(target="NONMEM", model="NONMEM_files/sim_logX.mod",
                data=list(dosing = list(list(AMT=1000,Time=0)),
                          manipulation = list(expression(PMA <- age_group_2_PMA(ID,age_group)),
                                              expression(WT <- get_weight(ID,PMA,probs = c(0.1,0.9))),
                                              expression(AMT <- AMT*(WT/70))
                          )
                )
                
  ),
  estimate=list(target="NONMEM", model="NONMEM_files/est_full_log_large2.mod")
)


results_LLM_adult_prior <- mbaod_simulate(cohorts=list(step_1_LLM_adult_prior), # anything after step_3 is the same as step_3
                                      ncohorts=1, # number of steps or cohorts in one AOD
                                      rep=31, #number of times to repeat the MBAOD simulation 
                                      name="LLM_power_eval_adult_prior", 
                                      description="25 steps, 1st step one group, steps 2-10 have 1 added group per step",
                                      seednr=1234,stop_crit_fun =stop_critX_power)







#LLM OD
step_1_LLM=list(
  design = list(
    groupsize = c(18,18),
    m=2,
    x = t(rbind(age_group=c(7,1))),
    xt = c(0.1, 2, 6, 12, 24)
  ),
  optimize=NULL,
  simulate=list(target="NONMEM", model="NONMEM_files/sim_logX.mod",
                data=list(dosing = list(list(AMT=1000,Time=0)),
                          manipulation = list(expression(PMA <- age_group_2_PMA(ID,age_group)),
                                              expression(WT <- get_weight(ID,PMA,probs = c(0.1,0.9))),
                                              expression(AMT <- AMT*(WT/70))
                          )
                )
                
  ),
  estimate=list(target="NONMEM", model="NONMEM_files/est_full_log_large2.mod")
)


results_LLM <- mbaod_simulate(cohorts=list(step_1_LLM), # anything after step_3 is the same as step_3
                                      ncohorts=1, # number of steps or cohorts in one AOD
                                      rep=100, #number of times to repeat the MBAOD simulation 
                                      name="LLM_power_eval", 
                                      description="25 steps, 1st step one group, steps 2-10 have 1 added group per step",
                                      seednr=123,stop_crit_fun =stop_critX_power)




#Propofol True ODSC
step_1_propofolT_odsc=list(
  design = list(
    groupsize = c(100,2,2,2,2,2),
    m=6,
    x = t(rbind(age_group=c(7,6,5,4,3,2))),
    xt = c(0.1, 2, 6, 12, 24)
  ),
  optimize=NULL,
  simulate=list(target="NONMEM", model="NONMEM_files/sim_propofol.mod",
                data=list(dosing = list(list(AMT=1000,Time=0)),
                          manipulation = list(expression(PMA <- age_group_2_PMA(ID,age_group)),
                                              expression(WT <- get_weight(ID,PMA,probs = c(0.1,0.9))),
                                              expression(AMT <- AMT*(WT/70))
                          )
                )
                
  ),
  estimate=list(target="NONMEM", model="NONMEM_files/est_full_propofol.mod")
)


results_propofolT_odsc <- mbaod_simulate(cohorts=list(step_1_propofolT_odsc), # anything after step_3 is the same as step_3
                                                ncohorts=1, # number of steps or cohorts in one AOD
                                                rep=100, #number of times to repeat the MBAOD simulation 
                                                name="PropofolT_power_eval_odsc", 
                                                description="25 steps, 1st step one group, steps 2-10 have 1 added group per step",
                                                seednr=123,stop_crit_fun =stop_critX_power)


step_1_gfrT_odsc=list(
  design = list(
    groupsize = c(100,2,2,2,2,2),
    m=6,
    x = t(rbind(age_group=c(7,6,5,4,3,2))),
    xt = c(0.1, 2, 6, 12, 24)
  ),
  optimize=NULL,
  simulate=list(target="NONMEM", model="NONMEM_files/sim_propofol.mod",
                data=list(dosing = list(list(AMT=1000,Time=0)),
                          manipulation = list(expression(PMA <- age_group_2_PMA(ID,age_group)),
                                              expression(WT <- get_weight(ID,PMA,probs = c(0.1,0.9))),
                                              expression(AMT <- AMT*(WT/70))
                          )
                )
                
  ),
  estimate=list(target="NONMEM", model="NONMEM_files/est_full_gfr.mod")
)


results_gfrT_odsc <- mbaod_simulate(cohorts=list(step_1_gfrT_odsc), # anything after step_3 is the same as step_3
                                         ncohorts=1, # number of steps or cohorts in one AOD
                                         rep=100, #number of times to repeat the MBAOD simulation 
                                         name="gfrT_power_eval_odsc", 
                                         description="25 steps, 1st step one group, steps 2-10 have 1 added group per step",
                                         seednr=1234,stop_crit_fun =stop_critX_power)

#Propofol True Scaling
step_1_propofolT_scaling=list(
  design = list(
    groupsize = c(100,6,7,6,6,6,6),
    m=7,
    x = t(rbind(age_group=c(7,6,5,4,3,2,1))),
    xt = c(0.1, 2, 6, 12, 24)
  ),
  optimize=NULL,
  simulate=list(target="NONMEM", model="NONMEM_files/sim_propofol.mod",
                data=list(dosing = list(list(AMT=1000,Time=0)),
                          manipulation = list(expression(PMA <- age_group_2_PMA(ID,age_group)),
                                              expression(WT <- get_weight(ID,PMA,probs = c(0.1,0.9))),
                                              expression(AMT <- AMT*(WT/70))
                          )
                )
                
  ),
  estimate=list(target="NONMEM", model="NONMEM_files/est_full_propofol.mod")
)


results_propofolT_scaling <- mbaod_simulate(cohorts=list(step_1_propofolT_scaling), # anything after step_3 is the same as step_3
                                          ncohorts=1, # number of steps or cohorts in one AOD
                                          rep=2, #number of times to repeat the MBAOD simulation 
                                          name="propofolT_power_eval_scaling", 
                                          description="25 steps, 1st step one group, steps 2-10 have 1 added group per step",
                                          seednr=1234,stop_crit_fun =stop_critX_power)





#Propofol True Adult Prior
step_1_propofolT_adult_prior=list(
  design = list(
    groupsize = c(100,6,6,6,6,6,6),
    m=7,
    x = t(rbind(age_group=c(7,6,5,4,3,2,1))),
    xt = c(0.1, 2, 6, 12, 24)
  ),
  optimize=NULL,
  simulate=list(target="NONMEM", model="NONMEM_files/sim_propofol.mod",
                data=list(dosing = list(list(AMT=1000,Time=0)),
                          manipulation = list(expression(PMA <- age_group_2_PMA(ID,age_group)),
                                              expression(WT <- get_weight(ID,PMA,probs = c(0.1,0.9))),
                                              expression(AMT <- AMT*(WT/70))
                          )
                )
                
  ),
  estimate=list(target="NONMEM", model="NONMEM_files/est_full_propofol.mod")
)


results_propofolT_adult_prior <- mbaod_simulate(cohorts=list(step_1_propofolT_adult_prior), # anything after step_3 is the same as step_3
                                          ncohorts=1, # number of steps or cohorts in one AOD
                                          rep=100, #number of times to repeat the MBAOD simulation 
                                          name="PropofolT_power_eval_adult_prior", 
                                          description="25 steps, 1st step one group, steps 2-10 have 1 added group per step",
                                          seednr=1234,stop_crit_fun =stop_critX_power)





#Propofol True OD


#Propofol GFR ODSC
step_1_gfrT_odsc=list(
  design = list(
    groupsize = c(100,2,2,2,2,2),
    m=6,
    x = t(rbind(age_group=c(7,6,5,4,3,2))),
    xt = c(0.1, 2, 6, 12, 24)
  ),
  optimize=NULL,
  simulate=list(target="NONMEM", model="NONMEM_files/sim_propofol.mod",
                data=list(dosing = list(list(AMT=1000,Time=0)),
                          manipulation = list(expression(PMA <- age_group_2_PMA(ID,age_group)),
                                              expression(WT <- get_weight(ID,PMA,probs = c(0.1,0.9))),
                                              expression(AMT <- AMT*(WT/70))
                          )
                )
                
  ),
  estimate=list(target="NONMEM", model="NONMEM_files/est_full_gfr.mod")
)


results_gfrT_odsc <- mbaod_simulate(cohorts=list(step_1_gfrT_odsc), # anything after step_3 is the same as step_3
                                    ncohorts=1, # number of steps or cohorts in one AOD
                                    rep=100, #number of times to repeat the MBAOD simulation 
                                    name="gfrT_power_eval_odsc", 
                                    description="25 steps, 1st step one group, steps 2-10 have 1 added group per step",
                                    seednr=1234,stop_crit_fun =stop_critX_power)
#Propofol GFR Scaling
step_1_gfrT_scaling=list(
  design = list(
    groupsize = c(100,6,7,6,6,6,6),
    m=7,
    x = t(rbind(age_group=c(7,6,5,4,3,2,1))),
    xt = c(0.1, 2, 6, 12, 24)
  ),
  optimize=NULL,
  simulate=list(target="NONMEM", model="NONMEM_files/sim_propofol.mod",
                data=list(dosing = list(list(AMT=1000,Time=0)),
                          manipulation = list(expression(PMA <- age_group_2_PMA(ID,age_group)),
                                              expression(WT <- get_weight(ID,PMA,probs = c(0.1,0.9))),
                                              expression(AMT <- AMT*(WT/70))
                          )
                )
                
  ),
  estimate=list(target="NONMEM", model="NONMEM_files/est_full_gfr.mod")
)


results_gfrT_scaling <- mbaod_simulate(cohorts=list(step_1_gfrT_scaling), # anything after step_3 is the same as step_3
                                                ncohorts=1, # number of steps or cohorts in one AOD
                                                rep=25, #number of times to repeat the MBAOD simulation 
                                                name="gfrT_power_eval_scaling2", 
                                                description="25 steps, 1st step one group, steps 2-10 have 1 added group per step",
                                                seednr=123,stop_crit_fun =stop_critX_power)

#Propofol GFR Adult Prior

step_1_gfrT_adult_prior=list(
  design = list(
    groupsize = c(100,6,6,6,6,6,6),
    m=7,
    x = t(rbind(age_group=c(7,6,5,4,3,2,1))),
    xt = c(0.1, 2, 6, 12, 24)
  ),
  optimize=NULL,
  simulate=list(target="NONMEM", model="NONMEM_files/sim_propofol.mod",
                data=list(dosing = list(list(AMT=1000,Time=0)),
                          manipulation = list(expression(PMA <- age_group_2_PMA(ID,age_group)),
                                              expression(WT <- get_weight(ID,PMA,probs = c(0.1,0.9))),
                                              expression(AMT <- AMT*(WT/70))
                          )
                )
                
  ),
  estimate=list(target="NONMEM", model="NONMEM_files/est_full_gfr.mod")
)


results_gfrT_adult_prior <- mbaod_simulate(cohorts=list(step_1_gfrT_adult_prior), # anything after step_3 is the same as step_3
                                                ncohorts=1, # number of steps or cohorts in one AOD
                                                rep=100, #number of times to repeat the MBAOD simulation 
                                                name="gfrT_power_eval_adult_prior2", 
                                                description="25 steps, 1st step one group, steps 2-10 have 1 added group per step",
                                                seednr=1234,stop_crit_fun =stop_critX_power)



#Propofol GFR OD

step_1_gfrT_OD=list(
  design = list(
    groupsize = c(100,12,12,12),
    m=4,
    x = t(rbind(age_group=c(7,6,3,2))),
    xt = c(0.1, 2, 6, 12, 24)
  ),
  optimize=NULL,
  simulate=list(target="NONMEM", model="NONMEM_files/sim_propofol.mod",
                data=list(dosing = list(list(AMT=1000,Time=0)),
                          manipulation = list(expression(PMA <- age_group_2_PMA(ID,age_group)),
                                              expression(WT <- get_weight(ID,PMA,probs = c(0.1,0.9))),
                                              expression(AMT <- AMT*(WT/70))
                          )
                )
                
  ),
  estimate=list(target="NONMEM", model="NONMEM_files/est_full_gfr.mod")
)


results_gfrT_OD <- mbaod_simulate(cohorts=list(step_1_gfrT_OD), # anything after step_3 is the same as step_3
                                           ncohorts=1, # number of steps or cohorts in one AOD
                                           rep=30, #number of times to repeat the MBAOD simulation 
                                           name="gfrT_power_eval_OD2", 
                                           description="25 steps, 1st step one group, steps 2-10 have 1 added group per step",
                                           seednr=12345,stop_crit_fun =stop_critX_power)


step_1_T_OD=list(
  design = list(
    groupsize = c(100,12,12,12),
    m=4,
    x = t(rbind(age_group=c(7,6,3,2))),
    xt = c(0.1, 2, 6, 12, 24)
  ),
  optimize=NULL,
  simulate=list(target="NONMEM", model="NONMEM_files/sim_propofol.mod",
                data=list(dosing = list(list(AMT=1000,Time=0)),
                          manipulation = list(expression(PMA <- age_group_2_PMA(ID,age_group)),
                                              expression(WT <- get_weight(ID,PMA,probs = c(0.1,0.9))),
                                              expression(AMT <- AMT*(WT/70))
                          )
                )
                
  ),
  estimate=list(target="NONMEM", model="NONMEM_files/est_full_propofol.mod")
)


results_T_OD <- mbaod_simulate(cohorts=list(step_1_T_OD), # anything after step_3 is the same as step_3
                                  ncohorts=1, # number of steps or cohorts in one AOD
                                  rep=100, #number of times to repeat the MBAOD simulation 
                                  name="Propofol_power_eval_OD2", 
                                  description="25 steps, 1st step one group, steps 2-10 have 1 added group per step",
                                  seednr=1234,stop_crit_fun =stop_critX_power)









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

############################################RUN THIS SECTION FOR MBAOD WITH SMALL Misspecification############################################
#Adults
step_1=list(
  design = list(
    groupsize = 100,
    x = c(age_group=c(7)),
    xt = c(0.1, 2, 6, 12, 24)
  ),
  optimize=NULL,
  simulate=list(target="NONMEM", model="NONMEM_files/sim_propofol.mod",
                data=list(dosing = list(list(AMT=1000,Time=0)),
                          manipulation = list(expression(PMA <- age_group_2_PMA(ID,age_group)),
                                              expression(WT <- get_weight(ID,PMA,probs = c(0.1,0.9))),
                                              expression(AMT <- AMT*(WT/70))
                          )
                )
                
  ),
  estimate=list(target="NONMEM", model="NONMEM_files/est_red_propofol.mod")
)

##Initial Design space for the age groups. subadults 12-18 y.o. and adults.
##If only adults and subadults are allowed, the information about TM50 is too sparse.
x.space <- cell(5,1)
x.space[,] <- list(c(2,3,4,5,6,7))
###

tm50 <- 3.651


#First cohort with children
step_2 = list(
  design = list(
    groupsize = 2,
    m=5,
    x = t(rbind(age_group=c(6,5,4,3,2))),
    xt = c(0.1, 2, 6, 12, 24)
  ),
  optimize=list(target="poped_R",
                model = list(
                  ff_file="PK.1.comp.maturation.ff",
                  fError_file="feps.add.prop",
                  fg_file="PK.1.comp.maturation.fg"
                ),
                design_space=list(minxt=0,
                                  maxxt=24,
                                  x_space = x.space
                ),
                parameters=list(
                  bpop=c(TM50=3.651, HILL=1.5261), # initial parameter guess not coming from previous step
                  manipulation=NULL # manipulation of initial parameters
                ),
                settings.db=list(
                  notfixed_sigma=c(1,0)
                ),
                settings.opt=list(
                  opt_xt=F,
                  opt_a=F,
                  opt_x=T,
                  bUseRandomSearch= 0,
                  bUseStochasticGradient = 0,
                  bUseBFGSMinimizer = 0,
                  bUseLineSearch = 1,
                  bUseExchangeAlgorithm=0,
                  EACriteria = 1,
                  compute_inv=T
                )
  ),
  simulate=list(target="NONMEM", model="NONMEM_files/sim_propofol.mod",
                data=list(
                  dosing = list(list(AMT=1000,Time=0)),
                  manipulation = list(expression(PMA <- age_group_2_PMA(ID,age_group)),
                                      expression(WT <- get_weight(ID,PMA,probs = c(0.1,0.9))),
                                      expression(AMT <- AMT*(WT/70))
                  )
                )                      
                
  ),
  estimate=list(target="NONMEM", model="NONMEM_files/est_full_propofol.mod")
)

x.space <- cell(2,1)
x.space[,] <- list(c(6,7))
step_3 <- step_2
step_3$optimize$parameters <- NULL
step_3$design$groupsize <- 1
step_3$design$x <- t(rbind(age_group=c(6,6)))
step_3$design$m <- 2
step_3$optimize$design_space$x_space<- x.space


results_mbaod_small <- mbaod_simulate(cohorts=list(step_1,step_2,step_3), # anything after step_3 is the same as step_3
                                      ncohorts=60, # number of steps or cohorts in one AOD
                                      rep=100, #number of times to repeat the MBAOD simulation 
                                      name="propofol_run_1234", lower_limit=0.6,higher_limit=1.4,
                                      description="25 steps, 1st step one group, steps 2-10 have 1 added group per step",
                                      seednr=1234, stop_crit_fun =stop_critX_2,run_commands="-retries=5", ci_corr=0,option=3)


