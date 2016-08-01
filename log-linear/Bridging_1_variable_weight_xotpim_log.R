
##Simulated bridging study of IV admin. drug X with true 1 comp kinetics from
# 25 y.o. adults with weight 70 to infants 3 monts -2 years old and
## simulated values (the truth) in this case will be c(TM50=100 (WEEKS),hill=1)
##Code examples for how to run MBAOD sim and regular OD with small misspecifications (TM50=75)
##Are presented below. TO run the MBAOD and OD with different or no misspecification in the initial guess
##Of the TM50 parameter, simply change the initial parameter guess where it is marked in the code.




# remove things from the global environment
rm(list=ls())

library(PopED)
library(mvtnorm)
library(xpose4)

# load the MBAOD package. Change to your relevant load call
devtools::load_all("C:/Users/erist836/Documents/GitHub/MBAOD/R")

# load the PopED model file
source("PopED_files/poped.mod.PK.1.comp.maturation_Xlog.R")
# load the weight estimator file
source("get_weight.R")
source("stop_critX_power.R")


############################################RUN THIS SECTION FOR MBAOD WITH SMALL Misspecification############################################
#Adults
step_1=list(
  design = list(
    groupsize = 100,
    x = c(age_group=c(7)),
    xt = c(0.1, 2, 6, 12, 24)
  ),
  optimize=NULL,
  simulate=list(target="NONMEM", model="NONMEM_files/sim_logX.mod",
                data=list(dosing = list(list(AMT=1000,Time=0)),
                          manipulation = list(expression(PMA <- age_group_2_PMA(ID,age_group)),
                                              expression(WT <- get_weight(ID,PMA,probs = c(0.1,0.9))),
                                              expression(AMT <- AMT*WT/70)
                                             )
                          )
                
  ),
  estimate=list(target="NONMEM", model="NONMEM_files/est_red_logX.mod")
)

x.space <- cell(1,1)
x.space[1,1] <- list(c(6,7))


#First cohort with children
step_2 = list(
  design = list(
    groupsize = 9,
    m=1,
    x = c(age_group=c(6)),
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
                  bpop=c(TM50=100), # initial parameter, guess true value = 100 CHANGE HERE IF YOU WANT A DIFFERENT GUESS
                  manipulation=NULL 
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
  simulate=list(target="NONMEM", model="NONMEM_files/sim_logX.mod",
                data=list(
                  dosing = list(list(AMT=1000,Time=0)),
                  manipulation = list(expression(PMA <- age_group_2_PMA(ID,age_group)),
                                      expression(WT <- get_weight(ID,PMA,probs = c(0.1,0.9))),
                                                 expression(AMT <- AMT*WT/70)
                                      )
                         )                      
                
  ),
  estimate=list(target="NONMEM", model="NONMEM_files/est_full_log_true.mod")
)


step_3 <- step_2
step_3$optimize$parameters <- NULL
step_3$design$groupsize <- 2


results_mbaod_small <- mbaod_simulate(cohorts=list(step_1,step_2,step_3), # anything after step_3 is the same as step_3
                              ncohorts=100, # number of steps or cohorts in one AOD
                              rep=100, #number of times to repeat the MBAOD simulation 
                              name="bridging_maturation_model_xoptim_restricted_group_small", 
                              description="25 steps, 1st step one group, steps 2-10 have 1 added group per step",
                              seednr=321,stop_crit_fun =stop_critX,zip_directories=F)



###########################Run this section for a regular OD with small misspecification#####################################

# remove things from the global environment
rm(list=ls())

#Adults
step_1=list(
  design = list(
    groupsize = 100,
    m=1,
    x = c(age_group=c(7)),
    xt = c(0.1, 2, 6, 12, 24)
  ),
  optimize=NULL,
  simulate=list(target="NONMEM", model="C:/Users/Eric Stromberg/Desktop/MBAOD_project/NONMEM_files/sim_logX.mod",
                data=list(dosing = list(list(AMT=1000,Time=0)),
                          manipulation = list(expression(PMA <- age_group_2_PMA(ID,age_group)),
                                              expression(WT <- get_weight(ID,PMA,probs = c(0.1,0.9))),
                                              expression(AMT <- AMT*WT/70)
                          )
                )
                
  ),
  estimate=list(target="NONMEM", model="C:/Users/Eric Stromberg/Desktop/MBAOD_project/NONMEM_files/est_red_logX.mod")
)

##Initial Design space for the age groups
x.space <- cell(1,1)
x.space[1,1] <- list(c(6,7))
###

#Children 6 months, Children 1.5 years old
step_2 = list(
  design = list(
    groupsize = 9,
    m=1,
    x = c(age_group=c(6)),
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
                  bpop=c(TM50=75),# initial parameter, guess true value = 100 CHANGE HERE IF YOU WANT A DIFFERENT GUESS
                  manipulation=NULL 
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
                  compute_inv=T,    
                  iFIMCalculationType=0
                )
  ),
  simulate=list(target="NONMEM", model="C:/Users/Eric Stromberg/Desktop/MBAOD_project/NONMEM_files/sim_logX.mod",
                data=list(
                  dosing = list(list(AMT=1000,Time=0)),
                  manipulation = list(expression(PMA <- age_group_2_PMA(ID,age_group)),
                                      expression(WT <- get_weight(ID,PMA,probs = c(0.1,0.9))),
                                      expression(AMT <- AMT*WT/70)
                  )
                )                      
                
  ),
  estimate=list(target="NONMEM", model="C:/Users/Eric Stromberg/Desktop/MBAOD_project/NONMEM_files/est_full_log_OD.mod")
)


step_3 <- step_2
step_3$optimize$parameters <- NULL
step_3$design$groupsize <- 2



results_all_OD_small <- mbaod_simulate(cohorts=list(step_1,step_2,step_3), # anything after step_3 is the same as step_3
                              ncohorts=100, # number of steps or cohorts in one AOD
                              rep=100, #number of times to repeat the MBAOD simulation 
                              name="bridging_maturation_model_Xlog_small_OD_simcrit", 
                              description="10 steps, 1st step one group, steps 2-10 have 2 groups per step",
                              seednr=321, stop_crit_fun = stop_crit_PopED)


##################END Regular OD Small Misspecification##############################################################










