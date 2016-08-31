
#set path to the directory where this file is located
setwd("C:/Users/erist836/Desktop/MBAOD2")

# remove things from the global environment
rm(list=ls())

library(PopED)
library(mvtnorm)
library(xpose4)

# load the MBAOD package, change to your relevant load call
devtools::load_all("C:/Users/erist836/Documents/GitHub/MBAOD/R")

# load the PopED model file
source("PopED_files/PKPD_example_1.R")

############################################RUN THIS SECTION FOR MBAOD WITH SMALL Misspecification############################################

step_1=list(
  design = list(
    groupsize = 4,
    m=4,
    a = list(c(AMT=20),c(AMT=40),c(AMT=80),c(AMT=160)),
    xt = c(0.5,2,3,6,24,36,72,120)
  ),
  optimize=NULL,
  simulate=list(target="NONMEM", model="NONMEM_files/ex1_sim.mod",
                data=list(dosing = list(list(AMT0=1,Time=0)),
                          manipulation = list(expression(AMT <- AMT0*AMT))
                )
                
  ),
  estimate=list(target="NONMEM", model="NONMEM_files/ex1_est_full.mod")
)




#First cohort with dose optimization
step_2 = list(
  design = list(
    groupsize = 1,
    m=4,
    a = list(c(AMT=20),c(AMT=40),c(AMT=80),c(AMT=160)),
    xt = c(0.5,2,3,6,24,36,72,120)#,
    #model_switch = c(2,2,2,2,2,2,2,2)
  ),
  optimize=list(target="poped_R",
                model = list(
                  ff_file="PKPD.ff",
                  fError_file="feps",
                  fg_file="PKPD.fg"
                ),
                design_space=list(mina=1,
                                  maxa=500
                ),
                parameters=list(
                  bpop=NULL,
                  d=NULL,
                  sigma=NULL,# initial parameter guess not coming from previous step
                  manipulation=NULL,
                  ds_index = c(1,1,1,0,0,0,0)# manipulation of initial parameters
                ),
                settings.db=list(
                  notfixed_sigma=c(1,0),
                  notfixed_bpop=c(0,0,0,1,1,1,1),
                  notfixed_d=c(0,0,1,1),
                  ds_index = c(1,1,1,0,0,0,0)
                ),
                settings.opt=list(
                  opt_xt=F,
                  opt_a=T,
                  opt_x=F,
                    
                  bUseRandomSearch= 1,
                  bUseStochasticGradient = 0,
                  bUseBFGSMinimizer = 0,
                  bUseLineSearch = 1,
                  bUseExchangeAlgorithm=0,
                  compute_inv=F
                )
  ),
  simulate=list(target="NONMEM", model="NONMEM_files/ex1_sim.mod",
                data=list(dosing = list(list(AMT0=1,Time=0)),
                          manipulation = list(expression(AMT <- AMT0*AMT))
                )
                
  ),
  estimate=list(target="NONMEM", model="NONMEM_files/ex1_est_full.mod")
)


results_mbaod_true <- mbaod_simulate(cohorts=list(step_1,step_2), # only the sim/est step in this case.
                                      ncohorts=10, # number of steps or cohorts in one AOD, if greater than number of included cohorts, last step is repeated
                                      rep=1, #number of times to repeat the MBAOD simulation 
                                      name="PKPD_devtest", 
                                      description="Test of Parameterization",
                                      seednr=123)




































#set path to the directory where this file is located
setwd("C:/Users/erist836/Desktop/MBAOD2")

# remove things from the global environment
rm(list=ls())

library(PopED)
library(mvtnorm)
library(xpose4)

# load the MBAOD package, change to your relevant load call
devtools::load_all("C:/Users/erist836/Documents/GitHub/MBAOD/R")

# load the PopED model file
source("PopED_files/PKPD_example_1.R")
source("stop_crit_dose2.R")
source("mbaod_simulate.R")
############################################RUN THIS SECTION FOR MBAOD WITH SMALL Misspecification############################################
# 
# step_1=list(
#   design = list(
#     groupsize = 8,
#     m=2,
#     a = list(c(AMT=0.00001),c(AMT=140)),
#     xt = c(0.5,2,3,6,24,36,72,120)
#   ),
#   optimize=NULL,
#   simulate=list(target="NONMEM", model="NONMEM_files/ex1_sim.mod",
#                 data=list(dosing = list(list(AMT0=1,Time=0)),
#                           manipulation = list(expression(AMT <- AMT0*AMT))
#                 )
#                 
#   ),
#   estimate=list(target="NONMEM", model="NONMEM_files/ex1_est_full.mod")
# )
# 
# 


#First cohort with dose optimization
step_1 = list(
  design = list(
    groupsize = 4,
    m=2,
    a = list(c(AMT=20),c(AMT=160)),
    xt = c(0.5,2,3,6,24,36,72,120)#,
    #model_switch = c(2,2,2,2,2,2,2,2)
  ),
  optimize=list(target="poped_R",
                model = list(
                  ff_file="PKPD.ff",
                  fError_file="feps",
                  fg_file="PKPD.fg"
                ),
                design_space=list(mina=0.0001,
                                  maxa=500
                ),
                parameters=list(sigma = diag(c(0.015,0.0001)),
                                bpop = c(CL = 0.15, V=8,Ka=1,Base= 1.5,EMAX=150,EC50=10.5, hill=3),
                                d=c(CL = 0.07, V=0.02, EMAX=0.0625, EC50=0.0625)
                  
                ), 
                settings.db=list(
                  notfixed_sigma=c(1,0),
                  notfixed_bpop=c(0,0,0,1,1,1,1),
                  notfixed_d=c(0,0,1,1)#,
                  #ds_index = c(1,1,1,0,0,0,0)
                ),
                settings.opt=list(
                  opt_xt=F,
                  opt_a=T,
                  opt_x=F,
                  iFIMCalculationType=0,
                  ofv_calc_type=4,
                  bUseRandomSearch= 1,
                  bUseStochasticGradient = 0,
                  bUseBFGSMinimizer = 0,
                  bUseLineSearch = 1,
                  bUseExchangeAlgorithm=0,
                  compute_inv=F,
                  d_switch = 1
                )
  ),  
  simulate=list(target="NONMEM", model="NONMEM_files/ex1_sim.mod",
                data=list(dosing = list(list(AMT0=1,Time=0)),
                          manipulation = list(expression(AMT <- AMT0*AMT))
                )
                
  ),  
  estimate=list(target="NONMEM", model="NONMEM_files/ex1_est_full_misspec_odsc.mod")
)  
   
step_2 <- step_1
step_2$design = list(
  groupsize = 2,
  m=1,
  a = list(c(AMT=300)),
  xt = c(0.5,2,3,6,24,36,72,120)#,
  #model_switch = c(2,2,2,2,2,2,2,2)
)
step_2$optimize$parameters <- NULL


results_mbaod_true <- mbaod_simulate(cohorts=list(step_1,step_2), # only the sim/est step in this case.
                                     ncohorts=50, # number of steps or cohorts in one AOD, if greater than number of included cohorts, last step is repeated
                                     rep=50, #number of times to repeat the MBAOD simulation 
                                     name="PKPD_D_2in1g_misspec_odsc1", ED_only_FIM =T, 
                                     description="Ds optimality of dose", stop_crit_fun = stop_dose,
                                     seednr=123)

effs <- rep(0,length(results_mbaod_true$iteration_1$stop_res))
effs2 <- effs
for(i in 1:length(results_mbaod_true$iteration_1$stop_res)){
effs[i]<- results_mbaod_true$iteration_1$stop_res[[paste("cohort_",i,sep="")]][[4]]
effs2[i]<- results_mbaod_true$iteration_1$stop_res[[paste("cohort_",i,sep="")]][[3]]

}
effs
effs2







#set path to the directory where this file is located
setwd("C:/Users/erist836/Desktop/MBAOD2")

# remove things from the global environment
rm(list=ls())

library(PopED)
library(mvtnorm)
library(xpose4)

# load the MBAOD package, change to your relevant load call
devtools::load_all("C:/Users/erist836/Documents/GitHub/MBAOD/R")

# load the PopED model file
source("PopED_files/PKPD_example_1.R")
source("stop_crit_dose2.R")
source("mbaod_simulate.R")

#First cohort with dose optimization
step_1 = list(
  design = list(
    groupsize = 4,
    m=2,
    a = list(c(AMT=20),c(AMT=160)),
    xt = c(0.5,2,3,6,24,36,72,120)#,
    #model_switch = c(2,2,2,2,2,2,2,2)
  ),
  optimize=list(target="poped_R",
                model = list(
                  ff_file="PKPD.ff",
                  fError_file="feps",
                  fg_file="PKPD.fg"
                ),
                design_space=list(mina=0.0001,
                                  maxa=500
                ),
                parameters=list(sigma = diag(c(0.015,0.0001)),
                                bpop = c(CL = 0.15, V=8,Ka=1,Base= 1.5,EMAX=150,EC50=10.5, hill=3),
                                d=c(CL = 0.07, V=0.02, EMAX=0.0625, EC50=0.0625)
                                
                ), 
                settings.db=list(
                  notfixed_sigma=c(1,0),
                  notfixed_bpop=c(0,0,0,1,1,1,1),
                  notfixed_d=c(0,0,1,1)#,
                  #ds_index = c(1,1,1,0,0,0,0)
                ),
                settings.opt=list(
                  opt_xt=F,
                  opt_a=T,
                  opt_x=F,
                  iFIMCalculationType=0,
                  ofv_calc_type=4,
                  bUseRandomSearch= 1,
                  bUseStochasticGradient = 0,
                  bUseBFGSMinimizer = 0,
                  bUseLineSearch = 1,
                  bUseExchangeAlgorithm=0,
                  compute_inv=F,
                  d_switch = 0,
                  ED_samp_size =10
                )
  ),  
  simulate=list(target="NONMEM", model="NONMEM_files/ex1_sim.mod",
                data=list(dosing = list(list(AMT0=1,Time=0)),
                          manipulation = list(expression(AMT <- AMT0*AMT))
                )
                
  ),  
  estimate=list(target="NONMEM", model="NONMEM_files/ex1_est_full.mod")
)  

step_2 <- step_1
step_2$design = list(
  groupsize = 2,
  m=1,
  a = list(c(AMT=300)),
  xt = c(0.5,2,3,6,24,36,72,120)#,
  #model_switch = c(2,2,2,2,2,2,2,2)
)
step_2$optimize$parameters <- NULL


results_mbaod_ED<- mbaod_simulate(cohorts=list(step_1,step_2), # only the sim/est step in this case.
                                     ncohorts=50, # number of steps or cohorts in one AOD, if greater than number of included cohorts, last step is repeated
                                     rep=50, #number of times to repeat the MBAOD simulation 
                                     name="PKPD_ED_misspec", 
                                     description="Ds optimality of dose", stop_crit_fun = stop_dose,
                                     seednr=1234)