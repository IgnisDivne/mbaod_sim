



# remove things from the global environment
rm(list=ls())

if(Sys.getenv("LOGNAME")=="ahooker"){
  devtools::load_all("~/Documents/_PROJECTS/PopED/repos/PopED")
} else {
  devtools::load_all("~/PopED")
}
library(mvtnorm)
library(xpose4)

# load the MBAOD package change to your relevant load call
if(Sys.getenv("LOGNAME")=="ahooker"){
  devtools::load_all("~/Documents/_PROJECTS/AOD/repos/MBAOD")
} else {
  devtools::load_all("~/MBAOD")
}


# load the PopED model file
source("PopED_files/pkpd_example_1.R")


source("stop_crit_dose2.R")
source("mbaod_simulate.R")


doses <- seq(0,500,0.5)
doses[1] <- 0.0001

a.space <- cell(2,1)
a.space[,] <- list(doses)
#First cohort with dose optimization
step_1 = list(
  design = list(
    groupsize = 4,
    m=2,
    a = list(c(AMT=20),c(AMT=160)),
    xt = c(0.5,3,60)#,
    #model_switch = c(2,2,2,2,2,2,2,2)
  ),
  optimize=list(target="poped_R",
                model = list(
                  ff_file="PKPD.ff",
                  fError_file="feps",
                  fg_file="PKPD.fg"
                ),
                design_space=list(a_space = a.space
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
                  parallel=F,
                  method=c("LS","ARS"), 
                  loop_methods=T,
                  eff_crit = 1.0001,
                  compute_inv=F,
                  iFIMCalculationType=0,
                  ofv_calc_type=4,
                  d_switch = 1
                )
  ),  
  simulate=list(target="NONMEM", model="NONMEM_files/ex1_sim.mod",
                data=list(dosing = list(list(AMT0=1,Time=0)),
                          manipulation = list(expression(AMT <- AMT0*AMT))
                )
                
  ),  
  estimate=list(target="NONMEM", model="NONMEM_files/ex1_est_full_misspec.mod")
)  

step_2 <- step_1
step_2$design = list(
  groupsize = 2,
  m=1,
  a = list(c(AMT=300)),
  xt = c(0.5,3,60)#,
  #model_switch = c(2,2,2,2,2,2,2,2)
)

a.space2 <- cell(1,1)
a.space2[,] <- list(doses)

step_2$optimize$parameters <- NULL
step_2$optimize$design_space$a_space <- a.space2


results_mbaod_true <- mbaod_simulate(cohorts=list(step_1,step_2), # only the sim/est step in this case.
                                     ncohorts=50, # number of steps or cohorts in one AOD, if greater than number of included cohorts, last step is repeated
                                     rep=100, #number of times to repeat the MBAOD simulation 
                                     name="PKPD_D_2in1g_misspec_odsc1", ED_only_FIM =F, 
                                     description="Ds optimality of dose", stop_crit_fun = stop_dose,
                                     seednr=123)

# effs <- rep(0,length(results_mbaod_true$iteration_1$stop_res))
# effs2 <- effs
# for(i in 1:length(results_mbaod_true$iteration_1$stop_res)){
# effs[i]<- results_mbaod_true$iteration_1$stop_res[[paste("cohort_",i,sep="")]][[4]]
# effs2[i]<- results_mbaod_true$iteration_1$stop_res[[paste("cohort_",i,sep="")]][[3]]
# 
# }
# effs
# effs2





