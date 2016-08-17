# 
# 
# 
# # remove things from the global environment
# rm(list=ls())
# setwd("Dir of the files")
# library(PopED)
# library(mvtnorm)
# library(xpose4)
# 
# # load the PopED model file
# source("PopED_files/poped.mod.PK.1.comp.maturation_real.R")
# # load the weight estimator file
# source("get_weight.R")
# 
# 
# 
# 
# 
# 
# #Adults
# step_1=list(
#   design = list(
#     groupsize = 100,
#     x = c(age_group=c(7)),
#     xt = c(0.1, 2, 6, 12, 24)
#   ),
#   optimize=NULL,
#   simulate=list(target="NONMEM", model="NONMEM_files/sim_propofol.mod",
#                 data=list(dosing = list(list(AMT=1000,Time=0)),
#                           manipulation = list(expression(PMA <- age_group_2_PMA(ID,age_group)),
#                                               expression(WT <- get_weight(ID,PMA,probs = c(0.1,0.9))),
#                                               expression(AMT <- AMT*(WT/70))
#                           )
#                 )
#                 
#   ),
#   estimate=list(target="NONMEM", model="NONMEM_files/est_red_propofol.mod")
# )
# 
# ##Initial Design space for the age groups. subadults 12-18 y.o. and adults.
# ##If only adults and subadults are allowed, the information about TM50 is too sparse.
# x.space <- cell(5,1)
# x.space[,] <- list(c(2,3,4,5,6,7))
# ###
# 
# #First cohort with children
# step_2 = list(
#   design = list(
#     groupsize = 2,
#     m=5,
#     x = t(rbind(age_group=c(6,5,4,3,2))),
#     xt = c(0.1, 2, 6, 12, 24)
#   ),
#   optimize=list(target="poped_R",
#                 model = list(
#                   ff_file="PK.1.comp.maturation.ff",
#                   fError_file="feps.add.prop",
#                   fg_file="PK.1.comp.maturation.fg"
#                 ),
#                 design_space=list(minxt=0,
#                                   maxxt=24,
#                                   x_space = x.space
#                 ),
#                 parameters=list(
#                   bpop=c(TM50=3.862, HILL=1.224), # initial parameter guess not coming from previous step
#                   manipulation=NULL # manipulation of initial parameters
#                 ),
#                 settings.db=list(
#                   notfixed_sigma=c(1,0)
#                 ),
#                 settings.opt=list(
#                   opt_xt=F,
#                   opt_a=F,
#                   opt_x=T,
#                   bUseRandomSearch= 0,
#                   bUseStochasticGradient = 0,
#                   bUseBFGSMinimizer = 0,
#                   bUseLineSearch = 1,
#                   bUseExchangeAlgorithm=0,
#                   EACriteria = 1,
#                   compute_inv=T
#                 )
#   ),
#   simulate=list(target="NONMEM", model="NONMEM_files/sim_propofol.mod",
#                 data=list(dosing = list(list(AMT=1000,Time=0)),
#                           manipulation = list(expression(PMA <- age_group_2_PMA(ID,age_group)),
#                                               expression(WT <- get_weight(ID,PMA,probs = c(0.1,0.9))),
#                                               expression(AMT <- AMT*(WT/70))
#                           )
#                 )
#                 
#   ),
#   estimate=list(target="NONMEM", model="NONMEM_files/est_full_gfr_OD.mod")
# )
# 
# x.space <- cell(2,1)
# x.space[,] <- list(c(2,3,4,5,6,7))
# step_3 <- step_2
# step_3$design$groupsize <- 1
# step_3$optimize$parameters <- NULL
# step_3$design$x <- t(rbind(age_group=c(6,6)))
# step_3$design$m <- 2
# step_3$optimize$design_space$x_space<- x.space
# 
# 
# results_OD_GFR <- mbaod_simulate(cohorts=list(step_1,step_2,step_3), # anything after step_3 is the same as step_3
#                                  ncohorts=100, # number of steps or cohorts in one AOD
#                                  rep=34, #number of times to repeat the MBAOD simulation 
#                                  name="gfr_run2_OD_crit", lower_limit=0.6,higher_limit=1.4,
#                                  description="25 steps, 1st step one group, steps 2-10 have 1 added group per step",
#                                  seednr=1234, stop_crit_fun =stop_critX_2,run_commands="-retries=5", #123 1-66, 
#                                  ci_corr=0,option=3)
# 






#set path to the directory where this file is located
setwd("C:/Users/erist836/Desktop/MBAOD_project")
# remove things from the global environment
rm(list=ls())

library(PopED)
library(mvtnorm)
library(xpose4)

# load the PopED model file
source("PopED_files/poped.mod.PK.1.comp.maturation_real2.R")
# load the weight estimator file
source("get_weight.R")

x.space <- cell(6,1)
x.space[,] <- list(c(1,2,3,4,5,6))

propofol_t.db <- create.poped.database(ff_file="PK.1.comp.maturation.ff", 
                                       fError_file="feps.add.prop",
                                       fg_file="PK.1.comp.maturation.fg",
                                       groupsize = 6,
                                       m=6,
                                       sigma = diag(c(0.015,0.0001)),
                                       bpop = c(CL = 1, V=3, TM50=38.5, HILL=4.6),
                                       d=c(CL=0.05,V=0.05),
                                       notfixed_bpop = c(1,1,1,1),
                                       notfixed_d = c(1,1),
                                       notfixed_sigma = c(1,0),
                                       xt=c(0.1, 2, 6, 12, 24),
                                       minxt = c(0,0,0,0,0),
                                       maxxt = c(24,24,24,24,24),
                                       x=t(rbind(age_group=c(1,2,3,4,5,6))),
                                       nx=1,
                                       discrete_x = x.space)  


ls_step_size <- 1#ifelse(fast,3,50)
propofol_T_OD_res <- poped_optimize(propofol_t.db,opt_xt=0,opt_a=0, opt_x=1,
                                    bUseRandomSearch= 0,
                                    bUseStochasticGradient = 0,
                                    bUseBFGSMinimizer = 0,
                                    bUseLineSearch = 1,
                                    ls_step_size=ls_step_size)   #gs= 12 12 12   x=6 3 2



x.space <- cell(6,1)
x.space[,] <- list(c(1,2,3,4,5,6))

propofol_m.db <- create.poped.database(ff_file="PK.1.comp.maturation.ff", 
                                       fError_file="feps.add.prop",
                                       fg_file="PK.1.comp.maturation.fg",
                                       groupsize = 6,
                                       m=6,
                                       sigma = diag(c(0.015,0.0001)),
                                       bpop = c(CL = 1, V=3, TM50=47.6, HILL=3.4),
                                       d=c(CL=0.05,V=0.05),
                                       notfixed_bpop = c(1,1,1,1),
                                       notfixed_d = c(1,1),
                                       notfixed_sigma = c(1,0),
                                       xt=c(0.1, 2, 6, 12, 24),
                                       minxt = c(0,0,0,0,0),
                                       maxxt = c(24,24,24,24,24),
                                       x=t(rbind(age_group=c(1,2,3,4,5,6))),
                                       nx=1,
                                       discrete_x = x.space) 


propofol_M_OD_res <- poped_optimize(propofol_m.db,opt_xt=0,opt_a=0, opt_x=1,
                                    bUseRandomSearch= 0,
                                    bUseStochasticGradient = 0,
                                    bUseBFGSMinimizer = 0,
                                    bUseLineSearch = 1,
                                    ls_step_size=ls_step_size)



























#set path to the directory where this file is located
setwd("C:/Users/erist836/Desktop/MBAOD_project")
# remove things from the global environment
rm(list=ls())

library(PopED)
library(mvtnorm)
library(xpose4)

# load the PopED model file
# load the weight estimator file

PKPD.fg <- function(x,a,bpop,b,bocc){
  ## -- parameter definition function
  
  parameters=c( CL   = bpop[1]*exp(b[1]),
                V    = bpop[2]*exp(b[2]),
                Ka   = bpop[3],#*exp(b[3]),
                Base = bpop[4],#*exp(b[3]),
                EMAX = bpop[5]*exp(b[3]),
                EC50 = bpop[6]*exp(b[4]),
                hill = bpop[7],#*exp(b[5]),
                Dose = a[1]
  )
  return( parameters )
}


PKPD.ff <- function(model_switch,xt,parameters,poped.db){
  with(as.list(parameters),{
    
    y=xt
    Ke=CL/V
    
    C = ((Dose*Ka*Ke)/(CL*(Ka-Ke)))*(exp(-Ke*xt)-exp(-Ka*xt))
    y = Base+((C^hill*EMAX)/(EC50^hill + C^hill))
    
    
    
    return(list( y= y,poped.db=poped.db))
  })
}



feps <- function(model_switch,xt,parameters,epsi,poped.db){
  ## -- Residual Error function
  ## -- Proportional PK + additive PD
  returnArgs <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db)) 
  y <- returnArgs[[1]]
  poped.db <- returnArgs[[2]]
  
  
  y <- y*(1+epsi[,1])+epsi[,2]
  
  
  
  
  return(list( y= y,poped.db =poped.db )) 
}



poped.db <- create.poped.database(ff_file="PKPD.ff",    
                                  fError_file="feps",
                                  fg_file="PKPD.fg",
                                  groupsize = 10,
                                  m=1,
                                  xt=c(0.5,2,3,6,24,36,72,120),
                                  sigma = diag(c(0.015,0.0001)),
                                  bpop = c(CL = 0.15, V=8,Ka=1,Base= 1,EMAX=100,EC50=7, hill=2),
                                  d=c(CL = 0.07, V=0.02, EMAX=0.0625, EC50=0.0625),
                                  notfixed_bpop = c(0,0,0,1,1,1,1),
                                  notfixed_d = c(0,0,1,1),
                                  notfixed_sigma = c(1,0),
                                  minxt = 0,
                                  maxxt = 121,
                                  a= 160,
                                  mina= 0,
                                  ds_index = c(1,1,1,0,0,0,0), 
                                  maxa=500,
                                  bUseRandomSearch= 1,
                                  bUseStochasticGradient = 0,
                                  bUseExchangeAlgorithm = 0, 
                                  EAStepSize=1,
                                  bUseBFGSMinimizer = 0,
                                  ofv_calc_type = 6,
                                  bUseLineSearch = 1
                                  
)  



propofol_T_OD_res <- poped_optim(pkpd.db,opt_a=T,opt_xt = F,compute_inv=F, eff_crit = 1.0001, methods=c("ARS", ) loop_methods = T,ofv_calc_type = 6)   


pkpd_Res <- poped_optimize(pkpd.db,opt_a = T,compute_inf=F
)  

plot_model_prediction(poped.db) 
