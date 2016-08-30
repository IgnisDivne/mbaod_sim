rm(list=ls())

#set path to the directory where this file is located
if(Sys.getenv("LOGNAME")=="ahooker"){
  setwd("~/Documents/_PROJECTS/AOD/repos/mbaod_sim/non-linear")
} else {
  setwd("C:/Users/erist836/Documents/GitHub/mbaod_sim/non-linear")  
}

# remove things from the global environment
rm(list=ls())

if(Sys.getenv("LOGNAME")=="ahooker"){
  devtools::load_all("~/Documents/_PROJECTS/PopED/repos/PopED")
} else {
  library(PopED)
}
library(mvtnorm)

PK.1.comp.maturation.fg <- function(x,a,bpop,b,bocc){
  ## -- parameter definition function
  
  parameters=c( CL=exp(bpop[1])*exp(b[1]),
                V=exp(bpop[2])*exp(b[2]),
                TM50=exp(bpop[3]),
                HILL=exp(bpop[4]),
                age_group = a[1],
                PMA = age_group_2_mean_PMA(a[1]),
                WT = get_median_weight(age_group_2_mean_PMA(a[1])))
  return( parameters )
}


PK.1.comp.maturation.ff <- function(model_switch,xt,parameters,poped.db){
  with(as.list(parameters),{
    y=xt
    
    CL=CL*(WT/70)^0.75*((PMA^HILL)/(PMA^HILL + TM50^HILL))
    V=V*(WT/70)
    
    
    DOSE=1000*(WT/70)
    y = DOSE/V*exp(-CL/V*xt)
    
    return(list( y= y,poped.db=poped.db))
  })
}


# load the weight estimator file
source("get_weight.R")

propofol_t.db <- create.poped.database(ff_file="PK.1.comp.maturation.ff", 
                                       fError_file="feps.add.prop",
                                       fg_file="PK.1.comp.maturation.fg",
                                       groupsize = 6,
                                       m=6,
                                       sigma = diag(c(0.015,0.0001)),
                                       bpop = c(CL = 1, V=3, TM50=log(38.5), HILL=log(4.6)),
                                       d=c(CL=0.05,V=0.05),
                                       #notfixed_bpop = c(1,1,1,1),
                                       #notfixed_d = c(1,1),
                                       notfixed_sigma = c(1,0),
                                       xt=c(0.1, 2, 6, 12, 24),
                                       minxt = 0,
                                       maxxt = 24,
                                       #x=t(rbind(age_group=c(1,2,3,4,5,6))),
                                       #nx=1,
                                       #discrete_x = x.space,
                                       a=rbind(1,2,3,4,5,6),
                                       discrete_a = list(c(1,2,3,4,5,6))
                                       )

##  create plot of model without variability 
plot_model_prediction(propofol_t.db)


##  create plot of model with variability 
#plot_model_prediction(propofol_t.db,IPRED=T,DV=T,separate.groups=T)

## evaluate initial design
# FIM <- evaluate.fim(propofol_t.db) 
# FIM
# det(FIM)
# get_rse(FIM,propofol_t.db)
# 
# ls_step_size <- 1#ifelse(fast,3,50)
# propofol_T_OD_res <- poped_optimize(propofol_t.db,opt_xt=0,opt_a=0, opt_x=1,
#                                     bUseRandomSearch= 0,
#                                     bUseStochasticGradient = 0,
#                                     bUseBFGSMinimizer = 0,
#                                     bUseLineSearch = 1,
#                                     ls_step_size=ls_step_size)   #gs= 12 12 12   x=6 3 2
# 
# get_rse(propofol_T_OD_res$fmf,propofol_T_OD_res$poped.db)
# plot_model_prediction(propofol_T_OD_res$poped.db)
# 

# slow optimization
propofol_T_OD_res <- poped_optim(propofol_t.db,opt_xt=0,opt_a=T,parallel = T)


# faster optimization
propofol_T_OD_res_fast <- poped_optim(propofol_t.db,opt_xt=0,opt_a=T,parallel = T,
                                      method=c("LS"), loop_methods=T)


get_rse(propofol_T_OD_res_fast$FIM,propofol_T_OD_res_fast$poped.db)
plot_model_prediction(propofol_T_OD_res_fast$poped.db)
