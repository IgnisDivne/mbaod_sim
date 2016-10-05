

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



pkpd.db <- create.poped.database(ff_file="PKPD.ff",    
                                 fError_file="feps",
                                 fg_file="PKPD.fg",
                                 groupsize = 10,
                                 m=1,
                                 xt=c(0.5,3,60),
                                 sigma = diag(c(0.015,0.0001)),
                                 bpop = c(CL = 0.15, V=8,Ka=1.5,Base= 1,EMAX=150,EC50=10.5, hill=3),
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



propofol_T_OD_res <- poped_optim(pkpd.db,opt_a=T,opt_xt = F,compute_inv=F, eff_crit = 1.0001, methods=c("ARS","LS"), loop_methods = T,ofv_calc_type = 4,d)   

xt=c(0.5,3,60)
plot_model_prediction(poped.db=pkpd.db) 
