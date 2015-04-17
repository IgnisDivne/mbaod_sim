
##Simulated bridging study of IV admin. drug X with true 1 comp kinetics from
# 25 y.o. adults with weight 70 to infants 3 monts -2 years old and
## simulated values (the truth) in this case will be c(TM50=80 (WEEKS),hill=2.4)

##Added and testing function which allows for adding children from the agegroup one below the lowest group which has passed the
##stoppign criteria...see stop_critX.R and stop_crit_POPED.R. This seems to reduce the total number of children.

#set path to the directory where this file is located
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


#Adults
step_1=list(
  design = list(
    groupsize = 100,
    x = c(age_group=c(7)),#6=adults ##t(rbind(PMA=c(1346)),WT=c(70))),  #AGE in PMA (weeks) assuming normal 39 week pregnancy
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

##Initial Design space for the age groups. subadults 12-18 y.o. and adults.
##If only adults and subadults are allowed, the information about TM50 is too sparse.
x.space <- cell(1,1)
x.space[1,1] <- list(c(6,7))
###

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
                  bpop=c(TM50=75), # initial parameter guess not coming from previous step
                  manipulation=NULL # manipulation of initial parameters
                ),
                settings.db=NULL,
                settings.opt=list(
                  opt_xt=T,
                  opt_a=F,
                  opt_x=T,
                  bUseRandomSearch= 1,
                  bUseStochasticGradient = 0,
                  bUseBFGSMinimizer = 0,
                  bUseLineSearch = 1,
                  compute_inv=T
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
  estimate=list(target="NONMEM", model="C:/Users/Eric Stromberg/Desktop/MBAOD_project/NONMEM_files/est_full_log_small.mod")
)


step_3 <- step_2
step_3$optimize$parameters <- NULL
step_3$design$groupsize <- 2


results_mbaod <- mbaod_simulate(cohorts=list(step_1,step_2,step_3), # anything after step_3 is the same as step_3
                              ncohorts=5, # number of steps or cohorts in one AOD
                              rep=1, #number of times to repeat the MBAOD simulation 
                              name="bridging_maturation_model_xoptim_restricted_group", 
                              description="25 steps, 1st step one group, steps 2-10 have 1 added group per step",
                              seednr=321, stop_crit_fun =stop_critX)


############################################################################

get_age_from_median_weight(get_median_weight(age_group_2_mean_PMA(c(1,5,10,18))))

###########################Regular OD#####################################

#set path to the directory where this file is located
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
source("stop_crit_PopED.R")


#Adults
step_1=list(
  design = list(
    groupsize = 100,
    m=1,
    x = c(age_group=c(7)),#6=adults ##t(rbind(PMA=c(1346)),WT=c(70))),  #AGE in PMA (weeks) assuming normal 39 week pregnancy
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
                  bpop=c(TM50=75),# initial parameter guess not coming from previous step
                  manipulation=NULL # manipulation of initial parameters
                ),
                settings.db=NULL,
                settings.opt=list(
                  opt_xt=F,
                  opt_a=F,
                  opt_x=T,
                  bUseRandomSearch= 1,
                  bUseStochasticGradient = 0,
                  bUseBFGSMinimizer = 0,
                  bUseLineSearch = 1,
                  compute_inv=T
                )
  ),
  simulate=NULL,
  estimate=NULL
)

step_3 <- step_2
step_3$optimize$parameters <- list(
                                    bpop=c(CL = 1.122,V=5.098,TM50=75),
                                    d=c(0.0950218, 0.0612376),
                                    sigma=matrix(c(0.0144145,0,0,0.0001),nrow=2,byrow=T),
                                    manipulation=NULL)
step_3$design$groupsize <- 1



results_all_OD_2 <- mbaod_simulate(cohorts=list(step_1,step_2,step_3), # anything after step_3 is the same as step_3
                              ncohorts=3, # number of steps or cohorts in one AOD
                              rep=1, #number of times to repeat the MBAOD simulation 
                              name="bridging_maturation_model_Xlog_OD", 
                              description="10 steps, 1st step one group, steps 2-10 have 2 groups per step",
                              seednr=321, stop_crit_fun = stop_crit_PopED)


# ## now make some plots
# ## here are some examples
# 
# library(grid)
library(ggplot2)
#library(reshape2)
full_design_list <- lapply(mapply("[[",design_list,design_name,SIMPLIFY=FALSE),function(x){do.call(create_design,x)})

#################################
######## optimized designs
#################################

#all_designs
design_list <- results_all[grep("^iteration",names(results_all))]

all_designs <- combine_designs(design_list,design_name = "final_design")

model = list(
  ff_file="PK.1.comp.maturation.ff",
  fError_file="feps.add.prop",
  fg_file="PK.1.comp.maturation.fg"
)

parameters_true=list(
  bpop=c(CL=1,V=20,EMAX=2,EC50=25,HILL=5),
  d=c(0.05,0.05),
  sigma=c(0.015,0.0015)
)

poped.db <- do.call(create.poped.database,c(all_designs,model,parameters_true))

plot1 <- plot_model_prediction(poped.db,y_lab="Concentration")
plot1 + theme(legend.position="none")



#################################
######## PARAMETER ESTIMATES 
#################################

true_values <- c(thetas=c(1,20,2,25,5),
                 omegas=sqrt(c(0.05,0.05)),
                 sigmas=sqrt(c(0.015,0.0015)))

plot_parameter_estimates(results_all,true_values)

#################################
######## VPC of IPRED from estimated models and true model
#################################

design_1 = list(
  groupsize = 200,
  m=1,
  a   = 35,
  xt = c(0.5,1,2,3,6,12,24)
)

design_2 = list(
  groupsize = 200,
  m=4,
  a   = rbind(10, 35, 55, 70),
  xt = c(0.5,1,2,3,6,12,24)
)

model = list(
  ff_file="PK.1.comp.maturation.ff",
  fError_file="feps.add.prop",
  fg_file="PK.1.comp.maturation.fg"
)

parameters_true=list(
  bpop=c(CL=1,V=20,EMAX=2,EC50=25,HILL=5),
  d=c(0.05,0.05),
  sigma=c(0.015,0.0015)
)

mbaod_vpc(design_1, 
          model, 
          parameters_true, 
          results_all)

mbaod_vpc(design_2, 
          model, 
          parameters_true, 
          results_all, 
          separate.groups=T)


# #################################
# ######## Clearance plots (specific for this problem) -- visualization of WT choices
# #################################

CL_mod <- function(params=list(BASE=1,EMAX=2,E50=25,HILL=5),IDV){
  with(params,{
    vals <- BASE+ (EMAX*IDV^HILL)/(E50^HILL + IDV^HILL)
    return(vals)
  })
}

df <- data.frame(WT=0:70)
df$CL=CL_mod(IDV=df$WT)

#all_designs
design_list <- results_all[grep("^iteration",names(results_all))]
all_designs <- combine_designs(design_list,design_name = "final_design")

df.2 <- data.frame(all_designs$a)
df.2$CL=CL_mod(IDV=df.2$WT)
nrep <- length(grep("^iteration",names(results_all)))
ncohort <- size(df.2,1)/nrep
df.2$Cohort=as.factor(rep(1:ncohort,nrep))  

p <- ggplot(data=df, aes(x=WT,y=CL))
p <- p+geom_line()
p+geom_point(data=df.2,aes(color=Cohort),size=4)

