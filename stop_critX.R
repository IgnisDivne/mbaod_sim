##Stopping function for MBAOD bridging study with 2 scaling parameters (Fixed WT hill on CL)
stop_critX <- function(cohort_num,cohort_res,option=2, nsim=100000, age_space = c(1,2,3,4,5,6),allow_next=T){
  if(cohort_num>1){
  print(paste("--::::Checking Stopping Criteria for Cohort", cohort_num,"::::--"))
  results_est <- cohort_res$est_result
  
  omega <- results_est$omega
  omega <- matrix_from_triangle(omega)
  seomegas <- results_est$seomegas
  seomegas <- matrix_from_triangle(seomegas)
  
  
  thetacl   <- results_est$thetas[1] #logcl
  thetav    <- results_est$thetas[2]
  thetaTM50 <- results_est$thetas[3]
  #thetaHILL <- results_est$thetas[4]
  omegacl<-omega[1,1]
  omegav <-omega[2,2]
  
  SEthetacl <- results_est$sethetas[1]
  SEthetav  <- results_est$sethetas[2]
  SETM50    <- results_est$sethetas[3]
  
  corr <- results_est$cov_mat[1,3]
  lcl_params <- c("THETA1","THETA3")#,"THETA4")
  
  covmat <- cohort_res$est_result$cov_mat
  cov_lcl <- covmat[rownames(covmat) %in% lcl_params,colnames(covmat) %in% lcl_params]
  cov_lcl <- as.matrix(cov_lcl)
  
  SEomegacl <- seomegas[1,1] 
  SEomegav  <- seomegas[2,2]
  
#Option 1 - Simulation based
#Simulatea large number of individuals with uncertainty based
#on the latest design and parameter estimates and parameter uncertainty
#Get parameters used for each individual and calculate geometrical mean and 95% CI
#for each parameter
#Compares the 95% CI to the range of 60% of 140% of geometrical mean for each parameter
#If the CI is within the mean range for all parameters and power is 80% or more,
#the function stops the MBAOD after this cohort. 
 
if(option==1){  
  
  CI_tvcl <-rnorm(nsim,thetacl,SEthetacl)
  CI_tvv <- rnorm(nsim,thetav,SEthetav)  
  #CI of scaling parameters
  CI_tvTM50 <- rnorm(nsim,thetaTM50,SETM50)
  CI_tvHILL <- rnorm(nsim,thetaHILL,SEHILL)
  #97.5th percentile value of omega = w^2 according to standard errors
  omcl <- rnorm(nsim,omegacl,SEomegacl)
  omv <- rnorm(nsim,omegav,SEomegav)
  
  #97.5th percentile etas from eta ~ N(0,w^2) 
  etacl <- rnorm(nsim,0,sqrt(omcl))
  etav <-rnorm(nsim,0,sqrt(omv))
  #CI of the individual parameters
  CI_cl <- CI_tvcl*exp(etacl)
  CI_V <- CI_tvv*exp(etav)
  
  #CI of the individually scaled parameters
  iscl <- CI_cl*(WT/70)^0.75*(AGE^CI_tvHILL/(AGE^CI_tvHILL+CI_tvTM50^CI_tvHILL))
  isV <-  CI_V*(WT/70)
  
  
  CI_sCL <- quantile(iscl, probs=c(0.025,0.975))
  interval_iscl <-  c(0.6*gm_mean(iscl),1.4*gm_mean(iscl))
  
  CI_sV <- quantile(isV, probs=c(0.025,0.975))
  interval_isv <-  c(0.6*gm_mean(isV),1.4*gm_mean(isV))
  
  print("CI of individual scaled Clearance")
  print(CI_sCL)
  print("60%-140% interval of geometric mean")
  print(interval_iscl)
  print("CI of individual scaled V")
  print(CI_sV)
  print("60%-140% interval of geometric mean")
  print(interval_isv)
  
  if(CI_sCL[1] > interval_iscl[1] & CI_sCL[2] < interval_iscl[2] 
     & CI_sV[1] > interval_isv[1] & CI_sV[2] < interval_isv[2]){
    print("Parameter CI within 60%-140% of geometric mean?:")
    return(TRUE)
  }else{
    print("Parameter CI within 60%-140% of geometric mean?:")
   return(FALSE) 
  }
  
}
  
#Option 2
#From the estimates of thetas, etas and eps and their respective standard errors:
#Calculate the 2.5th 92.5th percentiles of what is possible as described by Wang et al.
#Compare the "CI" to 60% and 140% of geometrical mean 
#(Theta?)
if(option==2){
  
  ###This block calculates the total number of individuals in each age group
  
  gs      <- cohort_res$opt_result$opt_output$poped.db$downsized.design$groupsize
  n_adults <- gs[1]
  nsub <- sum(gs) - n_adults
  WT <- c(6.906,8.4,9.29,11.11,13.72,15.96,18.14,21.15,23.96,26.74,31.59,36.03,40.53,47.06,51.908,58.02,62.78,66.96,69.11,71.48,73.49)
  AGE <- get_age_from_median_weight(WT)
  AGE_WT <- data.frame(AGE,WT)
  age_groups <- data.frame(group = 1:length(unique(AGE$AGE)),AGE =unique(AGE$AGE))

#setup the vector keeping track of groups which has passed  
group_pass  <- rep(0,length(age_space))

#Storing all calculated CIs in vectors for plotting later
CL_CIs <- matrix(nrow=length(AGE_WT$AGE),ncol=4)
colnames(CL_CIs) <- c("Age Group","PMA","2.5th","97.5th")
V_CIs <- CL_CIs
k=1 #to keep track of vector posistion

#For all pediatric groups, group 1 is adults, hence i in 2:length(AGE_WT)
  for(i in 1:length(age_space)){
    cat("\n")
    print(paste("___ Pediatric Sub Population",i,": Children of age",age_groups$AGE[age_space[i]],"___"))
    ####For all children ages in that group
    sub_group <- subset(AGE_WT,AGE_WT$AGE==age_groups$AGE[age_space[i]])
    sub_counter <- 0
    for(j in 1:length(sub_group$AGE)){
    
      ##scaled LogCL SE calculation 
      thetap1=thetacl #log
      theta2=thetaTM50
      theta3=0.75#thetaHILL
      
      
      #input variance-covariance matrix for the 3 parameters above

        #define the weight and age combination for SE estimation
        wt=sub_group$WT[j]
        age=sub_group$PMA[j]
  
        #define covariate model in R

        f <- function(x,y,z) x + log(wt/70)*z + log(age/(age+y))
        df<- deriv(body(f), c("x","y"))
        x=thetap1
        y=theta2
        z=theta3
      
        out=eval(df)
        dfout=attr(out,"gradient")
        varlcl=dfout%*%cov_lcl%*%t(dfout)
        SElcl=sqrt(varlcl)
  
        #df <- nsub -num_params. is num_params 8? for this model 4 thetas, 2 omegas 2 sigmas
        CI_CL=c(exp(-qt(0.975,nsub-6)*SElcl),exp(qt(0.975,nsub-6)*SElcl))
        CI_V=c(exp(-qt(0.975,nsub-6)*SEthetav),exp(qt(0.975,nsub-6)*SEthetav))
        
    


      
      
      
      if(CI_CL[1] > 0.6 & CI_CL[2] < 1.4 
         & CI_V[1] > 0.6 & CI_V[2] < 1.4){
        sub_counter <- sub_counter +1
      }else{
        print(" All Parameter CIs NOT within 60%-140% (0.6-1.4) of geometric mean: ")
        print(paste("Scaled logCL for age",sub_group$AGE[j],":"))
        print(CI_CL)
        print(paste("Scaled logV for age",sub_group$AGE[j],":"))
        print(CI_V)
       
  
      }
      
      CL_CIs[k,] <- c(i,age,CI_CL)
      V_CIs[k,] <- c(i,age,CI_V)
      k=k+1
      
  } #end for all ages in this age group
  
  if(length(sub_group$AGE) == sub_counter){ #If parameter precision inall pediatric subgroups have been reached
   
    print(paste("Stopping criteria reached for sub population",i)) 
    group_pass[i] <- 1
  }else{
    print(paste("Stopping criteria NOT reached for population",i)) 
  }
    

  }###end for all age groups##############
  if(sum(group_pass) == length(age_space)){ #If parameter precision inall pediatric subgroups have been reached
    cat("\n")
    print("Stopping criteria reached for all groups, stopping MBAOD") 
    stop_MBAOD <- TRUE
  }else{
    cat("\n")
    print("Stopping criteria NOT reached for all groups, starting next cohort")
    stop_MBAOD <- FALSE
    
  }



  xspace <- group_pass*age_space
  #Unlocks the oldest group for which the stopping critera has not been fufulled
  #E.g. for the MBAOD to allow adding children of group 4, the older children in group 5 must first have a good estimate
  #Added this after the design got "stuck"
  if(allow_next == T) xspace[min(xspace[xspace>1])-1] <- min(xspace[xspace>1])-1

  new_xspace <- unique(c(xspace[xspace != 0],6,7))
  return(list(stop_MBAOD,new_xspace,CL_CIs,V_CIs))

}



#Option 3 Calculating CI of all parameters based on parameter estimates and  SE of estimates


if(option==3){  
  gs      <- cohort_res$opt_result$opt_output$poped.db$downsized.design$groupsize
  n_adults <- gs[1]
  nsub <- sum(gs) - n_adults
  WT <- c(6.906,8.4,9.29,11.11,13.72,15.96,18.14,21.15,23.96,26.74,31.59,36.03,40.53,47.06,51.908,58.02,62.78,66.96,69.11,71.48,73.49)
  AGE <- get_age_from_median_weight(WT)
  AGE_WT <- data.frame(AGE,WT)
  age_groups <- data.frame(group = 1:length(unique(AGE$AGE)),AGE =unique(AGE$AGE))
  
  
  stop_MBAOD <- 0
  
  #Simulate nsim typical values given the standard error SE
  cl <-thetacl + rt(nsim,nsim-6)* SEthetacl
  V <- thetav + rt(nsim,nsim-6)*SEthetav  
 TM50 <- thetaTM50 + rt(nsim,nsim-6)*SETM50
  
  #sim_tvHILL <- rnorm(nsim,thetaHILL,SEHILL)
  #nsim values of omega = w^2 according to standard errors
 # omcl <- rnorm(nsim,omegacl,SEomegacl)
  #omv <- rnorm(nsim,omegav,SEomegav)
  
  #nsim etas from  N(mean=0,std = w) 
  #etacl <- rnorm(nsim,0,sqrt(omcl))
 # etav <-rnorm(nsim,0,sqrt(omv))
  #nsim individual values of the parameters
  sim_cl <- sim_tvcl #+etacl
  sim_V <- sim_tvv #+etav
  #For all pediatric subgroups (group 1 = adults)

  #For all pediatric groups, group 1 is adults, hence i in 2:length(AGE_WT)
  for(i in 1:length(age_groups$group)){
    cat("\n")
    print(paste("___ Pediatric Sub Population",i,": Children of age",age_groups$AGE[i],"___"))
    ####For all children ages in that group
    sub_group <- subset(AGE_WT,AGE_WT$AGE==age_groups$AGE[i])
    sub_counter <- 0
    for(j in 1:length(sub_group$AGE)){
      
      ##scaled LogCL SE calculation 
      wt=sub_group$WT[j]
      age=sub_group$PMA[j]
      #nsim individually scaled parameters 
      logiscl <- cl + 0.75*log(wt/70)+log(age/(age+TM50))
      logisV <-  V+log(wt/70)
      
      #standard error of the scaled logparameters
      SE_lcl <- sd(logiscl)
      SE_lv <- sd(logisV)
      
  
      #df <- nsub -num_params. is num_params 8? for this model 4 thetas, 2 omegas 2 sigmas
      CI_CL2 =c(exp(-qt(0.975,nsub-6)*SE_lcl),exp(qt(0.975,nsub-6)*SE_lcl))
      CI_V2=c(exp(-qt(0.975,nsub-6)*SE_lv),exp(qt(0.975,nsub-6)*SE_lv))
      
      if(CI_CL[1] > 0.6 & CI_CL[2] < 1.4 
         & CI_V[1] > 0.6 & CI_V[2] < 1.4){
        stop_MBAOD <- stop_MBAOD +1
        sub_counter <- sub_counter +1
      }else{
        print(" All Parameter CIs NOT within 60%-140% (0.6-1.4) of geometric mean: ")
        print(paste("Scaled logCL for age",sub_group$AGE[j],":"))
        print(CI_CL)
        print(paste("Scaled logV for age",sub_group$AGE[j],":"))
        print(CI_V)
      }
    } #end for all ages in this age group
    
    if(length(sub_group$AGE) == sub_counter){ #If parameter precision inall pediatric subgroups have been reached
      
      print(paste("Stopping criteria reached for sub population",i)) 
      
    }else{
      print(paste("Stopping criteria NOT reached for population",i)) 
    }
    
    
  }###end for all age groups##############
  
  if(stop_MBAOD == length(AGE_WT[,1])){ #If parameter precision inall pediatric subgroups have been reached
    cat("\n")
    print("Stopping criteria reached for all groups, stopping MBAOD") 
    return(TRUE)
  }else{
    cat("\n")
    print("Stopping criteria NOT reached for all groups, starting next cohort")
    return(FALSE)
  }
}



#Option 4
#Bootstrap the simulated data and get the CI and geom. mean for the individually scaled
#parameters and check that they are within the range of the geometrical mean

print("--------:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::-------")
}else{
  stop_mbaod <- FALSE
  new_xspace <- NULL
  return(list(stop_mbaod,new_xspace))
}

}


gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}