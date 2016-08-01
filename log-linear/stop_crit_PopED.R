##Stopping function for MBAOD bridging study for use with regular OD where the uncertianty calculations
##are based on the predicted FIM
stop_crit_PopED <- function(cohort_num,cohort_res,option=2,nsim=10000, CL_thetas=c(1,3),V_thetas=2,sim_params=loglin_size_mat_scaling, age_space = c(1,2,3,4,5,6),allow_next=2,power=F){
  print(paste("--::::Checking Stopping Criteria for Cohort", cohort_num,"::::--"))



  library(mvtnorm)

  
  ###This block calculates the total number of individuals in each age group
  gs      <- cohort_res$opt_result$opt_output$poped.db$design$groupsize
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
  
  #Storing all calculated powers in vectors for plotting later
  v_p_CL <- matrix(nrow=length(AGE_WT$AGE),ncol=3)
  colnames(v_p_CL) <- c("Age Group","PMA","Power")
  v_p_V <- v_p_CL
  k=1 #to keep track of vector posistion
  
  fim_mat <- cohort_res$opt_result$opt_output$fmf
  thetas<- cohort_res$opt_result$opt_output$poped.db$design$bpop[,2]
  omegas <- cohort_res$opt_result$opt_output$poped.db$design$d
  sigmas <- diag(cohort_res$opt_result$opt_output$poped.db$design$sigma)
  
  
  if(cohort_num>1){  
    
    #Option 1 - Simulation based
    #Simulatea large number of individuals with uncertainty based
    #on the latest design and parameter estimates and parameter uncertainty
    #Get parameters used for each individual and calculate geometrical mean and 95% CI
    #for each parameter
    #Compares the 95% CI to the range of 60% of 140% of geometrical mean for each parameter
    #If the CI is within the mean range for all parameters and power is 80% or more,
    #the function stops the MBAOD after this cohort. 
    
    if(option==1){    
      #   thetas <- c(paste("THETA",sort(c(CL_thetas,V_thetas)),sep=""))
      #   covmat <- covmat[rownames(covmat) %in% thetas,colnames(covmat) %in% thetas]
      #   covmat <- as.matrix(covmat)
      covmat <- 1/fim_mat
      covmat[covmat==Inf] <- 0
      covmat <- covmat[rowSums(covmat)!=0,colSums(covmat)!=0]
      thetas <- thetas[thetas!=0]
      sigmas <- sigmas[sigmas>0.0001]
      omegas <- omegas[omegas!=0]
      covmat <- as.matrix(covmat)
      
      print(covmat)
      
      print(c(thetas,sigmas,omegas))
      #For all pediatric groups, group 1 is adults, hence i in 2:length(AGE_WT)
      params <- rmvnorm(n=nsim, mean=c(thetas,sigmas,omegas),sigma=covmat)#,omegas,sigmas), sigma=covmat)
      for(i in 1:length(age_space)){
        cat("\n")
        print(paste("___ Pediatric Sub Population",i,": Children of age",age_groups$AGE[age_space[i]],"___"))
        ####For all children ages in that group
        sub_group <- subset(AGE_WT,AGE_WT$AGE==age_groups$AGE[age_space[i]])
        sub_counter <- 0
        
        for(j in 1:length(sub_group$AGE)){
          
          wt=sub_group$WT[j]
          age=sub_group$PMA[j]
          params_df <- sim_params(params,age,wt)
          
          SElcl=sd(params_df$LCL)
          SELV=sd(params_df$LV)
          
          
          if(is.na(SElcl)){
            print(paste("Clearance for Children of age",sub_group$AGE[j],"approaches Zero due to parameter estimates or SE of estimates"))
            
            print(paste("Thetas:",cohort_res$est_result$thetas))
            print(paste("SE Thetas:",cohort_res$est_result$sethetas))
          }else if(is.na(SELV)){
            print(paste("Volume of Distribution for Children of age",sub_group$AGE[j],"approaches Zero due to parameter estimates or SE of estimates"))
            print(paste("Thetas:",cohort_res$est_result$thetas))
            print(paste("SE Thetas:",cohort_res$est_result$sethetas))
          }else{
            
            CI_CL=c(exp(-qt(0.975,nsub-6)*SElcl),exp(qt(0.975,nsub-6)*SElcl))
            CI_V=c(exp(-qt(0.975,nsub-6)*SELV),exp(qt(0.975,nsub-6)*SELV))
            p_CL <- get_power_crit(SElcl, nsub-6)
            p_V <- get_power_crit(SELV, nsub-6)
            
            if(power==T){
              if(p_CL>=0.8 & p_V > 0.8){
                sub_counter <- sub_counter +1
              }else{
                print(" All Parameters Power of CI not >=80%: ")
                print(paste("Power of Scaled logCL for age",sub_group$AGE[j],":",p_CL*10,"%"))
                print(paste("Power of Scaled logV for age",sub_group$AGE[j],":",p_V*10,"%"))
                
              }
              
            }else{
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
            }
            
            CL_CIs[k,] <- c(i,age,CI_CL)
            V_CIs[k,] <- c(i,age,CI_V)
            v_p_CL[k,] <-c(i,age,p_CL)
            v_p_V[k,] <-c(i,age,p_V)
            
            k=k+1
          }
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
      if(allow_next == T & sum(xspace)>0) xspace[min(xspace[xspace>1])-1] <- min(xspace[xspace>1])-1
      
      new_xspace <- unique(c(xspace[xspace != 0],6,7))
      return(list(stop_MBAOD,new_xspace,CL_CIs,V_CIs))
      
    }
    
    
    #Option 2
    #From the estimates of thetas, etas and eps and their respective standard errors:
    #Calculate the 2.5th 92.5th percentiles of what is possible as described by Wang et al.
    #Compare the "CI" to 60% and 140% of geometrical mean 
    #(Theta?)
    if(option==2){
      ##Subset the COVMAT for delta9 method
      covmat <- 1/fim_mat
      covmat[covmat==Inf] <- 0
      print(covmat)
      #lcl_thetas <- c(paste("THETA",CL_thetas,sep=""))
      cov_lcl <- covmat[CL_thetas,CL_thetas]
      cov_lcl <- as.matrix(cov_lcl)
      
      #For all pediatric groups, group 1 is adults, hence i in 2:length(AGE_WT)
      for(i in 1:length(age_space)){
        cat("\n")
        print(paste("___ Pediatric Sub Population",i,": Children of age",age_groups$AGE[age_space[i]],"___"))
        ####For all children ages in that group
        sub_group <- subset(AGE_WT,AGE_WT$AGE==age_groups$AGE[age_space[i]])
        sub_counter <- 0
        for(j in 1:length(sub_group$AGE)){
          
          ##scaled LogCL SE calculation 
          base=cohort_res$est_result$thetas[1] #log
          TM50=cohort_res$est_result$thetas[3]
          HILL=0.75#thetaHILL
          
          SE_V <- cohort_res$est_result$sethetas[2]
          
          #input variance-covariance matrix for the 3 parameters above
          
          #define the weight and age combination for SE estimation
          wt=sub_group$WT[j]
          age=sub_group$PMA[j]
          
          #define covariate model in R
          
          f <- function(x,y,z) x + log(wt/70)*z + log(age/(age+y))
          df<- deriv(body(f), c("x","y"))
          x=base
          y=TM50
          z=HILL
          
          
          out=eval(df)
          dfout=attr(out,"gradient")
          varlcl=dfout%*%cov_lcl%*%t(dfout)
          SElcl=sqrt(varlcl)
          
          #df <- nsub -num_params. is num_params 8? for this model 4 thetas, 2 omegas 2 sigmas
          CI_CL=c(exp(-qt(0.975,nsub-6)*SElcl),exp(qt(0.975,nsub-6)*SElcl))
          CI_V=c(exp(-qt(0.975,nsub-6)*SE_V),exp(qt(0.975,nsub-6)*SE_V))
          
          
          

          
          if(power==T){
            p_CL <- get_power_crit(SElcl, nsub-6)
            p_V <- get_power_crit(SE_V, nsub-6)
            if(p_CL>=0.8 & p_V > 0.8){
              sub_counter <- sub_counter +1
            }else{
              print(" All Parameters Power of CI not >=80%: ")
              print(paste("Power of Scaled logCL for age",sub_group$AGE[j],":",p_CL*10,"%"))
              print(paste("Power of Scaled logV for age",sub_group$AGE[j],":",p_V*10,"%"))
            }
            v_p_CL[k,] <-c(i,age,p_CL)
            v_p_V[k,] <-c(i,age,p_V)
          }else{
            
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
      if(allow_next == 1 & sum(xspace)>0) xspace[min(xspace[xspace>1])-1] <- min(xspace[xspace>1])-1
      
      if(allow_next==2) xspace <- age_space
      
      new_xspace <- unique(c(xspace[xspace != 0],6,7))
      return(list(stop_MBAOD,new_xspace,CL_CIs,V_CIs))
      
    }
    
    
    print("--------:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::-------")
  }else{
    stop_mbaod <- FALSE
    new_xspace <- unlist(cohort_res$opt_result$opt_output$poped.db$design$discrete_x[cohort_num])
    return(list(stop_mbaod,new_xspace))
  }
  
}


gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}



get_power_crit <- function(sdlogX,nsub){
  f<-function(s,n,sigma){
    2*((n-1)/2/sigma**2)**((n-1)/2)/gamma(0.5*(n-1))*s**(n-2)*exp(-(n-1)*s**2/2/sigma**2)}
  result<-rep(0,nsub)
  sd<-sdlogX #standard deviation of logX from adult data
  
  n<-nsub
  tv<-qt(0.975, n)
  sup<-log(1.4)*sqrt(n)/tv
  g<-function(x){f(x,n,sd)}
  poweri<-integrate(g, 0, sup, subdivisions=100)
  power<-poweri[[1]]
  
  return(power)
}

loglin_size_mat_scaling <- function(params,age,wt){
  
  base<- params[,1]#log
  TM50<- params[,3]
  HILL<-0.75#thetaHILL
  V <- params[,2]
  
  LCL <-  base + log(wt/70)*HILL + log(age)-log(age+TM50)
  LV  <-  V + log(wt/70)*1
  
  df <- data.frame(LCL,LV)
  df <- df[is.finite(df$LCL) & is.finite(df$LV), ]
  return(df)
  
  
}

size_mat_scaling <- function(params,age,wt){
  
  base<- params[,1]
  TM50<- params[,3]
  HILL<-params[,4]
  V <- params[,2]
  
  CL <-  base * ((wt/70)^0.75) * (age^HILL)/(age^HILL+TM50^HILL)
  V  <-  V *(wt/70)
  
  df <- data.frame(LCL = CL,LV=V)
  df <- log(df[is.finite(df$LCL) & is.finite(df$LV), ])
  
  return(df)
}

emax <- function(params,age,wt){
  
  base<- params[,1]#log
  V <- params[,2]
  EMAX<- params[,3]
  EC50<- params[,4]
  HILL<-params[,5]
  
  
  CL=base+(EMAX*wt**HILL)/(EC50**HILL+wt**HILL)
  V=V*(wt/70)
  
  return(data.frame(LCL = log(CL),LV=log(V)))
}