##Stopping function for MBAOD bridging study with 2 scaling parameters (Fixed WT hill on CL)
stop_critX_2 <- function(cohort_num,cohort_res,option=3,nsim=1000, CL_thetas=c(1,3,4),V_thetas=2,ci_corr=0,
                         alpha=0.05,sim_params=size_mat_scaling, age_space = c(1,2,3,4,5,6),allow_next=T,
                         power=F,lower_limit=0.6,higher_limit=1.4, use_FIM=F){
  library(mvtnorm)
  print(paste("--::::Checking Stopping Criteria for Cohort", cohort_num,"::::--"))
  
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
  thetas <- cohort_res$est_result$thetas
  omegas <- unlist(cohort_res$est_result$omega)
  sigmas <- unlist(cohort_res$est_result$sigma)
  
  thetas <- thetas
  sigmas <- sigmas[sigmas>0.0001]
  omegas <- omegas[omegas!=0]
  
  
  df  <- nsub-6
  
  print(use_FIM)
  
  if(use_FIM==T){
    inv_FIM <- 1/cohort_res$opt_result$opt_output$fmf
    
    thetas<- cohort_res$opt_result$opt_output$poped.db$parameters$bpop[,2]
    print(thetas)
    omegas <- cohort_res$opt_result$opt_output$poped.db$parameters$d[,2]
    sigmas <- diag(cohort_res$opt_result$opt_output$poped.db$parameters$sigma)
    
    
    cohort_res$est_result$cov_mat <- "inv_FIM"
    
    
  }
  if(!is.null(cohort_res$est_result$cov_mat)  & cohort_num>1){  
    
    if(use_FIM==F){
      covmat <- cohort_res$est_result$cov_mat
      theta_rows <- c(paste("THETA",sort(c(CL_thetas,V_thetas)),sep=""))
      covmat <- covmat[rownames(covmat) %in% theta_rows,colnames(covmat) %in% theta_rows]
      covmat <- as.matrix(covmat)
      print(covmat)
    }else{
      covmat <- inv_FIM[sort(c(CL_thetas,V_thetas)), sort(c(CL_thetas,V_thetas))]
      covmat <- as.matrix(covmat)
      print(covmat)
    }  
    
    
  for (z in 2:length(covmat[1,])){
    if(sum(covmat[z,])==0){
      print(paste("Variance of Theta", z,"in the COV is exactly zero, Failing Stopping criterion due to estimation issues and adding another cohort." ))
      cohort_res$est_result$cov_mat <- NULL
    }
    
  }
  
}
  
  
  if(!is.null(cohort_res$est_result$cov_mat)  & cohort_num>1){  
    

    
    print(paste("You have selected option ",option,".",sep=""))
    
    
    
    
    if(option==3){    
      #   thetas <- c(paste("THETA",sort(c(CL_thetas,V_thetas)),sep=""))
      #   covmat <- covmat[rownames(covmat) %in% thetas,colnames(covmat) %in% thetas]
      #   covmat <- as.matrix(covmat)
      
      
      
      print("Entering Option 3")
      

      scale_mat  <- covmat*((df-2)/df)
      
      print(scale_mat)
      #For all pediatric groups, group 1 is adults, hence i in 2:length(AGE_WT)
      
      CLCI <- data.frame(low=1:1000,high=1:1000)    
      VCI <- data.frame(low=1:1000,high=1:1000) 
      for(i in 1:length(age_space)){
        cat("\n")
        print(paste("___ Pediatric Sub Population",i,": Children of age",age_groups$AGE[age_space[i]],"___"))
        ####For all children ages in that group
        sub_group <- subset(AGE_WT,AGE_WT$AGE==age_groups$AGE[age_space[i]])
        sub_counter <- 0
        
        if(ci_corr==1){ ##Bonferroni correction of CI
          alpha_corr <- alpha/(cohort_num-1)
          alpha <- alpha_corr
        }
        for(j in 1:length(sub_group$AGE)){
          
          wt=sub_group$WT[j]
          age=sub_group$PMA[j]
          for(samp in 1:1000){
            params <- rmvt(n=nsim, df=df, delta=thetas,sigma=scale_mat)
            params_df <- sim_params(params,age,wt)
            
            CLCI[samp,] <- quantile(params_df$CL, probs=c(alpha/2,1-alpha/2),na.rm=T)
            CLCI[samp,] <- CLCI[samp,]/gm_mean(params_df$CL)
            
            VCI[samp,] <- quantile(params_df$V, probs=c(alpha/2,1-alpha/2),na.rm=T)
            VCI[samp,] <- VCI[samp,]/gm_mean(params_df$V)
            
          }
          CL_CI <- c(median(CLCI[,1],na.rm = T), median(CLCI[,2],na.rm = T))
          
          V_CI <- c(median(VCI[,1],na.rm = T), median(VCI[,2],na.rm = T))
          
          if(is.na(CL_CI[1]) | is.na(CL_CI[2]) | CL_CI[1]==Inf| CL_CI[2]==Inf){
            print(paste("Clearance for Children of age",sub_group$AGE[j],"approaches Zero due to parameter estimates or SE of estimates"))
            
            print(paste("Thetas:",cohort_res$est_result$thetas))
            print(paste("SE Thetas:",cohort_res$est_result$sethetas))
          }else if(is.na(V_CI[1]) | is.na(V_CI[2]) | V_CI[1]==Inf| V_CI[2]==Inf){
            print(paste("Volume of Distribution for Children of age",sub_group$AGE[j],"approaches Zero due to parameter estimates or SE of estimates"))
            print(paste("Thetas:",cohort_res$est_result$thetas))
            print(paste("SE Thetas:",cohort_res$est_result$sethetas))
          }else{
            CI_CL <- CL_CI
            CI_V <- V_CI
            
            if(power==T){
              p_CL <- get_power_crit(SElcl, df)
              p_V <- get_power_crit(SELV, df)
              
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
              if(CI_CL[1] > lower_limit & CI_CL[2] < higher_limit 
                 & CI_V[1] > lower_limit & CI_V[2] < higher_limit){
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
      if(allow_next == T & sum(xspace)==0) xspace <- c(6)
      new_xspace <- unique(c(xspace[xspace != 0],6))
      return(list(stop_MBAOD,new_xspace,CL_CIs,V_CIs))
      
    }
    
    print("--------:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::-------")
  }else{
    stop_mbaod <- FALSE
    new_xspace <- unlist(cohort_res$opt_result$opt_output$poped.db$design_space$discrete_x[length(cohort_res$opt_result$opt_output$poped.db$design_space$discrete_x)])
    new_xspace <- unique(c(new_xspace,6))
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
  #params <- subset(params, params[,3]>0)
  
  
  #etacl<-params[,4]
  #etav<-params[,5]
  #TM50[TM50<0] <- 0
  
  base<- params[,1]# + etacl#log
  TM50<- exp(params[,3])
  
  V <- params[,2] #+ etav
  
  age <- rep(age,length(base))
  wt <- rep(wt,length(base))
  
  Hill=params[,4]
  
  LCL <-  base + log(wt/70)*0.75 + log(age^Hill/(age^Hill+TM50^Hill))
  LV  <-  V + log(wt/70)*1
  
  df <- data.frame(LCL,LV)
  
  return(df)
  
}

size_mat_scaling <- function(params,age,wt){
  
  base<-exp(params[,1])
  TM50<- exp(params[,3])
  HILL<-exp(params[,4])
  V <- exp(params[,2])
  
  
  
  CL <-  base * (wt/70)^0.75 * (age^HILL)/(age^HILL+TM50^HILL)
  V  <-  V *(wt/70)
  
  df <- data.frame(CL,V)
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























