##Stopping function for MBAOD bridging study with 2 scaling parameters (Fixed WT hill on CL)
stop_dose <- function(cohort_num,cohort_res,design.db=NULL,opt_setting=NULL,cohort_dir="",nsim=1000,samples=100,ci_corr=0,
                      alpha=0.05,sim_eff=PKPD.ff,sim_params=PKPD.fg,
                      power=F,lower_limit=0.6,higher_limit=1.4, use_FIM=F){
  library(mvtnorm)
  library(reshape)
  print(paste("--::::Checking Stopping Criteria for Cohort", cohort_num,"::::--"))
  
  ###This block calculates the total number of individuals in each age group
  
  if(cohort_num==1){
    xt = cohort_res$design$xt
    dose    <- unname(unlist(cohort_res$design$a))
    
    if(length(dose)>length(cohort_res$design$groupsize)){
      gs      <- rep(cohort_res$design$groupsize,length(dose))
    }else{
      gs <- cohort_res$design$groupsize
    }
    
  }else{ 
    xt = cohort_res$design$xt[1,]
    gs      <- cohort_res$opt_result$opt_output$poped.db$design$groupsize
    dose    <- cohort_res$opt_result$opt_output$poped.db$design$a
    
  }
  nsub    <- sum(gs)  
  #setup the vector keeping track of groups which has passed  
  group_pass  <- rep(0,length(gs)) 
  
  
  #Storing all calculated CIs in vectors for plotting later
  Eff_CIs <- matrix(0,nrow=length(gs),ncol=2+2*length(xt))
  colnames(Eff_CIs) <- c("Group", "Dose", paste(c("Lower t=","Higher t="),rep(xt,each=2)))  
  
  k=1 #to keep track of vector posistion 
  thetas <- cohort_res$est_result$thetas 
  omegas <- unlist(cohort_res$est_result$omega) 
  sigmas <- unlist(cohort_res$est_result$sigma) 
  
  thetas <- thetas[thetas!=0]
  sigmas <- sigmas[sigmas>0]
  omegas <- omegas[omegas!=0]
  
  
  
  not_fixed_thetas <- sum(cohort_res$opt_result$opt_output$poped.db$parameters$notfixed_bpop)
  not_fixed_etas <- sum(cohort_res$opt_result$opt_output$poped.db$parameters$notfixed_d)
  not_fixed_sigmas <- sum(cohort_res$opt_result$opt_output$poped.db$parameters$notfixed_sigma)
  
  df  <- nsub- not_fixed_thetas - (not_fixed_etas + not_fixed_sigmas)/2
  
  
  
  if(use_FIM==T){
    inv_FIM <- 1/cohort_res$opt_result$opt_output$fmf
    thetas<- cohort_res$opt_result$opt_output$poped.db$parameters$bpop[,2]
    covmat <- matrix(0,nrow=length(thetas),ncol=length(thetas))
    covmat[4:7,4:7] <- inv_FIM[1:sum(not_fixed_thetas),1:sum(not_fixed_thetas)]

    omegas <- cohort_res$opt_result$opt_output$poped.db$parameters$d[,2]
    sigmas <- diag(cohort_res$opt_result$opt_output$poped.db$parameters$sigma)
    
    
    cohort_res$est_result$cov_mat <- "inv_FIM"
  }
  
  if(use_FIM==F){
    covmat <- cohort_res$est_result$cov_mat
    theta_rows <- c(paste("THETA",1:7,sep=""))
    covmat <- covmat[rownames(covmat) %in% theta_rows,colnames(covmat) %in% theta_rows]
    covmat <- as.matrix(covmat)
    print(covmat)
  }else{
    #covmat <- inv_FIM[sort(c(CL_thetas,V_thetas)), sort(c(CL_thetas,V_thetas))]
    covmat <- as.matrix(covmat)
    print(covmat)
  }
  
  fail_flag<-F
  #   for (z in 1:length(covmat[1,])){
  #     if(sum(covmat[z,])==0){
  #       print(paste("Variance of Theta", z,"is zero, failing stopping criteria and adding new cohort" ))
  #       fail_flag<-T
  #     }
  #   }
  
  
  if(!is.null(cohort_res$est_result$cov_mat)  & fail_flag==F){  
    
    
    
    
    #   thetas <- c(paste("THETA",sort(c(CL_thetas,V_thetas)),sep=""))
    #   covmat <- covmat[rownames(covmat) %in% thetas,colnames(covmat) %in% thetas]
    #   covmat <- as.matrix(covmat)
    
    
    scale_mat  <- covmat*((df-2)/df)
    
    
    
    for(i in 1:length(dose)){
      cat("\n")
      print(paste("___: Dose arm",i,": ", dose[i],"mg :___"))
      ####For all children ages in that group
      
      if(ci_corr==1){ ##Bonferroni correction of CI
        alpha_corr <- alpha/(cohort_num-1)
        alpha <- alpha_corr
      }
      
      samp_CI <- matrix(0,nrow=samples,ncol=length(xt)*2)
      param_thetas <- sim_params(x=NULL,a=dose[i],bpop = thetas,b=rep(0,length(thetas)),bocc = NULL)
      typical_pred <- sim_eff(model_switch = NULL,xt = xt,parameters = param_thetas, poped.db=NULL)$y
      
      for(samp in 1:samples){
        
        params <- rmvt(n=nsim, df=df, delta=thetas,sigma=scale_mat)
        if(is.na(params[1,1])==T){
          params <- rmvt(n=nsim, df=df, delta=thetas,sigma=diag(diag(scale_mat),7,7))
        }
        eff_df <- matrix(0,nrow=nsim, ncol=length(xt))
        colnames(eff_df) <-c(xt) 
        
        for(j in 1:nsim){
          param_tmp <- sim_params(x=NULL,a=dose[i],bpop = params[j,],b=rep(0,length(thetas)),bocc = NULL)
          eff_df[j,] <- sim_eff(model_switch = NULL,xt = xt,parameters = param_tmp, poped.db=NULL)$y
        }
        eff_CI <- apply(eff_df,2,quantile, probs=c(alpha/2,1-alpha/2),na.rm=T)/matrix(rep(typical_pred,2),nrow = 2,byrow =T)
        samp_CI[samp,] <- melt(eff_CI)$value
        
      } 
      
      Eff_CIs[i,] <- c(i,dose[i],apply(samp_CI,2,median))
      
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
        
        
        if(sum(Eff_CIs[i,paste("Lower t=",xt)]>lower_limit)==length(xt)& 
           sum(Eff_CIs[i,paste("Higher t=",xt)]<higher_limit)==length(xt)){
          print(paste(" All Effect CIs within ",lower_limit*100,"%-",higher_limit*100,"%  of Typical Prediction: ",sep = ""))
          group_pass[i] <-1
        }else{
          print(paste(" All Effect CIs NOT within ",lower_limit*100,"%-",higher_limit*100,"%  of Typical Prediction: ",sep = ""))
          #print(Eff_CIs[i,])
          group_pass[i] <-0
        }
      }
      
      #end for all ages in this age group
      
    }###end for all age groups##############
    
    if(sum(group_pass) == length(gs)){ #If parameter precision inall pediatric subgroups have been reached
      cat("\n")
      print("Stopping criteria reached for all groups, stopping MBAOD") 
      stop_MBAOD <- TRUE
    }else{
      cat("\n")
      print("Stopping criteria NOT reached for all groups, starting next cohort")
      stop_MBAOD <- FALSE
      
    }
    ##Part to generate the true D-optimal design to calculate the efficiency later
    if (!is.null(cohort_res$opt_result)){
      true.db <- cohort_res$opt_result$opt_output$poped.db
      true.db$parameters$bpop[,2][4:7] <- c(1,100,7,2)
      true.db$parameters$d[,2][3:4] <- c(0.0625,0.0625)
      true.db$parameters$sigma[1,1] <- 0.015
      
      true.db$parameters$param.pt.val$bpop <- cbind(true.db$parameters$bpop[,2])
      true.db$parameters$param.pt.val$d[3,3] <-0.0625 
      true.db$parameters$param.pt.val$d[4,4] <-0.0625
      true.db$parameters$param.pt.val$sigma[1,1] <- 0.015
      
      
      design.db$parameters$bpop <- cbind(rep(0,7), c(0.15,8,1,1,100,7,2),rep(0,7))
      design.db$parameters$d[,2][3:4] <- c(0.0625,0.0625)
      design.db$parameters$sigma[1,1] <- 0.015
      
      n_p <- sum(not_fixed_thetas +not_fixed_etas + not_fixed_sigmas)
      opt_setting$d_switch <-1
      true_res <- do.call(poped_optim,
                          c(poped.db=list(design.db),
                            opt_setting,
                            out_file=file.path(cohort_dir,"PopED_ouput_true.txt")))

           
      guess_fim <- det(evaluate.fim(true.db,fim.calc.type = 1)) 
      true_fim  <- exp(true_res$ofv)#log(det(evaluate.fim(true_res$poped.db,fim.calc.type = 1))) 
      design_efficiency1 <-(log(guess_fim)/log(true_fim))^(1/n_p) 
      design_efficiency2 <- ofv_criterion(guess_fim, n_p, true.db)/ofv_criterion(true_fim, n_p, true.db)

      
    }else{ 
      design_efficiency1<-0
      design_efficiency2<-0
      true_res = NA
    }
    
    
    
    
    return(list(stop_MBAOD,Eff_CIs,design_efficiency1,design_efficiency2,true_res))
    
    
    
    print("--------:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::-------")
  }else{
    stop_mbaod <- FALSE
    return(list(stop_MBAOD,Eff_CIs,design_efficiency1,design_efficiency2,true_res))
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
