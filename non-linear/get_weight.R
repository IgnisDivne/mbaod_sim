#Function which takes the PMA and returns a weight from the distribution of weights 
#the PMA is associated with according to CDC DATA
get_weight <- function(ID,PMA, probs = c(0.025,0.975),median=FALSE){

weight_all <- 1:length(PMA)
PMA_all <- PMA

df <- data.frame(ID,PMA_all)
df_u <- unique(df)

PMA <- df_u$PMA
weight <-1: length(PMA)

  for (i in 1:length(PMA)){  
    #if 0<3 months
    if (PMA[i] < 53.05){
      weight[i] <- rtnorm(1,5.1,0.840,probs,median)
    
    } #3-<6 month
    else if (PMA[i] >= 53.05 & PMA[i] < 66.1){
      weight[i] <- rtnorm(1,6.906,0.788,probs,median)

    }#6-<9 month
    else if (PMA[i] >= 66.1 & PMA[i] < 79.15){
      weight[i] <- rtnorm(1,8.4,0.697,probs,median)
  
    }#9-<12 month
    else if (PMA[i] >= 79.15 & PMA[i] < 92.2){
      weight[i] <- rtnorm(1,9.29,0.878,probs,median)
     
    } #1 year
    else if (PMA[i] >= 92.2 & PMA[i] < 144.4){
      weight[i] <- rtnorm(1,11.11,1.181,probs,median)
      
    } #2 year
    else if (PMA[i] >= 144.4 & PMA[i] < 196.6){
      weight[i] <- rtnorm(1,13.72,1.333,probs,median)
   
    } #3 year
    else if (PMA[i] >= 196.6 & PMA[i] < 248.8){
      weight[i] <- rtnorm(1,15.96,1.666,probs,median)
    
    }  #4 year
    else if (PMA[i] >= 248.8 & PMA[i] < 301){
      weight[i] <- rtnorm(1,18.14,1.87,probs,median)
   
    }#5 year
    else if (PMA[i] >= 301 & PMA[i] < 353.3){
      weight[i] <- rtnorm(1,21.15,2.599,probs,median)
      
    }#6 year
    else if (PMA[i] >= 353.3 & PMA[i] < 405.5){
      weight[i] <- rtnorm(1,23.96,3.40,probs,median)
      
    }#7 year
    else if (PMA[i] >= 405.5 & PMA[i] < 457.7){
      weight[i] <- rtnorm(1,26.75,3.819,probs,median)
      
    }#8 year
    else if (PMA[i] >= 457.7 & PMA[i] < 509.9){
      weight[i] <- rtnorm(1,31.59,5.27,probs,median)
      
    }#9 year
    else if (PMA[i] >= 509.9 & PMA[i] < 562.1){
      weight[i] <- rtnorm(1,36.03,6.05,probs,median)
      
    }#10 year
    else if (PMA[i] >= 562.1 & PMA[i] < 614.3){
      weight[i] <- rtnorm(1,40.53,7.08,probs,median)
      
    }#11 year
    else if (PMA[i] >= 614.3 & PMA[i] < 666.5){
      weight[i] <- rtnorm(1,47.06,8.825,probs,median)
      
    }#12 year
    else if (PMA[i] >= 666.5 & PMA[i] < 718.7){
      weight[i] <- rtnorm(1,51.91,10.46,probs,median)
      
    }#13year
    else if (PMA[i] >= 718.7 & PMA[i] < 770.9){
      weight[i] <- rtnorm(1,58.02,9.67,probs,median)
      
    }#14year
    else if (PMA[i] >= 770.9 & PMA[i] <823.1 ){
      weight[i] <- rtnorm(1,62.78,11.76,probs,median)
      
    }#15year
    else if (PMA[i] >= 823.1 & PMA[i] < 875.3){
      weight[i] <- rtnorm(1,66.96,11.09,probs,median)
      
    }#16year
    else if (PMA[i] >= 875.3 & PMA[i] < 927.5){
      weight[i] <- rtnorm(1,69.11,10.16,probs,median)
      
    }#17year
    else if (PMA[i] >= 927.5 & PMA[i] < 979.7){
      weight[i] <- rtnorm(1,71.48,10.83,probs,median)
      
    }#18year
    else if (PMA[i] >= 979.7 & PMA[i] < 1031.9){
      weight[i] <- rtnorm(1,73.49,11.11,probs,median)
      
    }#20-29
    else if ( PMA[i] >= 1031.9){
      weight[i] <- rtnorm(1,78.367,12.76,probs,median)
     
    }
  }
  
df_u$weight <- weight
df_wt <- merge(df,df_u,by="ID")
df_wt <-df_wt[order(df_wt$ID),]
return(df_wt$weight)

}

get_median_weight <- function(PMA){
  
  weight <- 1:length(PMA)
  for (i in 1:length(PMA)){  
    #if 0<3 months
    if (PMA[i] < 53.05){
      weight[i] <- 5.1
      
    } #3<6 month
    else if (PMA[i] >= 53.05 & PMA[i] < 66.1){
      weight[i] <- 6.906
      
    }#6<9 month
    else if (PMA[i] >= 66.1 & PMA[i] < 79.15){
      weight[i] <- 8.4
      
    }#9<12 month
    else if (PMA[i] >= 79.15 & PMA[i] < 92.2){
      weight[i] <- 9.29
      
    } #1 year
    else if (PMA[i] >= 92.2 & PMA[i] < 144.4){
      weight[i] <- 11.11
      
    } #2 year
    else if (PMA[i] >= 144.4 & PMA[i] < 196.6){
      weight[i] <- 13.72
      
    } #3 year
    else if (PMA[i] >= 196.6 & PMA[i] < 248.8){
      weight[i] <- 15.96
      
    }  #4 year
    else if (PMA[i] >= 248.8 & PMA[i] < 301){
      weight[i] <- 18.14
      
    }#5 year
    else if (PMA[i] >= 301 & PMA[i] < 353.3){
      weight[i] <- 21.15
      
    }#6 year
    else if (PMA[i] >= 353.3 & PMA[i] < 405.5){
      weight[i] <- 23.96
      
    }#7 year
    else if (PMA[i] >= 405.5 & PMA[i] < 457.7){
      weight[i] <- 26.74
      
    }#8 year
    else if (PMA[i] >= 457.7 & PMA[i] < 509.9){
      weight[i] <- 31.59
      
    }#9 year
    else if (PMA[i] >= 509.9 & PMA[i] < 562.1){
      weight[i] <- 36.03
      
    }#10year
    else if (PMA[i] >= 562.1 & PMA[i] < 614.3){
      weight[i] <- 40.53
      
    }#11 year
    else if (PMA[i] >= 614.3 & PMA[i] < 666.5){
      weight[i] <- 47.06
      
    }#12year
    else if (PMA[i] >= 666.5 & PMA[i] < 718.8){
      weight[i] <- 51.908
      
    }##13year
    else if (PMA[i] >= 718.7 & PMA[i] < 770.9){
      weight[i] <- 58.02
      
    }#14year
    else if (PMA[i] >= 770.9 & PMA[i] <823.1 ){
      weight[i] <- 62.78
      
    }#15year
    else if (PMA[i] >= 823.1 & PMA[i] < 875.3){
      weight[i] <- 66.96
      
    }#16year
    else if (PMA[i] >= 875.3 & PMA[i] < 927.5){
      weight[i] <- 69.11
      
    }#17year
    else if (PMA[i] >= 927.5 & PMA[i] < 979.7){
      weight[i] <- 71.48
      
    }#18year
    else if (PMA[i] >= 979.7 & PMA[i] < 1031.9){
      weight[i] <- 73.49
      
    }#20-29
    else if ( PMA[i] >= 1031.9){
      weight[i] <- 78.367
      
    }
  }
  
  return(weight)
}

get_age_from_median_weight <- function(WT){
  
  AGE <- 1:length(WT)
  PMA <- 1:length(WT)
  
  for (i in 1:length(WT)){  
    #if 0<2 months
    if (WT[i] == 5.1){
      AGE[i] <- "0  to <2 months"
      PMA[i] <- 44.35
    } #2<6 month
    else if (WT[i] == 6.906){
      AGE[i] <- "3 to <5 months"
      PMA[i] <- 57.2
    }#6<9 month
    else if (WT[i]== 8.4){
      AGE[i] <- "6 to <12 months"
      PMA[i] <- 76.55
    }#9<12 month
    else if (WT[i] == 9.29){
      AGE[i] <- "6 to <12 months"
      PMA[i] <- 85.15
    } #1 year
    else if (WT[i] == 11.11){
      AGE[i] <- "1 to <2 years"
      PMA[i] <- 118.3
    } #2 year
    else if (WT[i] == 13.72){
      AGE[i] <- "2 to <6 years"
      PMA[i] <- 170.5
      
    } #3 year
    else if (WT[i] == 15.96){
      AGE[i] <- "2 to <6 years"
      PMA[i] <- 222.7
    }  #4 year
    else if (WT[i] == 18.14){
      AGE[i] <- "2 to <6 years"
      PMA[i] <- 274.9
      
    }#5 year
    else if (WT[i] == 21.15){
      AGE[i] <- "2 to <6 years"
      PMA[i] <- 327.1
      
    }#6 year
    else if (WT[i] == 23.96){
      AGE[i] <- "6 to <12 years"
      PMA[i] <- 379.3
      
    }#7 year
    else if (WT[i] == 26.74){
      AGE[i] <- "6 to <12 years"
      PMA[i] <- 431.5
      
    }#8 year
    else if (WT[i] == 31.59){
      AGE[i] <- "6 to <12 years"
      PMA[i] <- 483.7
    }#9 year
    else if (WT[i] == 36.03){
      AGE[i] <- "6 to <12 years"
      PMA[i] <- 535.9
      
    }#10year
    else if (WT[i] == 40.53){
      AGE[i] <- "6 to <12 years"
      PMA[i] <- 588.1
      
    }#11 year
    else if (WT[i] == 47.06){
      AGE[i] <- "6 to <12 years"
      PMA[i] <- 640.3
      
    }#12year
    else if (WT[i] == 51.908){
      AGE[i] <- "12 to 18 years"
      PMA[i] <- 692.5
      
    }#13 year
    else if (WT[i] ==   58.02){
      AGE[i] <- "12 to 18 years"
      PMA[i] <- 744.7
    
      }#14 year
    else if (WT[i] ==   62.78){
      AGE[i] <- "12 to 18 years"
      PMA[i] <- 744.7 +52.2
    
      }#15 year
    else if (WT[i] ==   66.96){
      AGE[i] <- "12 to 18 years"
      PMA[i] <- 744.7 + 52.2 + 52.2
    
      }#16 year
    else if (WT[i] ==   69.11){
      AGE[i] <- "12 to 18 years"
      PMA[i] <- 744.7 + 52.2 + 52.2 + 52.2
    
      }#17 year
    else if (WT[i] ==   71.48){
      AGE[i] <- "12 to 18 years"
      PMA[i] <- 744.7 + 4*52.2
    
    }#18 year
    else if (WT[i] ==   73.49){
      AGE[i] <- "12 to 18 years"
      PMA[i] <- 744.7 + 5*52.2
    }#20-29
    else if ( WT[i] == 78.367){
      AGE[i] <- "Adult"
      PMA[i] <- 1298.9
    }
  }
  AGE_df <- data.frame(AGE,PMA)
  return(AGE_df)
}

age_group_2_PMA <- function(ID,age_group,maxmin =F){
 # age_group <- age_group +1

  PMA_all <- 1:length(age_group)
  age_group_all <- age_group
  df <- data.frame(ID,age_group_all)
  df_u <- unique(df)
  age_group <- df_u$age_group

  
 
  #age_group=0          1       2        3    4     5     6    7
           #default, 3-<6mo,6-<12mo, 1-<2y,2-<6y,6-<12y,12-<18,20-29y
  ages_min <- c(53.05,66.1,92.2,144.4,353.3,666.5,1031.4)
  ages_max <- c(ages_min[2:length(ages_min)], 1513.8)
  df_u$PMA <- runif(length(age_group),ages_min[age_group],ages_max[age_group])
  df_out <- merge(df,df_u,by="ID")
  df_out <- df_out[order(df_out$ID),]
  return(df_out$PMA)
}
  
age_group_2_mean_PMA <- function(age_group){
  #age_group=0       1       2        3    4     5     6      7
          #default,3-<6mo,6-<12mo, 1-<2y,2-<6y,6-<12y,12-19y20-29y
age_group <- age_group +1

  ages_min <- c(1031.4,53.05,66.1,92.2,144.4,353.3,666.5,1031.4)
  ages_max <- c(ages_min[2:length(ages_min)], 1513.8)
  
  ages_mean <- rowMeans(data.frame(ages_min,ages_max))
  PMA <- ages_mean[age_group]

  return(PMA)
}



rtnorm <- function(n, mean, sd, probs,median=FALSE){
  if(median == FALSE){
    a <- quantile(rnorm(10000,mean,sd),probs = probs)
    round(qnorm(runif(n, pnorm(a[1], mean, sd), pnorm(a[2], mean, sd)), mean, sd),digits=1)
  }else{
    round(quantile(rnorm(100000,mean,sd),probs=0.5),digits=1)
  }
  
}



get_age_group_from_median_weight <- function(WT){
  
  AGE<- 1:length(WT)

  
  for (i in 1:length(WT)){  

    #2<6 month
    if(WT[i] == 6.906){
      AGE[i] <- 1
    }#6<9 month
    else if (WT[i]== 8.4){
      AGE[i] <- 2
    }#9<12 month
    else if (WT[i] == 9.29){
      AGE[i] <- 2
    } #1 year
    else if (WT[i] == 11.11){
      AGE[i] <- 3
    } #2 year
    else if (WT[i] == 13.72){
      AGE[i] <- 4 
    } #3 year
    else if (WT[i] == 15.96){
      AGE[i] <- 4
    }  #4 year
    else if (WT[i] == 18.14){
      AGE[i] <- 4
    }#5 year
    else if (WT[i] == 21.15){
      AGE[i] <- 4
    }#6 year
    else if (WT[i] == 23.96){
      AGE[i] <- 5
    }#7 year
    else if (WT[i] == 26.74){
      AGE[i] <- 5
    }#8 year
    else if (WT[i] == 31.59){
      AGE[i] <- 5
    }#9 year
    else if (WT[i] == 36.03){
      AGE[i] <- 5
    }#10year
    else if (WT[i] == 40.53){
      AGE[i] <- 5
    }#11 year
    else if (WT[i] == 47.06){
      AGE[i] <- 5
    }#12year
    else if (WT[i] == 51.908){
      AGE[i] <- 6
    }#13 year
    else if (WT[i] ==   58.02){
      AGE[i] <- 6
    }#14 year
    else if (WT[i] ==   62.78){
      AGE[i] <- 6
    }#15 year
    else if (WT[i] ==   66.96){
      AGE[i] <- 6
    }#16 year
    else if (WT[i] ==   69.11){
      AGE[i] <- 6
    }#17 year
    else if (WT[i] ==   71.48){
      AGE[i] <- 6
    }#18 year
    else if (WT[i] ==   73.49){
      AGE[i] <- 6
    }#20-29
    else if ( WT[i] == 78.367){
      AGE[i] <- 7
    }
  }
  return(AGE)
}


                     