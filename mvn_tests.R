covmat <- read.table("est.cov"
                     ,skip=1,header = T,row.names=1,check.names=T)


thetas <- c(1.04,19.9,2.19,3.66)
omegas <- c(0.0528,0,0.0533)
sigmas <- c(0.0149,0,0.0001)

library(stats)
library(mvtnorm)
nsub=11
nsim=1000
CL_thetas=c(1,3,4)
V_thetas=c(2)
thetas <- c(paste("THETA",sort(c(CL_thetas,V_thetas)),sep=""))
covmat <- covmat[rownames(covmat) %in% thetas,colnames(covmat) %in% thetas]
covmat <- as.matrix(covmat)

sigma2 <-cov2cor(covmat)
Sigma = covmat/(nsub-6/(nsub-6-2))

params <- data.frame(rmvt(n=nsim, df=nsub-6, delta=thetas, sigma=covmat))
params2 <- data.frame(rmvnorm(n=nsim, mean=c(thetas,sigmas,omegas), sigma=covmat))
params3 <- data.frame(rmvt(n=nsim, df=nsub-6, delta=thetas, sigma=Sigma))



df <- size_mat_scaling(params,age=1300,wt=78)



df <- df[is.finite(df$LCL) & is.finite(df$LV), ]


cohort_res <- aod_res$cohort_2
omegas <- unlist(cohort_res$est_result$omega)
sigmas <- unlist(cohort_res$est_result$sigma)
mean=c(1.016,2.983,49.79,8.279,omegas,sigmas)




cohort_res <- aod_res$cohort_3


