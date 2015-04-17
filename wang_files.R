f<-function(s,n,sigma){2*((n-1)/2/sigma**2)**((n-1)/2)/
                         gamma(0.5*(n-1))*s**(n-2)*exp(-(n-1)*s**2/2/sigma**2)}
nmin<-4
nmax<-25
nsub<-nmax-nmin+1
result<-rep(0,nsub)
sd<-0.4 #standard deviation of logX from adult data
for (i in 1:nsub) {
  n<-nmin+i-1
  tv<-qt(0.975, n-1)
  sup<-log(1.4)*sqrt(n)/tv
  g<-function(x){f(x,n,sd)}
  poweri<-integrate(g, 0, sup, subdivisions=100)
  result[i]<-poweri[[1]]
}
power<-data.frame(nsub=c(nmin:nmax), power=result)
power




#Assume covariate model CL=theta1*(wt/70)** theta2*age/(age+theta3) and
#Log of CL was modelled as LCL=thetap1+theta2*(wt/70)+log(age/(age+theta3))
#input parameter estimates
thetap1=3.7421
theta2=1.0078
theta3=4.8422
#input variance-covariance matrix for the 3 parameters above
covp3=matrix(
  c(
    0.29810, 0.05782, 1.27120,
    0.05782, 0.02921, 0.02073,
    1.27120, 0.02073, 8.42210
  ),nrow=3, byrow=T
)
#define the weight and age combination for SE estimation
wt=50
age=12
#define covariate model in R
f <- function(x,y,z) x + y*log(wt/70) + log(age/(age+z))
df<- deriv(body(f), c("x","y","z"))
x=thetap1
y=theta2
z=theta3
out=eval(df)
dfout=attr(out,"gradient")
varlcl=dfout%*%covp3%*%t(dfout)
#SE of LCL 0.09436884
SElcl=sqrt(varlcl)


c(exp(-qt(0.975,1000)*SElcl),exp(qt(0.975,1000)*SElcl))

# rep_limit<- function(j=2){
#   
#   if(j >1){
#     return(TRUE)
#   }else{
#     return(FALSE)
#   }
#   
# }