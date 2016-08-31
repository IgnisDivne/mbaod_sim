library(PopED)

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
   # MS <- model_switch
    y=xt
    Ke=CL/V
    
    C = ((Dose*Ka*Ke)/(CL*(Ka-Ke)))*(exp(-Ke*xt)-exp(-Ka*xt))
    y = Base+((C^hill*EMAX)/(EC50^hill + C^hill))
    
    #y[MS==1] = C[MS==1]
   # y[MS==2] = EFF[MS==2]
    
    return(list( y= y,poped.db=poped.db))
  })
}



feps <- function(model_switch,xt,parameters,epsi,poped.db){
  ## -- Residual Error function
  ## -- Proportional PK + additive PD
  returnArgs <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db)) 
  y <- returnArgs[[1]]
  poped.db <- returnArgs[[2]]
  
 # MS <- model_switch
  
  y <- y*(1+epsi[,1])+epsi[,2]
  #add.err <- y
  
  
 # y[MS==1] = add.err[MS==1]
 # y[MS==2] = prop.err[MS==2]
  
  return(list( y= y,poped.db =poped.db )) 
}


