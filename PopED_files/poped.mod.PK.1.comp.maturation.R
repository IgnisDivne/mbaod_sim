library(PopED)

PK.1.comp.maturation.fg <- function(x,a,bpop,b,bocc){
  ## -- parameter definition function 
  parameters=c( CL=bpop[1]*exp(b[1]),
                V=bpop[2]*exp(b[2]),
                TM50=bpop[3],
                HILL=bpop[4],
                PMA = a[1],
                WT = a[2])
  return( parameters )
}


PK.1.comp.maturation.ff <- function(model_switch,xt,parameters,poped.db){
  with(as.list(parameters),{
    y=xt
    
    CL=CL*(WT**0.75)*((1/70)**0.75)*(PMA**HILL)/(PMA**HILL + TM50**HILL)
    V=V*(WT/70)
    DOSE=1000*(WT/70)
    y = DOSE/V*exp(-CL/V*xt) 
       
    return(list( y= y,poped.db=poped.db))
  })
}

