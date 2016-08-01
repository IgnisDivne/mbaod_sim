library(PopED)

PK.1.comp.maturation.fg <- function(x,a,bpop,b,bocc){
  ## -- parameter definition function

  parameters=c( CL=bpop[1]+(b[1]),
                V=bpop[2]+(b[2]),
                TM50=bpop[3],
                #HILL=bpop[4],
                age_group = x[1],
                PMA = age_group_2_mean_PMA(x[1]),
                WT = get_median_weight(age_group_2_mean_PMA(x[1])))
  return( parameters )
}


PK.1.comp.maturation.ff <- function(model_switch,xt,parameters,poped.db){
  with(as.list(parameters),{
    y=xt
    HILL <- 0.75

    LCL=CL+HILL*log(WT/70)+log((PMA)/(PMA + TM50))
    LV=V+log(WT/70)


    DOSE=1000*(WT/70)
    V = exp(LV)
    CL = exp(LCL)
    y = DOSE/V*exp(-CL/V*xt)
       
    return(list( y= y,poped.db=poped.db))
  })
}

