library(PopED)

PK.1.comp.maturation.fg <- function(x,a,bpop,b,bocc){
  ## -- parameter definition function

  parameters=c( CL=exp(bpop[1])*exp(b[1]),
                V=exp(bpop[2])*exp(b[2]),
                TM50=exp(bpop[3]),
                HILL=exp(bpop[4]),
                age_group = a[1],
                PMA = age_group_2_mean_PMA(a[1]),
                WT = get_median_weight(age_group_2_mean_PMA(a[1])))
  return( parameters )
}


PK.1.comp.maturation.ff <- function(model_switch,xt,parameters,poped.db){
  with(as.list(parameters),{
    y=xt

    CL=CL*(WT/70)^0.75*((PMA^HILL)/(PMA^HILL + TM50^HILL))
    V=V*(WT/70)


    DOSE=1000*(WT/70)
    y = DOSE/V*exp(-CL/V*xt)
       
    return(list( y= y,poped.db=poped.db))
  })
}

