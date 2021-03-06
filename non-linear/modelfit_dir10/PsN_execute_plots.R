#START OF AUTO-GENERATED PREAMBLE, WILL BE OVERWRITTEN WHEN THIS FILE IS USED AS A TEMPLATE
#Created 2016-09-01 at 11:56

rplots.level <- 1
xpose.runno <- ''
toolname <- 'execute'
pdf.filename <- paste0('PsN_',toolname,'_plots.pdf')
pdf.title <- 'execute diagnostic plots est.mod'
working.directory<-'C:/Users/erist836/Documents/GitHub/mbaod_sim/non-linear/modelfit_dir10/'
raw.results.file <- 'raw_results_est.csv'
model.directory<-'C:/Users/erist836/Documents/GitHub/mbaod_sim/non-linear/propofol_run_run_dir_9/rep_1/cohort_2/'
model.filename<-'est.mod'
subset.variable<-NULL
mod.suffix <- '.mod'
mod.prefix <- 'est'
tab.suffix <- ''
theta.labels <- c('TVCL','TVV','TM50','Maturation HILL')
theta.fixed <- c(FALSE,FALSE,FALSE,FALSE)
omega.labels <- c('OMEGA(1,1)','OMEGA(2,2)')
omega.fixed <- c(FALSE,FALSE)
sigma.labels <- c('SIGMA(1,1)','SIGMA(2,2)')
sigma.fixed <- c(FALSE,TRUE)
n.eta <- 2
n.eps <- 2

pdf.filename <- paste0(mod.prefix,xpose.runno,'_plots.pdf')

setwd(working.directory)

############################################################################
#END OF AUTO-GENERATED PREAMBLE
#WHEN THIS FILE IS USED AS A TEMPLATE THIS LINE MUST LOOK EXACTLY LIKE THIS



library(xpose4)

xpdb<-xpose.data(xpose.runno,directory=model.directory,tab.suffix=tab.suffix,mod.prefix=mod.prefix,mod.suffix=mod.suffix)

pdf(file=pdf.filename,width=10,height=7,title=pdf.title)

#uncomment below to change the idv from TIME to something else such as TAD.
#Other xpose preferences could also be changed
#xpdb@Prefs@Xvardef$idv="TAD"
runsum(xpdb,dir=model.directory,
         modfile=paste(model.directory,model.filename,sep=""),
         listfile=paste(model.directory,sub(mod.suffix,".lst",model.filename),sep=""))
if (is.null(subset.variable)){
    print(basic.gof(xpdb))
    print(ranpar.hist(xpdb))
    print(ranpar.qq(xpdb))
    print(dv.preds.vs.idv(xpdb))
    print(dv.vs.idv(xpdb))
    print(ipred.vs.idv(xpdb))
    print(pred.vs.idv(xpdb))
    
}else{
    # change the subset variable from categorical to continuous or vice versa.
    # change.cat.cont(xpdb) <- c(subset.variable)
    print(basic.gof(xpdb,by=subset.variable,max.plots.per.page=1))
    print(ranpar.hist(xpdb,by=subset.variable,scales="free",max.plots.per.page=1))
    print(ranpar.qq(xpdb,by=subset.variable,max.plots.per.page=1))
    print(dv.preds.vs.idv(xpdb,by=subset.variable))
    print(dv.vs.idv(xpdb,by=subset.variable))
    print(ipred.vs.idv(xpdb,by=subset.variable))
    print(pred.vs.idv(xpdb,by=subset.variable))
}
  
if (rplots.level > 1){
  #individual plots of ten random IDs
  #find idcolumn
  idvar <- xvardef("id",xpdb)
  ten.random.ids<-sort(sample(unique(xpdb@Data[,idvar]),10,replace=FALSE))
  subset.string <- paste0(idvar,'==',paste(ten.random.ids,collapse=paste0(' | ',idvar,'==')))
  
  
  if (is.null(subset.variable)){
    print(ind.plots(xpdb,subset=subset.string))
  }else{
    for (flag in unique(xpdb@Data[,subset.variable])){
      print(ind.plots(xpdb,subset=paste0(subset.variable,'==',flag,' & (',subset.string,')')))
    }
  }
  
}

dev.off()


