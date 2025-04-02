# require(plyr) 
# require(WeightIt)
# require(cobalt)
# require(survey)
# require(data.table)
# require(sandwich)
require(brglm2)

loop.to.results.gen<- function(plasmode.data = simdata, 
                               analytic.data = full.data, 
                               demovars = demovars,
                               lab.list = labvars.original,
                               proxy.list = emp.cov.names,
                               nSim = nSim,
                               method = method.under.consider,
                               i=1){
  res0 <- NULL
  #try({
    options(warn=-1)
    start_time <- Sys.time()
    options(warn=-1)
    plasmodeData.i <- plyr::join(x = data.frame(idx=plasmode.data[,i],
                                                EVENT=plasmode.data[,i+nSim]), 
                                 y = analytic.data, by="idx", type="left")
    plasmodeData.i$outcome <- plasmodeData.i$EVENT
    
    spw0 <- create.weights(dat = plasmodeData.i, 
                           demovars=demovars, 
                           lab.list=lab.list, 
                           proxy.list.sel=NULL,
                           method.under.consider=method,
                           i=i)
    
    spw2 <- create.weights(dat = plasmodeData.i, 
                           demovars=demovars, 
                           lab.list=lab.list, 
                           proxy.list.sel=proxy.list,
                           method.under.consider=method,
                           i=i)
    
    res0 <- extract.res.f2(dat = plasmodeData.i, 
                       w2=spw0$w2, 
                       sw2=spw0$sw2, 
                       demovars=demovars, 
                       lab.list=lab.list, 
                       proxy.list.sel = NULL, 
                       method.under.consider= method, 
                       i=i)
    
    res <- extract.res.f2(dat = plasmodeData.i, 
                       w2=spw2$w2, 
                       sw2=spw2$sw2, 
                       demovars=demovars, 
                       lab.list=lab.list, 
                       proxy.list.sel = proxy.list, 
                       method.under.consider= method, 
                       i=i)
    res.both <- c(res0, res)
    saveRDS(res.both, file = paste0(method.under.consider,"/sim",
                               stringr::str_pad(i, 3, pad = "0"),
                               "res.Rds"))
    end_time <- Sys.time()
    run.timex <- end_time - start_time
    filename <- paste0(method.under.consider, "msg/iter", stringr::str_pad(i, 3, pad = "0"), ".txt")
    writeLines(as.character(run.timex), con = filename)
  #})
  return(res0)
}
