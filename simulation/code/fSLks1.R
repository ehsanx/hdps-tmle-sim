# library(plyr) 
# library(WeightIt)
# library(data.table)
# library(sandwich)
# library(dplyr)
# library(cobalt)
# library(survey)
# library(tableone)

# plasmode.data = simdata
# analytic.data = full.data
# demovars = demovars
# lab.list = labvars.original
# proxy.list = emp.cov.names
# nSim = nSim
# method = method.under.consider
# i=1

loop.to.results.gen<- function(plasmode.data = simdata, 
                               analytic.data = full.data, 
                               demovars = demovars,
                               lab.list = labvars.original,
                               proxy.list = emp.cov.names,
                               nSim = nSim,
                               method = method.under.consider,
                               i=1){
  res0 <- NULL
  try({
    options(warn=-1)
    start_time <- Sys.time()
    methodx <- "tmle1"
    plasmodeData.i <- plyr::join(x = data.frame(idx=plasmode.data[,i],
                                                EVENT=plasmode.data[,i+nSim]), 
                                 y = analytic.data, by="idx", type="left")
    plasmodeData.i$outcome <- plasmodeData.i$EVENT
    demoform <- paste0(demovars, collapse = "+")
    labform0 <- paste0(lab.list, collapse = "+")
    rhsformula <- paste0(c(demoform, labform0), collapse = "+")
    proxyform <- paste0(proxy.list, collapse = "+")
    rhsformula2 <- paste0(c(demoform, labform0, proxyform), collapse = "+")
    
    ps.formula <- as.formula(paste0("exposure", "~", rhsformula))
    ps.formula2 <- as.formula(paste0("exposure", "~", rhsformula2))
    
    spw <- create.weights.sl(dat = plasmodeData.i, 
                             methodx = methodx,
                             method.under.consider=method.under.consider,
                             i=i)
    
    res <- extract.res(dat = plasmodeData.i, 
                       w2=spw$w2, 
                       sw2=spw$sw2, 
                       demovars=demovars, 
                       lab.list=lab.list, 
                       proxy.list.sel = proxy.list, 
                       method.under.consider= method, 
                       i=i)
    
    resx <- res
    
    saveRDS(resx, file = paste0(method.under.consider,"/sim",
                               stringr::str_pad(i, 3, pad = "0"),
                               "res.Rds"))
    end_time <- Sys.time()
    run.timex <- end_time - start_time
    filename <- paste0(method.under.consider, "msg/iter", stringr::str_pad(i, 3, pad = "0"), ".txt")
    writeLines(as.character(run.timex), con = filename)
  })
  return(res0)
}
