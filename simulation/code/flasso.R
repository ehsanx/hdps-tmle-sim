# require(plyr) 
# require(dplyr) 
# require(WeightIt)
# require(data.table)
# require(sandwich)
# require(glmnet)
# require(cobalt)
# require(survey)
# require(SuperLearner)
# require(earth)
# require(xgboost)

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
      plasmodeData.i <- plyr::join(x = data.frame(idx=plasmode.data[,i],
                                                  EVENT=plasmode.data[,i+nSim]), 
                                   y = analytic.data, by="idx", type="left")
      plasmodeData.i$outcome <- plasmodeData.i$EVENT
      
      covarsTfull <- c(demovars, lab.list, proxy.list)
      Y.form <- as.formula(paste0(c("outcome~ exposure", covarsTfull), collapse = "+") )
      covar.mat <- model.matrix(Y.form, data = plasmodeData.i)[,-1]
      lasso.fit<-glmnet::cv.glmnet(y = plasmodeData.i$outcome, 
                                   x = covar.mat, 
                                   type.measure='mse',
                                   family="binomial",
                                   alpha = 1, nfolds = 5)
      coef.fit<-coef(lasso.fit,s='lambda.min',exact=TRUE)
      sel.variables<-row.names(coef.fit)[which(as.numeric(coef.fit)!=0)]
      proxy.list.sel <- proxy.list[proxy.list %in% sel.variables]
      
      saveRDS(proxy.list.sel, file = paste0("VarSel",method.under.consider,"/var",
                                            stringr::str_pad(i, 3, pad = "0"),
                                            "sel.Rds"))
      
      spw2 <- create.weights(dat = plasmodeData.i, 
                             demovars=demovars, 
                             lab.list=lab.list, 
                             proxy.list.sel=proxy.list.sel,
                             method.under.consider=method,
                             i=i)
      
      res <- extract.res(dat = plasmodeData.i, 
                         w2=spw2$w2, 
                         sw2=spw2$sw2, 
                         demovars=demovars, 
                         lab.list=lab.list, 
                         proxy.list.sel = proxy.list.sel, 
                         method.under.consider= method, 
                         i=i)
      
      saveRDS(res, file = paste0(method.under.consider,"/sim",
                                 stringr::str_pad(i, 3, pad = "0"),
                                 "res.Rds"))
      end_time <- Sys.time()
      run.timex <- end_time - start_time
      filename <- paste0(method.under.consider, "msg/iter", stringr::str_pad(i, 3, pad = "0"), ".txt")
      writeLines(as.character(run.timex), con = filename)
    })
  return(res0)
}
