require(autoCovariateSelection)
require(dplyr) 
require(WeightIt)
require(cobalt)
require(survey)
require(data.table)
require(sandwich)
require(glmnet)
require(plyr) 
require(tableone)

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
    methodx <- "hdpslassotmle"
    plasmodeData.i <- plyr::join(x = data.frame(idx=plasmode.data[,i],
                                                EVENT=plasmode.data[,i+nSim]), 
                                 y = analytic.data, by="idx", type="left")
    plasmodeData.i$outcome <- plasmodeData.i$EVENT
    plasmodeData.i$id <- 1:nrow(plasmodeData.i)
    proxy.data <- plasmodeData.i[,proxy.list]
    proxy.data$id <- 1:nrow(proxy.data)
    
    # Compute the sum for each column
    col_sums <- apply(proxy.data, 2, sum)
    
    # Find the names of columns with all 0 or all 1
    all_zero_or_one_colnames <- colnames(proxy.data)[col_sums == 0 | col_sums == nrow(proxy.data)]
    
    # Remove the columns with all 0 or all 1
    proxy.data.filtered <- proxy.data[, !(colnames(proxy.data) %in% all_zero_or_one_colnames)]
    
    out3 <- get_prioritised_covariates(df = proxy.data.filtered,
                                       patientIdVarname = "id", 
                                       exposureVector = plasmodeData.i$exposure,
                                       outcomeVector = plasmodeData.i$outcome,
                                       patientIdVector = plasmodeData.i$id, 
                                       k = 100)
    hdps.dim <- out3$autoselected_covariate_df
    hdps.data <- merge(plasmodeData.i[,c("id",
                                         "outcome", 
                                         "exposure", 
                                         demovars, 
                                         lab.list)], 
                       hdps.dim, by = "id")
    proxy.list.sel0 <- names(out3$autoselected_covariate_df[,-1])
    
    covarsTfull <- c(demovars, lab.list, proxy.list.sel0)
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

    fake_sim <- "fake_sim"

    saveRDS(fake_sim, file = paste0(method.under.consider,"/sim",
                                stringr::str_pad(i, 3, pad = "0"),
                                "res.Rds"))
    
    spw <- create.weights.sl(dat = plasmodeData.i, 
                             methodx = methodx,
                             method.under.consider=method.under.consider,
                             i=i)
    
    res <- extract.res(dat = plasmodeData.i, 
                       w2=spw$w2, 
                       sw2=spw$sw2, 
                       demovars=demovars, 
                       lab.list=lab.list, 
                       proxy.list.sel = proxy.list.sel, 
                       method.under.consider= method, 
                       i=i)
    res3 <- extract.res(dat = plasmodeData.i, 
                        w2=spw$w23, 
                        sw2=spw$sw23, 
                        demovars=demovars, 
                        lab.list=lab.list, 
                        proxy.list.sel = proxy.list.sel, 
                        method.under.consider= method, 
                        sln = 3,
                        i=i)
    
    resx <- c(res,res3)
    
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
