require(autoCovariateSelection)
require(dplyr) 
require(WeightIt)
require(cobalt)
require(survey)
require(data.table)
require(sandwich)
require(SuperLearner)
require(tmle)
require(earth)
require(glmnet)
require(nnet)
require(kernlab)
require(xgboost)

loop.to.results.gen<- function(plasmode.data = simdata, 
                               analytic.data = full.data, 
                               demovars = demovars,
                               lab.list = labvars.original,
                               proxy.list = emp.cov.names,
                               nSim = nSim,
                               method = method.under.consider,
                               i=1){
  res <- rep(NA,11)
  names(res) <- c("RD.psi", "RD.var.psi", "RD.CI1", "RD.CI2", "RD.pvalue", "OR.psi", 
                  "OR.CI1", "OR.CI2", "OR.pvalue", "OR.log.psi", "OR.var.log.psi")
  
  try({
    options(warn=-1)
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
    lasso.fit<-cv.glmnet(y = plasmodeData.i$outcome, 
                         x = covar.mat, 
                         type.measure='mse',
                         family="binomial",
                         alpha = 1, nfolds = 5)
    coef.fit<-coef(lasso.fit,s='lambda.min',exact=TRUE)
    sel.variables<-row.names(coef.fit)[which(coef.fit!=0)]
    proxy.list.sel <- proxy.list.sel0[proxy.list.sel0 %in% sel.variables]
    
    demoform <- paste0(demovars, collapse = "+")
    labform0 <- paste0(lab.list, collapse = "+")
    # rhsformula <- paste0(c(demoform, labform0), collapse = "+")
    proxyform <- paste0(proxy.list.sel, collapse = "+")
    rhsformula2 <- paste0(c(demoform, labform0, proxyform), collapse = "+")
    
    # ps.formula <- as.formula(paste0("exposure", "~", rhsformula))
    ps.formula2 <- as.formula(paste0("exposure", "~", rhsformula2))
    # W.out <- weightit(ps.formula, 
    #                   #stabilize = TRUE,
    #                   data = plasmodeData.i, 
    #                   estimand = "ATE",
    #                   method = "ps")
    
    SL.library <- c("SL.glm", # 1.375387 secs
                    "SL.earth", # 36.90798 secs
                    "SL.glmnet", # 35.69747 secs
                    #"SL.nnet", # 31.56851 secs
                    #"SL.ksvm",  # 2.746763 mins
                    "SL.xgboost") # 1.969915 mins
    
    W.out2 <- weightit(ps.formula2, 
                       data = plasmodeData.i, 
                       estimand = "ATE",
                       method = "super",
                       SL.library = SL.library)
    ObsData.noYA2 <- plasmodeData.i[,c(demovars, lab.list, proxy.list.sel)]
    
    ObsData.noYA2$idx <- NULL
    tmle.fit2 <- tmle::tmle(Y = plasmodeData.i$EVENT, 
                            A = plasmodeData.i$exposure, 
                            W = ObsData.noYA2, 
                            family = "binomial",
                            # V = 3,
                            V.Q = 3,
                            V.g = 3,
                            Q.SL.library = SL.library,
                            g1W = W.out2$ps)
    estRD.tmle2 <- unlist(tmle.fit2$estimates$ATE)
    n.list <- names(unlist(estRD.tmle2))
    n.listRD <- paste0("RD.",n.list)
    names(estRD.tmle2) <- n.listRD
    estOR.tmle2 <- unlist(tmle.fit2$estimates$OR)
    n.list2 <- names(unlist(estOR.tmle2))
    n.listOR <- paste0("OR.",n.list2)
    names(estOR.tmle2) <- n.listOR
    res <- c(estRD.tmle2, estOR.tmle2)
    
    ps2 <- W.out2$ps
    w2 <- W.out2$weights
    spw <- as.data.frame(cbind(ps2,w2))
    saveRDS(spw, file = paste0("weight",method.under.consider,"/pw",
                               stringr::str_pad(i, 3, pad = "0"),
                               "sl.RDS"))
    
    saveRDS(res, file = paste0(method.under.consider,"/sim",
                               stringr::str_pad(i, 3, pad = "0"),
                               "res.Rds"))
    
    saveRDS(proxy.list.sel, file = paste0("VarSel",method.under.consider,"/var",
                                          stringr::str_pad(i, 3, pad = "0"),
                                          "sel.Rds"))
    
  })
  return(res)
}
