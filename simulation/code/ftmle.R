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
require(gam)
require(kernlab)
require(xgboost)
# require(rJava)
# require(bartMachine)
# require(BayesTree)
require(dbarts)

SL.nnet.m <- function(...){
  SL.nnet(..., size=4)
}

loop.to.results.gen<- function(plasmode.data = simdata, 
                               analytic.data = full.data, 
                               demovars = demovars,
                               lab.list = labvars.original,
                               proxy.list = emp.cov.names,
                               nSim = nSim,
                               method = method.under.consider,
                               i=1){
  res <- rep(NA,22)
  names(res) <- c("RD.psi.t", "RD.var.psi.t", "RD.CI1.t", "RD.CI2.t", "RD.pvalue.t", 
                  "OR.psi.t", "OR.CI1.t", "OR.CI2.t", "OR.pvalue.t", "OR.log.psi.t", "OR.var.log.psi.t", 
                  "RD.psi.tp", "RD.var.psi.tp", "RD.CI1.tp", "RD.CI2.tp", "RD.pvalue.tp", 
                  "OR.psi.tp", "OR.CI1.tp", "OR.CI2.tp", "OR.pvalue.tp", "OR.log.psi.tp", "OR.var.log.psi.tp")
  try({
    options(warn=-1)
    plasmodeData.i <- plyr::join(x = data.frame(idx=plasmode.data[,i],
                                                EVENT=plasmode.data[,i+nSim]), 
                                 y = analytic.data, by="idx", type="left")
    plasmodeData.i$outcome <- plasmodeData.i$EVENT
    demoform <- paste0(demovars, collapse = "+")
    labform0 <- paste0(lab.list, collapse = "+")
    rhsformula <- paste0(c(demoform, labform0), collapse = "+")
    proxyform <- paste0(proxy.list, collapse = "+")
    if (proxyform == "") {
      rhsformula2 <- paste0(c(demoform, labform0), collapse = "+")
    } else {
      rhsformula2 <- paste0(c(demoform, labform0, proxyform), collapse = "+")
    }
    
    ps.formula <- as.formula(paste0("exposure", "~", rhsformula))
    ps.formula2 <- as.formula(paste0("exposure", "~", rhsformula2))
    
    SL.library <- c("SL.glm", # 1.375387 secs
                    "SL.earth", # 36.90798 secs
                    "SL.glmnet", # 35.69747 secs
                    #"SL.nnet", # 31.56851 secs
                    #"SL.ksvm",  # 2.746763 mins
                    "SL.xgboost") # 1.969915 mins
    W.out <- weightit(ps.formula, 
                      data = plasmodeData.i, 
                      estimand = "ATE",
                      method = "super",
                      SL.library = SL.library)
    ObsData.noYA <- plasmodeData.i[,c(demovars, lab.list)]
    
    tmle.fit <- tmle::tmle(Y = plasmodeData.i$EVENT, 
                           A = plasmodeData.i$exposure, 
                           W = ObsData.noYA, 
                           family = "binomial",
                           # V = 3,
                           V.Q = 3,
                           V.g = 3,
                           Q.SL.library = SL.library,
                           g1W = W.out$ps)
    estRD.tmle <- tmle.fit$estimates$ATE
    n.list <- names(unlist(estRD.tmle))
    n.listRD <- paste0("RD.",n.list)
    estRD.tmle <- as.vector(unlist(estRD.tmle))
    estOR.tmle <- tmle.fit$estimates$OR
    n.list2 <- names(unlist(estOR.tmle))
    n.listOR <- paste0("OR.",n.list2)
    estOR.tmle <- as.vector(unlist(estOR.tmle))
    est <- c(estRD.tmle, estOR.tmle)
    names(est) <- c(n.listRD, n.listOR)
    
    W.out2 <- weightit(ps.formula2, 
                       data = plasmodeData.i, 
                       estimand = "ATE",
                       method = "super",
                       SL.library = SL.library)
    ObsData.noYA2 <- plasmodeData.i[,c(demovars, lab.list, proxy.list)]
    
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
    estRD.tmle2 <- as.vector(unlist(tmle.fit2$estimates$ATE))
    estOR.tmle2 <- as.vector(unlist(tmle.fit2$estimates$OR))
    est2 <- c(estRD.tmle2, estOR.tmle2)
    names(est2) <- names(est)
    
    ps <- W.out$ps
    ps2 <- W.out2$ps
    w <- W.out$weights
    w2 <- W.out2$weights
    spw <- as.data.frame(cbind(ps,ps2,w,w2))
    saveRDS(spw, file = paste0("../weighttmle/pw",
                               stringr::str_pad(i, 3, pad = "0"),
                               "sl.RDS"))
    
    res <- c(est, est2)
    names(res) <- as.vector(outer(names(est),c("t","tp"), paste, sep="."))
    
    saveRDS(res, file = paste0("../",method.under.consider,"/sim",
                               stringr::str_pad(i, 3, pad = "0"),
                               "res.Rds"))
  })
  return(res)
}
