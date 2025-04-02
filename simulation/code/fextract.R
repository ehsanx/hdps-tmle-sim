get_smd_stats <- function(w, 
                          dat, 
                          demovars,
                          lab.list,
                          proxy.list = NULL) {
  cols_to_select <- c(demovars, lab.list)
  if (!is.null(proxy.list)) {
    cols_to_select <- c(cols_to_select, proxy.list)
  }
  cols_to_select <- c(cols_to_select, "exposure")
  
  dat_sub <- subset(dat, select = cols_to_select)
  dat_sub$w <- w
  
  dclus1<- survey::svydesign(id=~1,weights=~w,data=dat_sub)
  tab1 <- tableone::svyCreateTableOne(vars = cols_to_select,
                                      strata = "exposure", 
                                      data = dclus1)
  
  smds <- tableone::ExtractSmd(tab1)
  smds2 <- smds[!rownames(smds) %in% c("exposure", "w"), ]
  smds3 <- smds[!rownames(smds) %in% c("exposure", "w", proxy.list), ]
  
  prop_less_than_02 <- sum(smds2 < 0.2) / length(smds2)
  prop_less_than_01 <- sum(smds2 < 0.1) / length(smds2)
  imbalanced <- names(smds3[smds3>0.1])
  imbalancedb <- names(smds2[smds2>0.1])
  n.imbalanced <- length(imbalanced)
  smd_summary <- summary(smds2)
  w_summary <- summary(w)
  prop_less_than_02base <- sum(smds3 < 0.2) / length(smds3)
  prop_less_than_01base <- sum(smds3 < 0.1) / length(smds3)
  smd_summarybase <- summary(smds3)
  
  return(list(prop_less_than_02 = prop_less_than_02,
              prop_less_than_01 = prop_less_than_01,
              imbalanced = imbalanced,
              imbalancedb = imbalancedb,
              n.imbalanced = n.imbalanced,
              smd_summary = smd_summary,
              w_summary = w_summary,
              prop_less_than_02base= prop_less_than_02base,
              prop_less_than_01base = prop_less_than_01base,
              smd_summarybase = smd_summarybase))
}

extract.res <- function(dat = plasmodeData.i, 
                        w2, 
                        sw2, 
                        demovars, 
                        lab.list, 
                        proxy.list.sel = NULL, 
                        method.under.consider, 
                        sln = NULL,
                        i) {
  
  smd_list <- get_smd_stats(w=w2, dat=dat, demovars=demovars, lab.list=lab.list, proxy.list = proxy.list.sel)
  saveRDS(smd_list, file = paste0(method.under.consider,"SMD/W",sln,"smd",
                                  stringr::str_pad(i, 3, pad = "0"),
                                  "res.Rds"))
  # Adj 1: only imbalanced covariates and proxies
  imbalanced <- unique(smd_list$imbalanced)
  #imbalanced <- imbalanced[!(imbalanced %in% c(demovars, lab.list))]
  imbalancedform <- if (length(imbalanced) == 0) NULL else if (length(imbalanced) > 1) paste0(imbalanced, collapse = '+') else imbalanced
  demoform <- paste0(demovars, collapse = "+")
  labform0 <- paste0(lab.list, collapse = "+")
  rhsformula2x <- if (length(imbalanced) == 0) {
    paste0(c(demoform, labform0), collapse = "+")
  } else {
    paste0(c(demoform, labform0, imbalancedform), collapse = "+")
  }
  if (length(imbalanced) == 0) out.formula2x <- as.formula(paste0("EVENT", "~", "exposure"))
  if (length(imbalanced) > 0) out.formula2x <- as.formula(paste0("EVENT", "~", "exposure + ", imbalancedform))
  fit2 <- glm(out.formula2x,
              data = dat,weights = w2,
              family= gaussian(link = "identity"))
  estRDx2 <- summary(fit2)$coef["exposure",c("Estimate", "Std. Error", "Pr(>|t|)")]
  estRDx2[2] <- sqrt(sandwich::sandwich(fit2)[2,2])
  
  # Adj 2: no adjustment
  out.formula20 <- as.formula(paste0("EVENT", "~", "exposure"))
  fit20 <- glm(out.formula20,
               data = dat,weights = w2,
               family= gaussian(link = "identity"))
  estRDx20 <- summary(fit20)$coef["exposure",c("Estimate", "Std. Error", "Pr(>|t|)")]
  estRDx20[2] <- sqrt(sandwich::sandwich(fit20)[2,2])
  
  # Adj 3: all covariates and proxies
  proxyform <- paste0(proxy.list.sel, collapse = "+") 
  #rhsformula2 <- paste0(c(demoform, labform0, proxyform), collapse = "+")
  rhsformula2 <- paste0(c(demoform, labform0), collapse = "+")
  if (!is.null(proxy.list.sel)) {
    rhsformula2 <- paste0(rhsformula2, "+", proxyform)
  }
  
  out.formula20a <- as.formula(paste0("EVENT", "~", "exposure + ", rhsformula2))
  fit20a <- glm(out.formula20a,
                data = dat,weights = w2,
                family= gaussian(link = "identity"))
  estRDx20a <- summary(fit20a)$coef["exposure",c("Estimate", "Std. Error", "Pr(>|t|)")]
  estRDx20a[2] <- sqrt(sandwich::sandwich(fit20a)[2,2])
  
  # Adj 4: only imbalanced covariates
  imbalanced_subset <- imbalanced[imbalanced %in% union(demovars, lab.list)]
  imbalancedform <- if (length(imbalanced_subset) == 0) NULL else if (length(imbalanced_subset) > 1) paste0(imbalanced_subset, collapse = '+') else imbalanced_subset
  demoform <- paste0(demovars, collapse = "+")
  labform0 <- paste0(lab.list, collapse = "+")
  rhsformula2mx <- if (length(imbalanced_subset) == 0) {
    paste0(c(demoform, labform0), collapse = "+")
  } else {
    paste0(c(demoform, labform0, imbalancedform), collapse = "+")
  }
  if (length(imbalanced_subset) == 0) out.formula2mx <- as.formula(paste0("EVENT", "~", "exposure"))
  if (length(imbalanced_subset) > 0) out.formula2mx <- as.formula(paste0("EVENT", "~", "exposure + ", imbalancedform))
  fit2m <- glm(out.formula2mx,
               data = dat,weights = w2,
               family= gaussian(link = "identity"))
  estRDx2m <- summary(fit2m)$coef["exposure",c("Estimate", "Std. Error", "Pr(>|t|)")]
  estRDx2m[2] <- sqrt(sandwich::sandwich(fit2m)[2,2])
  
  # Adj 5: all baseline
  out.formula2b <- as.formula(paste0("EVENT", "~", "exposure + ", 
                                     paste0(c(demoform, labform0), collapse = "+")))
  fit2b <- glm(out.formula2b,
               data = dat,weights = w2,
               family= gaussian(link = "identity"))
  estRDx2b <- summary(fit2b)$coef["exposure",c("Estimate", "Std. Error", "Pr(>|t|)")]
  estRDx2b[2] <- sqrt(sandwich::sandwich(fit2b)[2,2])
  
  smd_list.sw <- get_smd_stats(w=sw2, dat=dat, demovars=demovars, lab.list=lab.list, proxy.list = proxy.list.sel)
  saveRDS(smd_list.sw, file = paste0(method.under.consider,"SMD/SW",sln,"smd",
                                     stringr::str_pad(i, 3, pad = "0"),
                                     "res.Rds"))
  # Adj 1: only imbalanced covariates and proxies
  imbalanced.sw <- unique(smd_list.sw$imbalanced)
  imbalancedform <- if (length(imbalanced.sw) == 0) NULL else if (length(imbalanced.sw) > 1) paste0(imbalanced.sw, collapse = '+') else imbalanced.sw
  rhsformula2s <- if (length(imbalanced.sw) == 0) {
    paste0(c(demoform, labform0), collapse = "+")
  } else {
    paste0(c(demoform, labform0, imbalancedform), collapse = "+")
  }
  
  if (length(imbalanced.sw) == 0) out.formula2s <- as.formula(paste0("EVENT", "~", "exposure"))
  if (length(imbalanced.sw) > 0) out.formula2s <- as.formula(paste0("EVENT", "~", "exposure + ", imbalancedform))
  fit2s <- glm(out.formula2s,
               data = dat,weights = sw2,
               family= gaussian(link = "identity"))
  estRDx2s <- summary(fit2s)$coef["exposure",c("Estimate", "Std. Error", "Pr(>|t|)")]
  estRDx2s[2] <- sqrt(sandwich::sandwich(fit2s)[2,2])
  
  # Adj 2: no adjustment
  out.formula2s0 <- as.formula(paste0("EVENT", "~", "exposure"))
  fit2s0 <- glm(out.formula2s0,
                data = dat,weights = sw2,
                family= gaussian(link = "identity"))
  estRDx2s0 <- summary(fit2s0)$coef["exposure",c("Estimate", "Std. Error", "Pr(>|t|)")]
  estRDx2s0[2] <- sqrt(sandwich::sandwich(fit2s0)[2,2])
  
  # Adj 3: all covariates and proxies
  out.formula2a <- as.formula(paste0("EVENT", "~", "exposure + ", rhsformula2))
  fit2a <- glm(out.formula2a,
               data = dat,weights = sw2,
               family= gaussian(link = "identity"))
  estRDx2a <- summary(fit2a)$coef["exposure",c("Estimate", "Std. Error", "Pr(>|t|)")]
  estRDx2a[2] <- sqrt(sandwich::sandwich(fit2a)[2,2])
  
  # Adj 4: only imbalanced covariates
  imbalanced.sw_subset <- imbalanced.sw[imbalanced.sw %in% union(demovars, lab.list)]
  imbalancedform <- if (length(imbalanced.sw_subset) == 0
                        ) NULL else if (length(imbalanced.sw_subset) > 1
                                        ) paste0(imbalanced.sw_subset, collapse = '+'
                                                 ) else imbalanced.sw_subset
  demoform <- paste0(demovars, collapse = "+")
  labform0 <- paste0(lab.list, collapse = "+")
  rhsformula2msx <- if (length(imbalanced.sw_subset) == 0) {
    paste0(c(demoform, labform0), collapse = "+")
  } else {
    paste0(c(demoform, labform0, imbalancedform), collapse = "+")
  }
  
  if (length(imbalanced.sw_subset) == 0) out.formula2msx <- as.formula(paste0("EVENT", "~", "exposure"))
  if (length(imbalanced.sw_subset) > 0) out.formula2msx <- as.formula(paste0("EVENT", "~", "exposure + ", imbalancedform))
  fit2ms <- glm(out.formula2msx,
               data = dat,weights = sw2,
               family= gaussian(link = "identity"))
  estRDx2ms <- summary(fit2ms)$coef["exposure",c("Estimate", "Std. Error", "Pr(>|t|)")]
  estRDx2ms[2] <- sqrt(sandwich::sandwich(fit2ms)[2,2])
  
  # Adj 5: all baseline
  fit2bs <- glm(out.formula2b,
               data = dat,weights = sw2,
               family= gaussian(link = "identity"))
  estRDx2bs <- summary(fit2bs)$coef["exposure",c("Estimate", "Std. Error", "Pr(>|t|)")]
  estRDx2bs[2] <- sqrt(sandwich::sandwich(fit2bs)[2,2])
  
  res <- c(estRDx2, estRDx20, estRDx20a, estRDx2m, estRDx2b,
           estRDx2s, estRDx2s0, estRDx2a, estRDx2ms, estRDx2bs)
  resnx <- c("wi", "w0", "wa", "wm", "wb",
             "swi", "sw0", "swa", "swm", "swb")
  if (is.null(proxy.list.sel)) {
    name.change <- "base"
  } else {
    name.change <- NULL
  }
  names(res) <- as.vector(outer(c("b","se","p"),
                                paste0(resnx, sln,name.change), 
                                paste, sep="-"))
  return(res)
}

extract.res.f2 <- function(dat = plasmodeData.i, 
                           w2, 
                           sw2, 
                           demovars, 
                           lab.list, 
                           proxy.list.sel = NULL, 
                           method.under.consider, 
                           sln = NULL,
                           i) {
  
  smd_list <- get_smd_stats(w=w2, dat=dat, demovars=demovars, lab.list=lab.list, proxy.list = proxy.list.sel)
  saveRDS(smd_list, file = paste0(method.under.consider,"SMD/W",sln,"smd",
                                  stringr::str_pad(i, 3, pad = "0"),
                                  "res.Rds"))
  # Adj 1: only imbalanced covariates and proxies
  imbalanced <- unique(smd_list$imbalanced)
  #imbalanced <- imbalanced[!(imbalanced %in% c(demovars, lab.list))]
  imbalancedform <- if (length(imbalanced) == 0) NULL else if (length(imbalanced) > 1) paste0(imbalanced, collapse = '+') else imbalanced
  demoform <- paste0(demovars, collapse = "+")
  labform0 <- paste0(lab.list, collapse = "+")
  rhsformula2x <- if (length(imbalanced) == 0) {
    paste0(c(demoform, labform0), collapse = "+")
  } else {
    paste0(c(demoform, labform0, imbalancedform), collapse = "+")
  }
  if (length(imbalanced) == 0) out.formula2x <- as.formula(paste0("EVENT", "~", "exposure"))
  if (length(imbalanced) > 0) out.formula2x <- as.formula(paste0("EVENT", "~", "exposure + ", imbalancedform))
  fit2 <- glm(out.formula2x,
              data = dat,weights = w2,
              family= gaussian(link = "identity"))
  estRDx2 <- extract.bse(fit2, dataf = dat, fx = out.formula2x)
  
  # Adj 2: no adjustment
  out.formula20 <- as.formula(paste0("EVENT", "~", "exposure"))
  fit20 <- glm(out.formula20,
               data = dat,weights = w2,
               family= gaussian(link = "identity"))
  estRDx20 <- extract.bse(fit20, dataf = dat, fx = out.formula20)
  
  # Adj 3: all covariates and proxies
  proxyform <- paste0(proxy.list.sel, collapse = "+") 
  #rhsformula2 <- paste0(c(demoform, labform0, proxyform), collapse = "+")
  rhsformula2 <- paste0(c(demoform, labform0), collapse = "+")
  if (!is.null(proxy.list.sel)) {
    rhsformula2 <- paste0(rhsformula2, "+", proxyform)
  }
  
  out.formula20a <- as.formula(paste0("EVENT", "~", "exposure + ", rhsformula2))
  fit20a <- glm(out.formula20a,
                data = dat,weights = w2,
                family= gaussian(link = "identity"))
  estRDx20a <- extract.bse(fit20a, dataf = dat, fx = out.formula20a)
  
  # Adj 4: only imbalanced covariates
  imbalanced_subset <- imbalanced[imbalanced %in% union(demovars, lab.list)]
  imbalancedform <- if (length(imbalanced_subset) == 0) NULL else if (length(imbalanced_subset) > 1) paste0(imbalanced_subset, collapse = '+') else imbalanced_subset
  demoform <- paste0(demovars, collapse = "+")
  labform0 <- paste0(lab.list, collapse = "+")
  rhsformula2mx <- if (length(imbalanced_subset) == 0) {
    paste0(c(demoform, labform0), collapse = "+")
  } else {
    paste0(c(demoform, labform0, imbalancedform), collapse = "+")
  }
  
  if (length(imbalanced_subset) == 0) out.formula2mx <- as.formula(paste0("EVENT", "~", "exposure"))
  if (length(imbalanced_subset) > 0) out.formula2mx <- as.formula(paste0("EVENT", "~", "exposure + ", imbalancedform))
  fit2m <- glm(out.formula2mx,
               data = dat,weights = w2,
               family= gaussian(link = "identity"))
  estRDx2m <- extract.bse(fit2m, dataf = dat, fx = out.formula2mx)
  
  # Adj 5: all baseline
  out.formula2b <- as.formula(paste0("EVENT", "~", "exposure + ", 
                                     paste0(c(demoform, labform0), collapse = "+")))
  fit2b <- glm(out.formula2b,
               data = dat,weights = w2,
               family= gaussian(link = "identity"))
  estRDx2b <- extract.bse(fit2b, dataf = dat, fx = out.formula2b)
  
  smd_list.sw <- get_smd_stats(w=sw2, dat=dat, demovars=demovars, lab.list=lab.list, proxy.list = proxy.list.sel)
  saveRDS(smd_list.sw, file = paste0(method.under.consider,"SMD/SW",sln,"smd",
                                     stringr::str_pad(i, 3, pad = "0"),
                                     "res.Rds"))
  # Adj 1: only imbalanced covariates and proxies
  imbalanced.sw <- unique(smd_list.sw$imbalanced)
  imbalancedform <- if (length(imbalanced.sw) == 0) NULL else if (length(imbalanced.sw) > 1) paste0(imbalanced.sw, collapse = '+') else imbalanced.sw
  rhsformula2s <- if (length(imbalanced.sw) == 0) {
    paste0(c(demoform, labform0), collapse = "+")
  } else {
    paste0(c(demoform, labform0, imbalancedform), collapse = "+")
  }
  
  if (length(imbalanced.sw) == 0) out.formula2s <- as.formula(paste0("EVENT", "~", "exposure"))
  if (length(imbalanced.sw) > 0) out.formula2s <- as.formula(paste0("EVENT", "~", "exposure + ", imbalancedform))
  fit2s <- glm(out.formula2s,
               data = dat,weights = sw2,
               family= gaussian(link = "identity"))
  estRDx2s <- extract.bse(fit2s, dataf = dat, fx = out.formula2s)
  
  # Adj 2: no adjustment
  out.formula2s0 <- as.formula(paste0("EVENT", "~", "exposure"))
  fit2s0 <- glm(out.formula2s0,
                data = dat,weights = sw2,
                family= gaussian(link = "identity"))
  estRDx2s0 <- extract.bse(fit2s0, dataf = dat, fx = out.formula2s0)
  
  # Adj 3: all covariates and proxies
  out.formula2a <- as.formula(paste0("EVENT", "~", "exposure + ", rhsformula2))
  fit2a <- glm(out.formula2a,
               data = dat,weights = sw2,
               family= gaussian(link = "identity"))
  estRDx2a <- extract.bse(fit2a, dataf = dat, fx = out.formula2a)
  
  # Adj 4: only imbalanced covariates
  imbalanced.sw_subset <- imbalanced.sw[imbalanced.sw %in% union(demovars, lab.list)]
  imbalancedform <- if (length(imbalanced.sw_subset) == 0
  ) NULL else if (length(imbalanced.sw_subset) > 1
  ) paste0(imbalanced.sw_subset, collapse = '+'
  ) else imbalanced.sw_subset
  demoform <- paste0(demovars, collapse = "+")
  labform0 <- paste0(lab.list, collapse = "+")
  rhsformula2msx <- if (length(imbalanced.sw_subset) == 0) {
    paste0(c(demoform, labform0), collapse = "+")
  } else {
    paste0(c(demoform, labform0, imbalancedform), collapse = "+")
  }
  if (length(imbalanced.sw_subset) == 0) out.formula2msx <- as.formula(paste0("EVENT", "~", "exposure"))
  if (length(imbalanced.sw_subset) > 0) out.formula2msx <- as.formula(paste0("EVENT", "~", "exposure + ", imbalancedform))
  fit2ms <- glm(out.formula2msx,
                data = dat,weights = sw2,
                family= gaussian(link = "identity"))
  estRDx2ms <- extract.bse(fit2ms, dataf = dat, fx = out.formula2msx)
  
  # Adj 5: all baseline
  fit2bs <- glm(out.formula2b,
                data = dat,weights = sw2,
                family= gaussian(link = "identity"))
  estRDx2bs <- extract.bse(fit2bs, dataf = dat, fx = out.formula2b)
  
  res <- c(estRDx2, estRDx20, estRDx20a, estRDx2m, estRDx2b,
           estRDx2s, estRDx2s0, estRDx2a, estRDx2ms, estRDx2bs)
  resnx <- c("wi", "w0", "wa", "wm", "wb",
             "swi", "sw0", "swa", "swm", "swb")
  if (is.null(proxy.list.sel)) {
    name.change <- "base"
  } else {
    name.change <- NULL
  }
  names(res) <- as.vector(outer(c("b","se","p"),
                                paste0(resnx, sln,name.change), 
                                paste, sep="-"))
  return(res)
}

extract.bse <- function(fit, fx, dataf) {
  fit <- update(fit, formula. = fx, method = "brglmFit", type = "correction", data=dataf)
  fit_summary <- summary(fit)
  est.summary <- fit_summary$coef["exposure", c("Estimate", "Std. Error", "Pr(>|z|)")]
  est.summary[2] <- sqrt(sandwich::sandwich(fit)[2, 2])
  names(est.summary) <- c("Estimate", "Std. Error", "Pr(>|t|)")
  return(est.summary)
}



create.weights <- function(dat = plasmodeData.i, 
                           demovars, 
                           lab.list, 
                           proxy.list.sel = NULL,
                           method.under.consider,
                           i) {
  demoform <- paste0(demovars, collapse = "+")
  labform0 <- paste0(lab.list, collapse = "+")
  proxyform <- paste0(proxy.list.sel, collapse = "+")
  rhsformula2 <- paste0(c(demoform, labform0), collapse = "+")
  if (!is.null(proxy.list.sel)) {
    rhsformula2 <- paste0(rhsformula2, "+", proxyform)
  }
  
  ps.formula2 <- as.formula(paste0("exposure", "~", rhsformula2))
  W.out2 <- WeightIt::weightit(ps.formula2,
                               data = dat,
                               estimand = "ATE",
                               method = "ps")
  ps2 <- W.out2$ps
  w2 <- W.out2$weights
  wt <- WeightIt::trim(W.out2$weights, at = 0.99)
  sw2 <- WeightIt::get_w_from_ps(W.out2$ps,
                                 treat = dat$exposure,
                                 estimand = "ATE",
                                 stabilize = TRUE)
  swt <- WeightIt::trim(sw2, at = 0.99)
  spw2 <- as.data.frame(cbind(ps2, w2, wt, sw2, swt))
  
  if (is.null(proxy.list.sel)) {
    name.change <- "base"
  } else {
    name.change <- NULL
  }
  
  saveRDS(spw2, file = paste0("weight",method.under.consider,"/pw",
                              stringr::str_pad(i, 3, pad = "0"),
                              name.change, ".Rds"))
  return(spw2)
}


create.weights.sl <- function(dat = plasmodeData.i, 
                              methodx = methodx,
                              method.under.consider,
                              i) {
  
  W.out <- readRDS(file = paste0("weight",methodx,"/pw",
                                 stringr::str_pad(i, 3, pad = "0"),
                                 "sl.RDS"))
  ps2 <- W.out$ps2
  w2 <- W.out$w2
  wt <- WeightIt::trim(w2, at = 0.99)
  sw2 <- WeightIt::get_w_from_ps(ps2,
                                 treat = dat$exposure,
                                 estimand = "ATE",
                                 stabilize = TRUE)
  swt <- WeightIt::trim(sw2, at = 0.99)
  spw2 <- as.data.frame(cbind(ps2, w2, wt, sw2, swt))
  
  W.out3 <- readRDS(file = paste0("weight",methodx,"S/pw",
                                  stringr::str_pad(i, 3, pad = "0"),
                                  "sl.RDS"))
  ps23 <- W.out3$ps2
  w23 <- W.out3$w2
  wt3 <- WeightIt::trim(w23, at = 0.99)
  sw23 <- WeightIt::get_w_from_ps(ps23,
                                 treat = dat$exposure,
                                 estimand = "ATE",
                                 stabilize = TRUE)
  swt3 <- WeightIt::trim(sw23, at = 0.99)
  spw23 <- as.data.frame(cbind(ps23, w23, wt3, sw23, swt3))
  
  
  spw <- as.data.frame(cbind(spw2, spw23))
  
  saveRDS(spw, file = paste0("weight",method.under.consider,"/pw",
                              stringr::str_pad(i, 3, pad = "0"),
                              ".Rds"))
  return(spw)
}

create.weights.sl0 <- function(dat = plasmodeData.i, 
                              methodx = methodx,
                              method.under.consider,
                              i) {
  
  W.out <- readRDS(file = paste0("weight",methodx,"/pw",
                                 stringr::str_pad(i, 3, pad = "0"),
                                 "sl.RDS"))
  ps2 <- W.out$ps
  w2 <- W.out$w
  wt <- WeightIt::trim(w2, at = 0.99)
  sw2 <- WeightIt::get_w_from_ps(ps2,
                                 treat = dat$exposure,
                                 estimand = "ATE",
                                 stabilize = TRUE)
  swt <- WeightIt::trim(sw2, at = 0.99)
  spw2 <- as.data.frame(cbind(ps2, w2, wt, sw2, swt))
  
  W.out3 <- readRDS(file = paste0("weight",methodx,"S/pw",
                                  stringr::str_pad(i, 3, pad = "0"),
                                  "sl.RDS"))
  ps23 <- W.out3$ps
  w23 <- W.out3$w
  wt3 <- WeightIt::trim(w23, at = 0.99)
  sw23 <- WeightIt::get_w_from_ps(ps23,
                                  treat = dat$exposure,
                                  estimand = "ATE",
                                  stabilize = TRUE)
  swt3 <- WeightIt::trim(sw23, at = 0.99)
  spw23 <- as.data.frame(cbind(ps23, w23, wt3, sw23, swt3))
  
  
  spw <- as.data.frame(cbind(spw2, spw23))
  
  saveRDS(spw, file = paste0("weight",method.under.consider,"/pwu",
                             stringr::str_pad(i, 3, pad = "0"),
                             ".Rds"))
  return(spw)
}

convert_dctmlefit_to_vector <- function(dctmlefit) {
  
  # Extract ATE values
  ate_values <- unlist(dctmlefit$ATE)
  
  # Extract weight values
  weight_values <- as.vector(unlist(c(dctmlefit$weight[1,-1],
                                      dctmlefit$weight[2,-1]))) 
  
  # Combine ATE and weight values into a single vector
  combined_vector <- c(ate_values, weight_values)
  
  # Attach informative names to the vector elements
  colnames_ate <- colnames(dctmlefit$ATE)
  colnames_weight <- colnames(dctmlefit$weight)[-1] 
  combined_names <- c(colnames_ate, paste0("x_", colnames_weight), paste0("y_", colnames_weight))
  names(combined_vector) <- combined_names
  
  return(combined_vector)
}


OR2RD <- function(OR,datax = full.data, p0 = NULL){
  if (is.null(p0)){
    #load(file = "~/GitHub/DR-hdps/sim/data/plasmode.RData")
    tab <- table(datax$exposure==0,datax$outcome)
    tab.background.exposure <- tab[rownames(tab) == TRUE,]
    tab.background.exposure.by.outcome <- tab.background.exposure/sum(tab.background.exposure)
    p0 <- as.numeric(tab.background.exposure.by.outcome["1"])
  }
  # p0 = risk of outcome in unexposed population
  O0 = p0/(1-p0)
  RD = O0*(OR-1)/( (1+OR*O0)*(1+O0) )
  return(list(RD=RD, p0 = p0))
}

remain.x <- function(pathx = "~/GitHub/DR-hdPS-sim/result/tmle/DCtmleS/"){
  setwd(pathx)
  nnx <- list.files(pattern = ".Rds", full.names=T)
  `%ni%`<-Negate(`%in%`)
  m <- gregexpr('[0-9]+',nnx)
  donen <- as.character(regmatches(nnx,m))
  alln <- stringr::str_pad(as.character(1:1000), 3, pad = "0")
  xxx <- as.numeric(as.character(alln[alln %ni% donen]))
  # x<-list.files(pattern = ".Rds", full.names=T)
  # setdiff(1:1000, as.numeric(gsub("\\D", "", x)))
  # x<-list.files(pattern = ".Rds", full.names=T)
  # head(x)
  # x2<-stringr::str_sub(x, start=5, -6)
  # head(x2)
  # y <- setdiff(1:1000, as.numeric(gsub("\\D", "", x2)))
  # length(y)
  # edit(y)
  return(xxx)
}
#1000-length(remain.x())

