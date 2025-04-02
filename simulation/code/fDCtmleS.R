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
# require(gam)
require(kernlab)
require(xgboost)
require(purrr)
# require(rJava)
# require(bartMachine)
# require(BayesTree)
# require(dbarts)
# library(Crossfit) 

# plasmode.data = simdata 
# analytic.data = full.data 
# demovars = demovars
# lab.list = labvars.original
# proxy.list = emp.cov.names
# nSim = nSim
# method = method.under.consider
# i=1

tmle_single_g1_p = function(data,
                            exposure,
                            outcome,
                            covarsT,
                            covarsO,
                            family.y,
                            learners,
                            control,
                            n_split,
                            rand_split,
                            gbounds,
                            Qbounds,
                            seed=146){
  
  
  #### step: 1.1 ############
  
  set.seed(seed)
  splits_p <- sample(rep(1:n_split, diff(floor(nrow(data) * c(0:n_split/n_split)))))
  
  if(family.y == "gaussian"){
    original.Y = outcome
    min.Y = min(data[, outcome])
    max.Y = max(data[, outcome])
    data$Y.star = (data[, outcome] - min.Y)/(max.Y-min.Y)
    outcome = "Y.star"
  }
  
  data_pp = data_p = data %>% mutate(s=splits_p)  %>% arrange(s)
  
  # Create nested dataset
  
  dat_nested_p <- data_p %>%
    group_by(s) %>%
    tidyr::nest()
  
  #### step: 1.2.1 ############
  
  # P-score model
  
  pi_fitter <- function(df){
    SuperLearner::SuperLearner(Y=as.matrix(df[, exposure]),
                               X=df[, covarsT],
                               family=binomial(),
                               SL.library=learners,
                               cvControl=control)
  }
  
  dat_nested_p <- dat_nested_p %>% mutate(pi_fit = map(data, pi_fitter))
  
  
  ########## 1.3.1 ##########################
  pi = H1 = H0 = list()
  k = dim(data_p)[2]
  for(i in 1:n_split){
    pi[[i]] = predict(dat_nested_p$pi_fit[[i]], newdata = data_pp[, covarsT])$pred
    pi[[i]] = ifelse(pi[[i]] < gbounds, gbounds, ifelse(pi[[i]] > (1-gbounds), (1-gbounds), pi[[i]]))
    H1[[i]] = data_pp[, exposure]/pi[[i]]
    H0[[i]] = (1 - data_pp[, exposure])/(1 - pi[[i]])
    data_p = suppressMessages(bind_cols(data_p, pi[[i]], H1[[i]], H0[[i]]))
    
  }
  
  names(data_p) <-   c(names(data_p)[1:k], paste0(rep(c("pi", "H1_", "H0_"),
                                                      times = n_split), rep(1:n_split, each = 3)))
  
  #Outcome model
  
  #### step: 1.2.2 ############
  
  if(family.y == "binomial"){
    mu_fitter <- function(df){
      SuperLearner::SuperLearner(Y=as.matrix(df[, outcome]),
                                 X=df[, c(exposure, covarsO)],
                                 family=binomial(),
                                 SL.library=learners,
                                 cvControl=control)
    }
  }
  
  if(family.y == "gaussian"){
    mu_fitter <- function(df){
      SuperLearner::SuperLearner(Y=as.matrix(df[, outcome]),
                                 X=df[, c(exposure, covarsO)],
                                 family="gaussian",
                                 SL.library=learners,
                                 cvControl=control)
    }
  }
  
  
  
  dat_nested_p <- dat_nested_p %>% mutate(mu_fit=map(data, mu_fitter))
  
  
  # Calculation of mu using each split
  
  dat1_p = dat0_p = data_pp
  
  dat1_p[, exposure] = 1
  
  dat0_p[, exposure] = 0
  
  k1 = dim(data_p)[2]
  
  mu_name = c(paste0(rep(c("mu_", "mu1_", "mu0_"), times = n_split), rep(1:n_split, each = 3)))
  mu = mu1 = mu0 = list()
  
  ########## 1.3.2 ##########################
  
  if(family.y == "binomial"){
    for(i in 1:n_split){
      mu[[i]] = predict(dat_nested_p$mu_fit[[i]], newdata = data_pp[, c(exposure, covarsO)], type = "response")$pred
      mu[[i]] = ifelse(mu[[i]] < Qbounds, Qbounds, ifelse(mu[[i]] > 1-Qbounds, 1-Qbounds, mu[[i]]))
      mu1[[i]] = predict(dat_nested_p$mu_fit[[i]], newdata = dat1_p[, c(exposure, covarsO)], type = "response")$pred
      mu1[[i]] = ifelse(mu1[[i]] < Qbounds, Qbounds, ifelse(mu1[[i]] > 1-Qbounds, 1-Qbounds, mu1[[i]]))
      mu0[[i]] = predict(dat_nested_p$mu_fit[[i]], newdata = dat0_p[, c(exposure, covarsO)], type = "response")$pred
      mu0[[i]] = ifelse(mu0[[i]] < Qbounds, Qbounds, ifelse(mu0[[i]] > 1-Qbounds, 1-Qbounds, mu0[[i]]))
      data_p = suppressMessages(bind_cols(data_p, mu[[i]], mu1[[i]], mu0[[i]]))
    }
  }
  if(family.y == "gaussian"){
    for(i in 1:n_split){
      mu[[i]] = predict(dat_nested_p$mu_fit[[i]], newdata = data_pp[, c(exposure, covarsO)], type = "response")$pred
      mu[[i]] = ifelse(mu[[i]] < Qbounds, Qbounds, ifelse(mu[[i]] > 1-Qbounds, 1-Qbounds, mu[[i]]))
      mu1[[i]] = predict(dat_nested_p$mu_fit[[i]], newdata = dat1_p[, c(exposure, covarsO)], type = "response")$pred
      mu1[[i]] = ifelse(mu1[[i]] < Qbounds, Qbounds, ifelse(mu1[[i]] > 1-Qbounds, 1-Qbounds, mu1[[i]]))
      mu0[[i]] = predict(dat_nested_p$mu_fit[[i]], newdata = dat0_p[, c(exposure, covarsO)], type = "response")$pred
      mu0[[i]] = ifelse(mu0[[i]] < Qbounds, Qbounds, ifelse(mu0[[i]] > 1-Qbounds, 1-Qbounds, mu0[[i]]))
      data_p = suppressMessages(bind_cols(data_p, mu[[i]], mu1[[i]], mu0[[i]]))
    }
  }
  
  
  names(data_p) = c(names(data_p)[1:k1], mu_name)
  
  
  pi_p <- suppressWarnings(data_p %>%
                             select(s, paste0("pi", 1:n_split)))
  
  mu_p <- suppressWarnings(data_p %>%
                             select(s, paste0("mu_", 1:n_split)))
  
  mu1_p <- suppressWarnings(data_p %>%
                              select(s, paste0("mu1_", 1:n_split)))
  
  mu0_p <- suppressWarnings(data_p %>%
                              select(s, paste0("mu0_", 1:n_split)))
  
  Y_p <-  suppressWarnings(data_p %>%
                             select(s, outcome))
  
  
  X_p <- suppressWarnings(data_p %>%
                            select(s, exposure))
  
  H_p <-  data_p %>%
    select(s, paste0(rep(c("H0_", "H1_"), times = n_split), rep(1:n_split, each = 2)))
  
  
  ################## Step: 1.3.3 #######################################
  
  epsilon = iid = list(); mu0_1 = mu1_1 = mu0_1s = mu1_1s = mu0_1ss = mu1_1ss = list()
  if(rand_split == TRUE){
    for(i in 1:n_split){
      iid[[i]] = sample(setdiff(1:n_split, i) , 2, replace = FALSE)
      
      h0 = H_p[ , paste0("H0_", iid[[i]][1])]
      
      h1 = H_p[ , paste0("H1_", iid[[i]][1])]
      
      muu = mu_p[ , paste0("mu_", iid[[i]][2])]
      
      y =  Y_p[ , outcome]
      
      dd <- data.frame(y=y, h0=h0, h1=h1, muu=muu)
      
      epsilon[[i]] <- coef(glm(dd[,1] ~ -1 + dd[,2] + dd[,3] + offset(qlogis(dd[,4])), family = binomial()))
      
      mu0_1s[[i]] = plogis(qlogis(pull(mu0_p, paste0("mu0_",  iid[[i]][2]))) + epsilon[[i]][1] / (1 - pull(pi_p, paste0("pi", iid[[i]][1]))))
      
      mu1_1s[[i]] = plogis(qlogis(pull(mu1_p, paste0("mu1_", iid[[i]][2]))) + epsilon[[i]][2] / pull(pi_p, paste0("pi", iid[[i]][1])))
      
    }
    
  }
  
  
  if(rand_split == FALSE){
    pi_id = c(2:n_split, 1); mu_id = c(3:n_split, 1, 2)
    for(i in 1:n_split){
      
      h0 = H_p[, paste0("H0_", pi_id[i])]
      
      h1 = H_p[, paste0("H1_", pi_id[i])]
      
      muu = mu_p[, paste0("mu_", mu_id[i])]
      
      y =  Y_p[, outcome]
      
      dd <- data.frame(y=y, h0=h0, h1=h1, muu=muu)
      
      epsilon[[i]] <- coef(glm(dd[,1] ~ -1 + dd[,2] + dd[,3] + offset(qlogis(dd[,4])), family = binomial()))
      
      mu0_1s[[i]] = plogis(qlogis(pull(mu0_p, paste0("mu0_",  mu_id[i]))) + epsilon[[i]][1] / (1 - pull(pi_p, paste0("pi", pi_id[i]))))
      
      mu1_1s[[i]] = plogis(qlogis(pull(mu1_p, paste0("mu1_", mu_id[i]))) + epsilon[[i]][2] / pull(pi_p, paste0("pi", pi_id[i])))
      
    }
  }
  
  mu0_1_p = suppressMessages(bind_cols(mu0_1s))
  mu1_1_p = suppressMessages(bind_cols(mu1_1s))
  
  dat_name = names(data_p)
  
  data_p = suppressMessages(bind_cols(data_p, mu0_1_p, mu1_1_p))
  
  
  
  names(data_p) = c(dat_name, paste0("mu0_1_", 1:n_split), paste0("mu1_1_", 1:n_split))
  
  if(family.y == "gaussian"){
    for(i in 1:n_split){
      data_p[ ,paste0("mu1_1_", i)] <- data_p[ ,paste0("mu1_1_", i)]*(max.Y-min.Y) + min.Y
      data_p[ ,paste0("mu0_1_", i)] <- data_p[ ,paste0("mu0_1_", i)]*(max.Y-min.Y) + min.Y
    }
    outcome = original.Y
  }
  
  
  r1 = r0 = rd = NULL
  r1_0 <- list()
  if1 = if0 = ifd = list()
  v1 = v0 = vd = NULL
  rdc <- NULL
  
  mu_summ <- list(); n <- NULL
  
  ############### step 1.3.4 - 1.3.5 ######################
  
  if(rand_split == FALSE){
    for(i in 1:n_split){
      data_p$ss <- data_p$s
      mu_summ[[i]] <- data_p %>% filter(s %in% i) %>%
        select(exposure, outcome, paste0("pi", i), paste0("mu1_1_", i), paste0("mu0_1_", i), ss)
      names(mu_summ[[i]]) <- c("x", "y", "pi", "mu1", "mu0", "s")
      n[i] <- nrow(mu_summ[[i]])
    }
  }
  
  
  if(rand_split == TRUE){
    
    for(i in 1:n_split){
      data_p$ss <- data_p$s
      mu_summ[[i]] <- data_p %>% filter(s %in% i) %>%
        select(exposure, outcome, paste0("pi", i), paste0("mu1_1_", i), paste0("mu0_1_", i), ss)
      names(mu_summ[[i]]) <- c("x", "y", "pi", "mu1", "mu0", "s")
      n[i] <- nrow(mu_summ[[i]])
    }
  }
  
  
  mu_sum <- bind_rows(mu_summ)
  mu_sum$ss <- rep(1:n_split, times = n)
  
  
  
  
  ############ step 1.6 ########################
  
  for(i in 1:n_split){
    rdc[i] <- mean(mu_sum[mu_sum$ss == i, ]$mu1) - mean(mu_sum[mu_sum$ss == i, ]$mu0)
    D1 = with(mu_sum[mu_sum$ss == i, ], (x/pi)*(y-mu1) + (mu1 - mean(mu1)))
    D0 = with(mu_sum[mu_sum$ss == i, ], ((1-x)/(1-pi))*(y-mu0) + (mu0 - mean(mu0)))
    EIC = D1-D0
    vd[i] <- var(EIC)/nrow(data)
  }
  
  
  rd = mean(rdc)
  var1 = mean(vd)
  
  res <- data.frame(rd=rd, var = var1)
  
  return(res)
  
}

DC_tmle_g1_k <- function(data,
                         exposure,
                         outcome,
                         covarsT,
                         covarsO,
                         family.y="binomial",
                         learners,
                         control,
                         n_split,
                         num_cf,
                         rand_split=FALSE,
                         gbounds = NULL,
                         Qbounds = 5e-04,
                         seed=146,
                         conf.level=0.95,
                         stat = "median"){
  
  runs <- list()
  # placeholder_output <- generate_placeholder_output(learners, n_split)
  #Run on num_cf splits
  set.seed(seed)
  cf_seed = sample(num_cf)
  n = nrow(data)
  if(is.null(gbounds)) gbounds = 5/sqrt(n)/log(n)
  #################### step 2 and 3 ########################
  
  for(cf in 1:num_cf){
    seed1 = cf_seed[cf]
    fit_result = try({
      tmle_single_g1_p(data,
                       exposure,
                       outcome,
                       covarsT,
                       covarsO,
                       family.y,
                       learners,
                       control,
                       n_split,
                       rand_split,
                       gbounds,
                       Qbounds,
                       seed=seed1)}, silent = TRUE)
    
    if (inherits(fit_result, "try-error")) {
      fit_sngle <-  data.frame(rd=NA, var = NA)
    } else {
      fit_sngle <- fit_result
    }
    
    runs[[cf]] <- fit_sngle
  }
  
  res = dplyr::bind_rows(runs)
  
  if(stat == "mean"){
    medians <- apply(res, 2, mean, na.rm = TRUE)
    res <- res %>% mutate(var0 = var + (rd - medians[1])^2)
    #  results <- apply(res, 2, mean, na.rm = TRUE)
  }
  if(stat == "median"){
    medians <- apply(res, 2, median, na.rm = TRUE)
    res <- res %>% mutate(var0 = var + (rd - medians[1])^2)
    #  results <- apply(res, 2, median, na.rm = TRUE)
  }
  t.value = qt((1-conf.level)/2, nrow(data), lower.tail = F)
  
  #l_ci = results[1] - t.value*sqrt(results[3])
  #u_ci = results[1] + t.value*sqrt(results[3])
  
  #res1 = tibble(ATE=results[1], se = sqrt(results[3]), lower.ci = l_ci, upper.ci = u_ci)
  
  return(res)
}

loop.to.results.gen<- function(plasmode.data = simdata, 
                               analytic.data = full.data, 
                               demovars = demovars,
                               lab.list = labvars.original,
                               proxy.list = emp.cov.names,
                               nSim = nSim,
                               method = method.under.consider,
                               i=1){
  res <- rep(NA,16)
  names(res) <- c("rd", "se", "lower.ci", "upper.ci", "x_SL.glm_All", "x_SL.earth_All", 
                  "x_SL.glmnet_All", "x_SL.ksvm_All", "x_SL.xgboost_All", "x_prev", 
                  "y_SL.glm_All", "y_SL.earth_All", "y_SL.glmnet_All", "y_SL.ksvm_All", 
                  "y_SL.xgboost_All", "y_prev")
  
  try({
    options(warn=-1)
    plasmodeData.i <- plyr::join(x = data.frame(idx=plasmode.data[,i],
                                                EVENT=plasmode.data[,i+nSim]), 
                                 y = analytic.data, by="idx", type="left")
    plasmodeData.i$outcome <- plasmodeData.i$EVENT
    
    SL.library <- c("SL.glm", 
                    "SL.earth", 
                    "SL.glmnet")
    
    n.fd <- nrow(plasmodeData.i)
    covarsT <- c(demovars, lab.list, proxy.list)
    covarsO <- c(demovars, lab.list, proxy.list)
    dctmlefit <- DC_tmle_g1_k(data=plasmodeData.i,
                               exposure="exposure",
                               outcome="outcome",
                               covarsT=covarsT,
                               covarsO=covarsO,
                               family.y="binomial",
                               learners=SL.library, # c("SL.glm.dcTMLE", "SL.mean.dcTMLE", "SL.gam4.dcTMLE", "SL.xgboost.dcTMLE"),
                               control=list(V = 3, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                               num_cf=25, # 100, # 5,
                               n_split=3,
                               rand_split=FALSE,
                               Qbounds = 0.0025,
                               gbounds = 5/sqrt(n.fd)/log(n.fd), # 0.0025,
                               seed=111,
                               conf.level=0.95,
                               stat = "median")
    res <- dctmlefit
    # colnames(res) <- c("rd", "se", "lower.ci", "upper.ci")
    saveRDS(res, file = paste0(method.under.consider,"/sim",
                               stringr::str_pad(i, 3, pad = "0"),
                               "res.Rds"))
  })
  return(res)
}