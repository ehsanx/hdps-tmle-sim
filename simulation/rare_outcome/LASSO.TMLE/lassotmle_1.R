require(parallel)
require(plyr)

method.under.consider <- "lassotmle1"
mainPath <- paste0(getwd(), "/")
if (! dir.exists(paste0(mainPath, paste0(method.under.consider,"/") ))){
  dir.create(paste0(mainPath, paste0(method.under.consider,"/") ), recursive = TRUE)
}
if (! dir.exists(paste0(mainPath, paste0("weight", method.under.consider,"/") ))){
  dir.create(paste0(mainPath, paste0("weight", method.under.consider,"/") ), recursive = TRUE)
}
if (! dir.exists(paste0(mainPath, paste0("VarSel", method.under.consider,"/") ))){
  dir.create(paste0(mainPath, paste0("VarSel", method.under.consider,"/") ), recursive = TRUE)
}
if (! dir.exists(paste0(mainPath, "result"))){
  dir.create(paste0(mainPath, "result"), recursive = TRUE)
}

load(file = "../../data/plasmode_outcome_rare.RData") 
rm(simdata.obj)
rm(simdata.obj0)
rm(basicdata)
rm(extendeddata)
numCores <- detectCores()
numCoresUsed <- numCores
iterStart <- 1
iterStop <- 20
cl <- makeCluster(numCoresUsed)
clusterExport(cl, "simdata", envir=.GlobalEnv)
clusterExport(cl, "full.data", envir=.GlobalEnv)
clusterExport(cl, "demovars", envir=.GlobalEnv)
clusterExport(cl, "labvars.original", envir=.GlobalEnv)
clusterExport(cl, "emp.cov.names", envir=.GlobalEnv)
clusterExport(cl, "nSim", envir=.GlobalEnv)
clusterExport(cl, "method.under.consider", envir=.GlobalEnv)
clusterEvalQ(cl,{
  #source the simulation code
  source(paste0("../../code/f",method.under.consider,".R"), local = TRUE)
})
start_time <- Sys.time()
results <- parLapply(cl, iterStart:iterStop, fun = function(j){ 
  loop.to.results.gen(plasmode.data = simdata, 
                      analytic.data = full.data, 
                      proxy.list = emp.cov.names,
                      lab.list = labvars.original,
                      demovars = demovars,
                      nSim = nSim,
                      method = method.under.consider,
                      i=j)
})
stopCluster(cl)

resultsDF <- data.frame()
for (i in iterStart:iterStop) {
  results_name <- paste0("sim", stringr::str_pad(i, 3, pad = "0"), "res.Rds")
  result <- readRDS(paste0(method.under.consider, "/", results_name))
  resultsDF <- rbind(resultsDF, result)
  if (i == iterStop) {
    colnames(resultsDF) <- names(result)
  }
}