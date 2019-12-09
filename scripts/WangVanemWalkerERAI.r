require(parallel)
require(glmnet)
require(forecast)
require(scoringRules)
require(Rcpp)
require(moments)

#options("width"=220)
#options("width"=180)

## Import functions to run the experiments below
source("functions.r")

n.training.test = 1

area.longSWH = 6:8 #1:181
area.latSWH = 6:8 # 1:68

mc.cores = 1 # min(length(area.longSWH), 60)

na.thresh = 500

runFourier = TRUE
if(runFourier) {    
    # Data for the seasonal models
    load("~/Documents/Data/HDwave/SWH_SLP_SLP.grad.subset.15.Rdata")
    # load("~/Documents/Data/HDwave/longlatSWHSLP.Rdata")

    spatial.neighborhoods = c(2) # c(5) # For the Vanem&Walker model (max, min, mean).
    models.fourier = c("Proposed model, AR order = 5, spatial neigb. = 2")
    n.models.fourier = length(models.fourier)

    # For the Fourier transforms
    m.fourier = 365.25*4 # Number of observations per year. Four observations per day.
    K.fourier = 2 # Number of Fourier periodics
    n.time = length(time.all)

    intercept.fourier = fourier(rep(1,n.time), K=K.fourier, m=m.fourier)
    colnames(intercept.fourier) = paste("intercept", colnames(intercept.fourier), sep = "_")

    training.test = list()
    part = list()

    part[[1]] = which(years.all %in% 2006:2014)
    part[[2]] = which(years.all %in% 2015)
    training.test[[1]] = part

    preddistrj <- mclapply(area.longSWH, get.preddistr, mc.cores = mc.cores, mc.preschedule=FALSE)
    save(preddistrj, file = "~/HDwaveData/Results/preddistrj.to.Hanne.Rdata")

    preddistr.mean = array(dim = c(length(area.longSWH), length(area.latSWH), length(training.test[[1]][[2]])))
    preddistr.sd = array(dim = c(length(area.longSWH), length(area.latSWH)))
    preddistr.lambda = array(dim = c(length(area.longSWH), length(area.latSWH)))
    
    for(j in 1:length(area.longSWH)) {
        cat("j =", j, "\n")        
        preddistr.mean[j,,] = preddistrj[[j]][[1]]
        preddistr.sd[j,] = preddistrj[[j]][[2]]
        preddistr.lambda[j,] = preddistrj[[j]][[3]]
    }
    save(preddistr.mean, preddistr.sd, preddistr.lambda, file = "/shared/preddistr.to.Hanne.Rdata")
}
