require(glmnet)
require(forecast)

## Import functions to run the experiments below
source("getPreddistr.r")
source("BoxCoxLambda.r")
source("BoxCoxLambdaKnown.r")
source("fourier.r")
source("qBoxCox.r")
source("pBoxCox.r")

## Load ERA-Interim data 
load("../data/SWH_SLP_SLP.grad.subset.7.Rdata")

longitudeSWH
latitudeSWH
range(years.all)
summary(SWH)
summary(SLP)
summary(SLP.grad)

## Which grid cells to predict
lonSWH = 4
latSWH = 4

## How many missing values is "good enough"
na.thresh = 500

## Define the model
spatial.neighborhoods = 2 
models.fourier = c("Proposed model, AR order = 5, spatial neigb. = ", spatial.neighborhoods, "\n")
n.models.fourier = length(models.fourier)

## Create Fourier terms 
m.fourier = 365.25*4 # Number of observations per year. Four observations per day.
K.fourier = 2 # Number of Fourier terms
n.time = length(time.all)
intercept.fourier = fourier(rep(1,n.time), K=K.fourier, m=m.fourier)
colnames(intercept.fourier) = paste("intercept", colnames(intercept.fourier), sep = "_")

## Split that data in a training set and a test set 
training.test = list()
training.test[[1]] = which(years.all %in% 2006:2014)
training.test[[2]] = which(years.all %in% 2015)

## Estimate model parameters and obtain predictive distributions
pred.dist  <- getPreddistr(j=lonSWH, k=latSWH, neig=spatial.neighborhoods)
print(pred.dist$fits)

## Save predictive distributions
pred.mean = pred.dist$pred.mean
pred.sd = pred.dist$pred.sd
pred.lambda = pred.dist$pred.lambda

## Get a vector of observations in the test period
obs  <- SWH[lonSWH, latSWH, training.test[[2]]]

## Calculate PIT values and plot a PIT histogram for the test set 
pit  <- pBoxCox(obs, pred.mean, pred.sd, pred.lambda)
hist(pit, freq=FALSE, nclass=10, col="gray", xlab="PIT value", 
     main=paste("Lon = ", longitudeSWH[lonSWH], ", Lat = ", latitudeSWH[latSWH]))
abline(a=1, b=0, lty=2)

## Calculate the mean absolute error over the test period
pred.median  <- qBoxCox(0.5, pred.mean, pred.sd, pred.lambda)
mae  <- mean(abs(obs - pred.median))
mae

## Plot the prediction and the observation in the first and last 100 time points of the test period
upper  <- qBoxCox(0.95, pred.mean, pred.sd, pred.lambda)
lower  <- qBoxCox(0.05, pred.mean, pred.sd, pred.lambda)

t.period  <- c(1:100)
t.upper  <- upper[t.period]
t.lower  <- lower[t.period]

plot(t.period, obs[t.period], type="l", xlab="Time point in test period", ylab="SWH", ylim=c(0,15),
     main=paste("Lon = ", longitudeSWH[lonSWH], ", Lat = ", latitudeSWH[latSWH]))
polygon(c(t.period, rev(t.period), t.period[1]), c(t.lower, rev(t.upper), t.lower[1]), col="gray90", border="NA")
lines(t.period, obs[t.period])
lines(t.period, pred.median[t.period], col="gray50", lty=2)
legend("topright", lty=c(1,2,1), lwd=c(1,1,4), col=c("black", "gray50", "gray90"),
       legend=c("Observation", "Predictive median", "90% prediction interval"))

t.period  <- c(1361:1460)
t.upper  <- upper[t.period]
t.lower  <- lower[t.period]

plot(t.period, obs[t.period], type="l", xlab="Time point in test period", ylab="SWH", ylim=c(0,15),
     main=paste("Lon = ", longitudeSWH[lonSWH], ", Lat = ", latitudeSWH[latSWH]))
polygon(c(t.period, rev(t.period), t.period[1]), c(t.lower, rev(t.upper), t.lower[1]), col="gray90", border="NA")
lines(t.period, obs[t.period])
lines(t.period, pred.median[t.period], col="gray50", lty=2)
legend("topright", lty=c(1,2,1), lwd=c(1,1,4), col=c("black", "gray50", "gray90"),
       legend=c("Observation", "Predictive median", "90% prediction interval"))

## Predictive trajectories for 10 time points
t.ind  <- c(40:49)
random.q  <- array(NA, dim=c(10,10))
for(i in 1:10) random.q[,i]  <- qBoxCox(runif(10), pred.mean[t.ind[i]], pred.sd, pred.lambda)
random.q
plot(t.ind, random.q[1,], type="l", col="gray50",
     xlab="Time point in test period", ylab="SWH", ylim=c(5,12),
     main=paste("Lon = ", longitudeSWH[lonSWH], ", Lat = ", latitudeSWH[latSWH]))
for(i in 2:10) lines(t.ind, random.q[i,], col="gray50")
lines(t.ind, obs[t.ind], col="black", lwd=2)

## Learn correlation from previous timepoints (last 100 time points in test period)
sample.q  <- array(NA, dim=c(10,10))
for(i in 1:10) sample.q[,i]  <- rank(random.q[,i])
sample.q
T  <- length(training.test[[2]])
h.ind  <- training.test[[2]][(T-99):T]
hist.obs  <- t(array(SWH[lonSWH, latSWH, h.ind], dim=c(10,10)))
hist.q  <- array(NA, dim=c(10,10))
for(i in 1:10) hist.q[,i]  <- rank(hist.obs[,i])
hist.q
sort.q  <- random.q
for(i in 1:10) sort.q[,i]  <- sort(random.q[,i])[hist.q[,i]]
plot(t.ind, sort.q[1,], type="l", col="gray50",
     xlab="Time point in test period", ylab="SWH", ylim=c(5,12),
     main=paste("Lon = ", longitudeSWH[lonSWH], ", Lat = ", latitudeSWH[latSWH]))
for(i in 2:10) lines(t.ind, sort.q[i,], col="gray50")
lines(t.ind, obs[t.ind], col="black", lwd=2)
