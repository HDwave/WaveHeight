# Performs Box-Cox transformation. The transformation parameter (lambda) is chosen to minimize deviation from the normal distribution (minumum sum of squared skewness and squared kurtosis)
# @param obs The observations to be transformed
# @return The observations transformed
BoxCoxLambda = function(obs) {
	#lambda = seq(-1.0, 1.0, 0.1)
	lambda = seq(0.0, 1.0, 0.1)
	#lambda = seq(-0.5, 1.0, 0.1)
	obs.trans = matrix(nrow = length(lambda), ncol = length(obs))
	normdev = rep(NA, length(lambda)) # Holds the amount of deviation from the normal distribution
    #indNA <- which(is.na(obs))
    #obs
	for(i in 1:length(lambda)) {
		if(lambda[i] == 0) {
			obs.trans[i,] = log(obs)
		} else {
			obs.trans[i,] = (obs^lambda[i]-1)/lambda[i]
		}
		normdev[i] = skewness(obs.trans[i,],na.rm = TRUE)^2 + 0.25*(kurtosis(obs.trans[i,],na.rm = TRUE))^2
	}
	return(list(data=obs.trans[which.min(normdev),],lambda = lambda[which.min(normdev)]))
}
