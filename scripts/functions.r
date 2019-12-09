merge.acf = function(acfj) {
    created = FALSE
    n.acf = 0
    for(i in 1:length(acfj)) {
        if(i%%1e5 == 0) {
            cat("i =", i, "i/length(acfj)*100 =", i/length(acfj)*100, "\n")
        }
        if(length(acfj[[i]]) >= 2) {
            if(created) {
                len1 = length(acf.)
                len2 = length(acfj[[i]])
                len = min(len1, len2)
                acf. = acf.[1:len]
                add = acfj[[i]][1:len]
                acf. = acf. + add
                n.acf = n.acf + 1
            } else {
                acf. = acfj[[i]]
                created = TRUE
                n.acf = n.acf + 1
            }
        }
    }
    acf. = acf./n.acf
    return(acf.)
}

# Computes autocorrelation for one traice in one dimension for three dimensional array
# j - The index to the two other directions
# td - Trace to compute autocorrelation
ac3 = function(j, td, min.trace.length) {
    dex = (1:3)[-td]
    dims = dim(d)
    n1 = dims[dex[1]]
    n2 = dims[dex[2]]
    
    if(j%%n2 == 0) {
        i2 = n2
        i1 = j/n2
     } else {
        i2 = j%%n2
        i1 = (j - j%%n2)/n2 + 1
    }

    #cat("j =", j, "i1 =", i1, "i2 =", i2, "dims =", dims, "\n")
    if(round(td) == 1) {
        trace = d[,i1,i2]
    } else if(round(td) == 2) {
        trace = d[i1,,i2]
    } else { # td == 3
        trace = d[i1,i2,]
    }
    
    if(sum(!is.na(trace)) > 0) {
        trace = as.numeric(na.contiguous(trace))
        if(length(trace) > min.trace.length) {
            lag.max = min(1e3, length(trace))
            ac = acf(trace, lag.max = lag.max)
            retur = c(ac$acf)
        } else {
            retur = NA
        }
    } else {
        retur = NA
    }
    return(retur)
}

pSupport = function(mean, sd, lambda) {
    if(lambda > 1e-6) {
        p = pnorm(-1/lambda, mean = mean, sd = sd)
    } else if(lambda < -1e-6) {
        p = 1-pnorm(-1/lambda, mean = mean, sd = sd)
    } else { # lognormal
        p = rep(0, length(mean))
    }
    return(p)
}


ppredbc = function(obs.bc, mean, sd, lambda) {
    if(round(lambda,2) < 0) {
        p = pnorm(obs.bc, mean = mean, sd = sd)/pnorm(-1/lambda, mean = mean, sd = sd)
    } else if(round(lambda,2) == 0) {
        p = pnorm(obs.bc, mean = mean, sd = sd)
    } else { # lambda > 0
        p = pnorm(obs.bc, mean = mean, sd = sd)/(1-pnorm(-1/lambda, mean = mean, sd = sd))
    }
    return(p)
}

qboxcox = function(p, mean, sd, lambda) {
    if(round(lambda,2) < 0) {
        q = (lambda*(sd*qnorm(p)+mean)+1)^(1/lambda)
    } else if(round(lambda,2) == 0) {
        q = exp(mean + sd*qnorm(p))
    } else { # lambda > 0
        T = (1/lambda + mean)/sd
        Vp = 1 - (1-p)*pnorm(T)
        q = (lambda*(sd*qnorm(Vp)+mean)+1)^(1/lambda)    
    }
    return(q)
}


dboxcox.old = function(x, mean, sd, lambda, log = FALSE) {
    val = rep(NA, length(x))
    if(log == FALSE) {
        if(round(lambda,2) < 0) {
            val[x>0] = 1/(sqrt(2*pi)*sd*pnorm(-1/(lambda*sd)-mean/sd))*exp(-1/(2*sd^2)* ((x[x>0]^lambda-1)/lambda - mean)^2)*abs(x[x>0]^(lambda-1))
            val[x<=0] = 0    
        } else if(round(lambda,2) == 0) {
            val = dlnorm(x, mean, sd, log = FALSE)
        } else { # lambda > 0
            val[x>0] = 1/(sqrt(2*pi)*sd*(1-pnorm(-1/(lambda*sd)-mean/sd)))*exp(-1/(2*sd^2)* ((x[x>0]^lambda-1)/lambda - mean)^2)*abs(x[x>0]^(lambda-1))
            val[x<=0] = 0
        }
    } else { # log == TRUE
        if(round(lambda,2) < 0) {
            val[x>0] = -log( sqrt(2*pi)*sd*pnorm(1/(lambda*sd)-mean/sd) ) - 1/(2*sd^2)* ((x[x>0]^lambda-1)/lambda - mean)^2 + log(abs(x[x>0]^(lambda-1)))
            val[x<=0] = -Inf
        } else if(round(lambda,2) == 0) {
            val = dlnorm(x, mean, sd, log = TRUE)
        } else { # lambda > 0
            val[x>0] = -log(sqrt(2*pi)*sd*(1-pnorm(-1/(lambda*sd)-mean/sd))) - 1/(2*sd^2)* ((x[x>0]^lambda-1)/lambda - mean)^2 + log(abs(x[x>0]^(lambda-1)))
            val[x<=0] = -Inf
        }    
    }
    return(val)
}



dboxcox = function(x, mean, sd, lambda, log = FALSE) {
    if(log == FALSE) {
        val = rep(0, length(x))
        if(round(lambda,2) < 0) {
            val[x>0] = x[x>0]^(lambda-1)*sd^(-1)*dnorm(((x[x>0]^lambda-1)/lambda-mean)/sd)*pnorm((-1/lambda-mean)/sd)^(-1)
        } else if(round(lambda,2) == 0) {
            val = dlnorm(x, mean, sd, log = FALSE)
        } else {
            val[x>0] = x[x>0]^(lambda-1)*sd^(-1)*dnorm(((x[x>0]^lambda-1)/lambda-mean)/sd)*pnorm((1/lambda+mean)/sd)^(-1)
        }
    } else {
        val = rep(-Inf, length(x))
        if(round(lambda,2) < 0) {
            val[x>0] = (lambda-1)*log(x[x>0]) - log(sd) + dnorm(((x[x>0]^lambda-1)/lambda-mean )/sd, log = TRUE) - log(pnorm( (-1/lambda - mean)/sd))
        } else if(round(lambda,2) == 0) {
            val = dlnorm(x, mean, sd, log = TRUE)
        } else {
            val[x>0] = (lambda-1)*log(x[x>0]) - log(sd) + dnorm(((x[x>0]^lambda-1)/lambda-mean)/sd, log = TRUE) - log(pnorm((1/lambda+mean)/sd))
        }
    }
    return(val)
}


# Computes the reliable index
# U - The rank values or PIT values normalized to the zero one intervall.
reliabilityIndex = function(U, n.bins) {
    U = U[ !is.na(U) ]
    n = length(U)
    if(n > 0) {
        seps = seq(0, 1, length.out = n.bins+1)
        probs = rep(NA, n.bins)
        for(i in 1:n.bins) {
            probs[i] = sum( U >= seps[i] & U <= seps[i+1] )
        }
        probs = probs/n
        unif.prob = rep(1/n.bins, n.bins)
        RI = mean(abs(probs - unif.prob)*100)
    } else {
        RI = NA
    }
    return(RI)
}

# Computes the reliable index squared
# U - The rank values or PIT values normalized to the zero one intervall.
reliabilityIndexSquare = function(U, n.bins) {
    U = U[ !is.na(U) ]
    n = length(U)
    if(n > 0) {
        seps = seq(0, 1, length.out = n.bins+1)
        probs = rep(NA, n.bins)
        for(i in 1:n.bins) {
            probs[i] = sum( U >= seps[i] & U <= seps[i+1] )
        }
        probs = probs/n
        unif.prob = rep(1/n.bins, n.bins)
        RI = mean((probs - unif.prob)^2*100)
    } else {
        RI = NA
    }
    return(RI)
}
   
fourier = function(x, K, m) {
    n = length(x)
    idx = 1:n
    fourier.x = matrix(nrow = n, ncol = 2*K)
    coln = rep(NA, 2*K)
    for(k in 1:K) {
        fourier.x[,2*k-1] = sin(2*pi*k*idx/m)*x
        fourier.x[,2*k] = cos(2*pi*k*idx/m)*x 
        coln[2*k-1] = paste("sin", k, sep = "")
        coln[2*k] = paste("cos", k, sep = "")
    }
    colnames(fourier.x) = coln
    return(fourier.x)
}

plotMap <- function(field,latitude,longitude,height = 7,width = 9,printToFile = FALSE,fname = "plot.pdf",cmap = "jet",landcolor = "wheat",seacolor = "aliceblue",xlimit = NULL,ylimit = NULL,cbarTitle = "Intensity",cbarLegendName = "Intensity",truncateMinMax = NULL){
  require(ggplot2)
  require(reshape)
  require(maps)
  require(PBSmapping)
  require(data.table)
  #lats <- c(39,42,57)
  #longs <- c(-66,-13.5,-36)
  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

  wrld <- map_data("world")
  setnames(wrld, c("X","Y","PID","POS","region","subregion"))

  if(is.null(xlimit)){
    xlimit <- c(-180,180)
    ylimit <- c(-90,90)
  } else {
    wrld = clipPolys(wrld, xlim=xlim,ylim=ylim, keepExtra=TRUE)
  }

  ##Truncate field
  if(!is.null(truncateMinMax)){
      field[which(field < truncateMinMax[1])] = truncateMinMax[1]
      field[which(field > truncateMinMax[2])] = truncateMinMax[2]
  }
  
  mpr <- melt(field[append(which(longitude > 180),which(longitude<=180)),])

  #theme_set(theme_bw())
  
  biglat <- array(NA,length(latitude)*length(longitude))
  biglon <- biglat
  value <- biglat
  counter = 1
  for(i in 1:length(latitude)){
    for(j in 1:length(longitude)){
      biglat[counter] <- latitude[i]
      biglon[counter] <- longitude[j]
      value[counter] <- mpr[counter,"value"]
      counter = counter + 1
    }
  }

  df <- data.frame(Longitude = biglon,Latitude = biglat,Intensity = value)

  if(printToFile){
    pdf(file = fname,height = height,width = width)
    p <- ggplot() + coord_map(xlim = xlimit,ylim = ylimit) + geom_raster(aes(x = Longitude,y = Latitude,fill = Intensity),data = df) + scale_fill_gradientn(colours = jet.colors(7),na.value = "transparent",name = cbarLegendName) + geom_polygon(fill = landcolor,col = landcolor,aes(x = X,y = Y,group = PID),wrld) + coord_equal() + theme(panel.background = element_rect(fill = seacolor))
    print(p)
    dev.off()
  } else {
    X11(height = 7,width = 9)

    if(cmap=="jet"){
        ggplot() + coord_map(xlim = xlimit,ylim = ylimit) + geom_raster(aes(x = Longitude,y = Latitude,fill = Intensity),data = df) + scale_fill_gradientn(colours = jet.colors(7),na.value = "transparent",name = cbarLegendName) + geom_polygon(fill = landcolor,col = landcolor,aes(x = X,y = Y,group = PID),wrld) + coord_equal() + theme(panel.background = element_rect(fill = seacolor))
        #+geom_point(aes(x = longs,y = lats),size = 4,colour = "purple",pch = 18)
    } else {
      ggplot() + geom_raster(aes(x = Longitude,y = Latitude,fill = Intensity),data = df) + scale_fill_gradientn(colours = c("#0000FFFF","#FFFFFFFF","#FF0000FF"),na.value = "transparent",name = cbarLegendName) + geom_polygon(fill = landcolor,col = landcolor,aes(x = X,y = Y,group = PID),wrld) + coord_equal() + xlim(xlimit)+ylim(ylimit)
    }
   
  }
}


InvBoxCox = function(obst, lambda) {
    if(is.na(lambda)) {
        obs = NA
    } else {
        if(lambda == 0) {
            obs = exp(obst)
        } else {
            obs = (lambda*obst + 1)^(1/lambda)
        }
    }
    return(obs)
}

BoxCoxLambdaKnown = function(obs, lambda) {
    obs[round(lambda,2) == 0] = log(obs[lambda == 0])
    obs[round(lambda,2) != 0] = (obs[lambda != 0]^lambda-1)/lambda
    return(obs)
}

BoxCoxLambdaKnown2 = function(obs, lambda) {
    if(round(lambda,3) == 0) {
        obs = log(obs)
    } else {
        obs = (obs^lambda - 1)/lambda
    }
    return(obs)
}


# Performs Box-Cox transformation. The transformation parameter (lambda) is chosen to minimize deviation from the normal distribution (minumum sum of squared skewness and squared kurtosis)
# @param obs The observations to be transformed
# @return The observations transformed
BoxCox = function(obs) {
    lambda = seq(0, 1, 0.1)
    obs.trans = matrix(nrow = length(lambda), ncol = length(obs))
    normdev = rep(NA, length(lambda)) # Holds the amount of deviation from the normal distribution
    for(i in 1:length(lambda)) {
        if(lambda[i] == 0) {
            obs.trans[i,] = log(obs)
        } else {
            obs.trans[i,] = (obs^lambda[i]-1)/lambda[i]
        }
        normdev[i] = skewness(obs.trans[i,],na.rm = TRUE)^2 + 0.25*(kurtosis(obs.trans[i,],na.rm = TRUE))^2
    }
    return(obs.trans[which.min(normdev) ,])
}

# Performs Box-Cox transformation. The transformation parameter (lambda) is chosen to minimize deviation from the normal distribution (minumum sum of squared skewness and squared kurtosis)
# @param obs The observations to be transformed
# @return The observations transformed
BoxCoxLambdaOnly = function(obs) {
	lambda = seq(0, 1, 0.1)
	obs.trans = matrix(nrow = length(lambda), ncol = length(obs))
	normdev = rep(NA, length(lambda)) # Holds the amount of deviation from the normal distribution
	for(i in 1:length(lambda)) {
		if(lambda[i] == 0) {
			obs.trans[i,] = log(obs)
		} else {
			obs.trans[i,] = (obs^lambda[i]-1)/lambda[i]
		}
		normdev[i] = skewness(obs.trans[i,],na.rm = TRUE)^2 + 0.25*(kurtosis(obs.trans[i,], na.rm = TRUE))^2
	}
	return(lambda[which.min(normdev)])
}

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

# Subtracts the mean
# @param x Data to be centered
# @return Centered data
standardize <- function(x){x-mean(x)}

# Performs F-test between two models
# @param m1 Nested model
# @param m2 Full model
# @return pval The p-value
Ftest <- function(m1,m2){

  RSS1 <- sum(m1$residuals^2)
  RSS2 <- sum(m2$residuals^2)
  m2par <- length(m2$coefficients)
  Nobs <- length(m2$residuals)
  acfs <- acf(m2$residuals)
  Neq <- Nobs/(1+2*acfs[1])

  Fenum <- RSS1 - RSS2
  Fdenom <- RSS2/(nss-m2par)
  Fobs <- Fenum/Fdenom
  pval <- 1 - pf(Fobs,1,Neq-m2par)
  
  return(pval)
}

# Performs Box-Cox transformation for every time trace on a field time object.
# @param cube The field time object
# @return The observations transformed
BoxCoxCubeLambdaOnly = function(data) {
  n.long = dim(data)[1]
  n.lat = dim(data)[2]
  newdata <- array(NA,dim = c(n.long,n.lat))
  for(i in 1:n.long) {
    for(j in 1:n.lat) {
      if(sum(is.na(data[i,j,]))/length(data[i,j,]) < 0.9) {
        dummy <- BoxCoxLambda(data[i,j,])
        newdata[i,j] = dummy$lambda
      }
    }
  }
  return(newdata)
}

# Performs Box-Cox transformation for every time trace on a field time object.
# @param cube The field time object
# @return The observations transformed
BoxCoxCube = function(data) {
	n.long = dim(data)[1]
	n.lat = dim(data)[2]
        n.time = dim(data)[3]
        newdata <- array(NA,dim = c(n.long,n.lat,n.time))
	for(i in 1:n.long) {
		for(j in 1:n.lat) {
			#if(sum(is.na(data[i,j,]))/length(data[i,j,]) < 0.9) {
				newdata[i,j,] = BoxCox(data[i,j,])
			#}
		}
	}
	return(newdata)
}

# Performs Box-Cox transformation for every time trace on a field time object.
# @param cube The field time object
# @return The observations transformed and lambda values
BoxCoxCubeLambda = function(data) {
    n.long <- dim(data)[1]
    n.lat <- dim(data)[2]
    n.time = dim(data)[3]
    lambdas <- array(NA,dim =c(n.long,n.lat))
    newdata <- array(NA,dim = c(n.long,n.lat,n.time))
    for(i in 1:n.long) {
        for(j in 1:n.lat) {
            if(sum(is.na(data[i,j,]))/length(data[i,j,]) < 0.8) {
                out <- BoxCoxLambda(data[i,j,])
                newdata[i,j,] = out$data
                lambdas[i,j] = out$lambda
            }
        }
    }
    return(list(bcData=newdata,lambdas=lambdas))
}


# Latitude part of paralellication of function below
#BoxCoxCubeLambdaLat <- function(i, n.lat, n.times, lambdaKnown = FALSE) {
#    if(lambdaKnown) {
#        for(j in 1:n.lat) {
#            if(sum(is.na(data.to.BC[i,j,])) == 0) {
#                
#        returnObj = bcData
#    } else {
#        datai = list()
#        lambdasi = rep(NA, n.lat)
#        cat("i =", i, "\n")
#        flush.console()
#        for(j in 1:n.lat) {
#            if(sum(is.na(data.to.BC[i,j,])) == 0) {
#                out <- BoxCoxLambda(data.to.BC[i,j,])
#                datai[[j]] = out$data
#                lambdasi[j] = out$lambda
#            } else {
#                datai[[j]] = rep(NA, n.times)
#            }
#        }
#        returnObj = list(bcDdatai = datai, lambdas=lambdasi)
#    }
#    return()
#}    

# Parallel version: Performs Box-Cox transformation for every time trace on a field time object.
# @param cube The field time object
# @return The observations transformed and lambda values
BoxCoxCubeLambdaPar = function(mc.cores = 1, lambdaKnown = FALSE) {
    n.long <- dim(data.to.BC)[1]
    n.lat <- dim(data.to.BC)[2]
    n.times <- dim(data.to.BC)[3]
    lambdas <- array(NA,dim =c(n.long,n.lat))
    i <- 1:n.long
    out <<- mclapply(i, BoxCoxCubeLambdaLat, n.lat = n.lat, n.times = n.times, lambdaKnown = lambdaKnown, mc.cores = mc.cores)
    
    if(lambdaKnown) {
        returnObj = out
    } else {
        bcData = array(dim = c(n.long, n.lat, n.times))
        lambdas = array(dim = c(n.long, n.lat))
        
        # Transform to suitable formats. Any faster way to transform to array object???
        for(i in 1:n.long) {
            for(j in 1:n.lat) {
                bcData[i,j,] = out[[i]][[1]][[j]]
            }
            lambdas[i,] = out[[i]][[2]]
        }
        returnObj = list(bcData = bcData, lambdas = lambdas)
    }
    return(returnObj)
}

# Performs Box-Cox transformation for every time trace on a field time object.
# @param cube The field time object
# @return The observations transformed
BoxCoxCubePar.foreach = function(data,use.max.cores = TRUE,numCores = 6) {
    require(foreach)
    require(doParallel)
    require(parallel)

    #If not number of cores specified, use maximum available cores
    if(use.max.cores){
        numCores <- detectCores()
    } else numCores <- numCores

    cl <- makeCluster(numCores)
    registerDoParallel(cl)
    
    n.long = dim(data)[1]
    n.lat = dim(data)[2]
    n.time = dim(data)[3]
    
    dataTrans <- aperm(data,c(2,1,3))
    #newdata <- array(NA,dim = c(n.lat,n.long,n.time))
    newdata <- foreach(i = 1:3) %dopar%{
        for(j in 1:5) {
            if(sum(is.na(dataTrans[j,i])) == 0) {
                to.newdata <- sum(BoxCox(dataTrans[j,i,]))
            }
        }
        to.newdata
    }
   return(newdata)
}

# Get anomalies
# @param field Find anomalies for this field
# @param baseline The baseline for evaluating the anomalies
# @param n.long Length number of longitude grid cells
# @param n.lat Length number of latitude grid cells
# @param n.time Number of observations in time
# @return The anomalies in every global position for each time
getAnom = function(field,baseline, n.long, n.lat, n.time) {
  anomalies <- array(NA,dim = c(n.long,n.lat,n.time))
  for (i in 1:n.time){
  anomalies[,,i] <- field[,,i]-baseline
}
  return(anomalies)
}

# Get Year Vector
# @param yr.init Start year
# @param n.time Number of time steps
# @param month Month vector
# @return A vector of years
getYear = function(yr.init, n.time,month) {
  yrs <- array(NA,n.time)
  for (i in 1:n.time){
    yrs[i] = yr.init
    if(month[i] == 12) yr.init = yr.init + 1
  }
  return(yrs)
}

# Computes the gradient in each position
# @param field Find gradient in every global position and time for this field
# @param latitude Latitude
# @param longitude Loongitude
# @return The gradient of the SLP in every global position for each time step
getGradients = function(field = newfield, longitude = longitude, latitude = latitude,time = time){
  require(raster)
  require(rasterVis)
  
  n.long <- length(longitude)
  n.lat <- length(latitude)
  n.time <- length(time)
  gradField <- array(NA,dim = c(n.long,n.lat,length(time)))

  for(i in 1:n.time){
    rasterField <- raster(ncol = n.long,nrow = n.lat,xmn = min(longitude),xmx = max(longitude),ymn = min(latitude),ymx = max(latitude))
    rasterField[] <- t(field[,,i])
    slope <- terrain(rasterField,unit = "degrees",neighbours = 4)
    slopeVal <- getValues(slope)
    slopeField <- matrix(slopeVal,nrow = length(longitude),ncol = length(latitude))
    gradField[,,i] <- slopeField
  }
  return(gradField)
}


# Fits regression models betwenn response (SWH) and covariates (SLP, gradient SLP, ...)
# @param SWH SWH in every global position over a given time span (training set)
# @param covariates List of covariates over the same global positions (SLP, gradient SLP, ...)
# @param n.long Length number of longitude grid cells
# @param n.lat Length number of latitude grid cells
# @param n.times Number of observations in time
# @param gradMethod The method used to compute the gradients
# @param model The (regression) model used to fit SWH and covariates
# @return The fitted models
fitModel = function(SWH, covariates, n.long, n.lat, n.times, model, other = NULL) {
    if(model == "WangSwail2006") {
        SLP = covariates$SLP
        SLP.grad = covariates$SLP.grad
        fit = list() # Store the regression model in every global position
        for(j in 1:n.long) {
            fitlat = list()
            for(k in 1:n.lat) {
                if( (sum(is.na(SWH[j,k,])) + sum(is.na(SLP[j,k,])) + sum(is.na(SLP.grad[j,k,]))) == 0) {
                    SLPjk = SLP[j,k,]  # Renaming to be able to call lm.predict. Any better way to do this?
                    gradjk = SLP.grad[j,k,]                    
                    fitlat[[k]] = lm(SWH[j,k,] ~ SLPjk + gradjk)
                }
            }
            fit[[j]] = fitlat
        }
    } else if(model == "Wangetal2012") {
        SLP = covariates$SLP
        SLP.grad = covariates$SLP.grad
        # Performs the analysis separated for four seasons. I.e. we end up with four linear regressions in each global position
        n.season = 4 # JFM, AMJ, JAS, OND
        seasons = list(c(1:3), c(4:6), c(7:9), c(10:12))
        months = rep(1:12, round(n.times/12))
        fit = list() # Store the regression models in every global position
        for(j in 1:n.long) {
            fitlat = list()
            for(k in 1:n.lat) {
                if( (sum(is.na(SWH[j,k,])) + sum(is.na(SLP[j,k,])) + sum(is.na(SLP.grad[j,k,]))) == 0) {
                    fitseason = list()
                    for(i in 1:n.season) {
                        idx.months = months %in% seasons[[i]]
                        SLPjk = SLP[j,k,idx.months]  # Renaming to be able to call lm.predict. Any better way to do this?
                        gradjk = BoxCox(SLP.grad[j,k,idx.months])
                        fitseason[[i]] = lm(BoxCox(SWH[j,k,idx.months]) ~ SLPjk + gradjk)
                    }
                    fitlat[[k]] = fitseason 
                }
            }
            fit[[j]] = fitlat
        }        
    } else {
        ...
    }
    return(fit)
}

# Run the latitude part of the model fit in Wang & Swail (2006)
fitWang06Lat <- function(j) {
    # Use the global versions veriables below??
    fitlat = list()
    for(k in 1:n.lat) {
        if( (sum(is.na(SWH.anom.training[j,k,])) + sum(is.na(SLP.anom.training[j,k,])) + sum(is.na(SLP.grad.training[j,k,]))) == 0) {
            SLPjk = SLP.anom.training[j,k,]  # 1)Renaming to be able to call lm.predict. Any better way to do this?
            gradjk = SLP.grad.training[j,k,]                    
            fitlat[[k]] = lm(SWH.anom.training[j,k,] ~ SLPjk + gradjk)
        } else {
            fitlat[[k]] = NA
        }
    }
    return(fitlat)
}

# Run the latitude part of the model fit in Wang et al. (2012)
fitWang12Lat <- function(j) {
    fitlat = list()
    for(k in 1:n.lat) {
        if( (sum(is.na(SWH.training.BC[j,k,])) + sum(is.na(SLP.anom.training[j,k,])) + sum(is.na(SLP.grad.sq.BC.training[j,k,]))) == 0) {
            fitseason = list()
            for(i in 1:n.seasons.wang12) {
                idx.months = months.wang12 %in% seasons.wang12[[i]]
                SLPjk = SLP.anom.training[j,k,idx.months]  # Renaming to be able to call lm.predict. Any better way to do this?
                gradjk = SLP.grad.sq.BC.training[j,k,idx.months]
                fitseason[[i]] = lm(SWH.training.BC[j,k,idx.months] ~ SLPjk + gradjk)
            }
            fitlat[[k]] = fitseason 
        } else {
            fitlat[[k]] = NA
        }    
    }
    return(fitlat)
}

# Parallel version: Fits regression models betwenn response (SWH) and covariates (SLP, gradient SLP, ...)
# @param SWH SWH in every global position over a given time span (training set)
# @param covariates List of covariates over the same global positions (SLP, gradient SLP, ...)
# @param n.long Length number of longitude grid cells
# @param n.lat Length number of latitude grid cells
# @param n.times Number of observations in time
# @param gradMethod The method used to compute the gradients
# @param model The (regression) model used to fit SWH and covariates
# @return The fitted models
fitModelPar = function(model, mc.cores = 1, other = NULL) {
    if(model == "WangSwail2006") {
        fit = list() # Store the regression model in every global position
        
        j <- 1:n.long
        fit <- mclapply(j, fitWang06Lat, mc.cores = mc.cores)
    } else if(model == "Wangetal2012") {
        # Performs the analysis separated for four seasons.wang12. I.e. we end up with four linear regressions in each global position

        j <- 1:n.long
        fit <- mclapply(j, fitWang12Lat, mc.cores = mc.cores)
    } else {
        ...
    }
    return(fit)
}


# Find area with land
# @param SWH SWH in every global position over a given time span 
# @param covariates List of covariates over the same global positions (SLP, gradient SLP, ...)
# @param n.long Length number of longitude grid cells
# @param n.lat Length number of latitude grid cells
# @param n.times Number of observations in time
# @return Factor matrix indicating areas of land
findLand = function(SWH, covariates, n.long, n.lat) {    
    thisIsLand <- array(NA,dim = c(n.long, n.lat))
    SLP = covariates$SLP
    SLP.grad = covariates$SLP.grad
    
    for(j in 1:n.long) {
        for(k in 1:n.lat) {
            if( (sum(is.na(SWH[j,k,])) + sum(is.na(SLP[j,k,])) + sum(is.na(SLP.grad[j,k,]))) == 0) {
                thisIsLand[j,k] = 0
            } else thisIsLand[j,k] = 1
        }
    }
    return(thisIsLand)
}

# Find proportion of area where both parameters in LM is statistically significant on a specified level
# @param Pmat Matrix of positions where parameters are significant
# @param land Matrix of positions where theer is land
# @return Proportion of sea area where both parameters in LM is statistically significant
findPropSign = function(Pmat,land) {

    #Find total area, need to remove the frame around which comes from the gradient
    totalArea <- dim(land)[1]*dim(land)[2]-2*dim(land)[1]-2*(dim(land)[2]-2)
    indLand <- length(which(land == 1))
    propSign <- length(which(Pmat==1))/(totalArea-indLand)

    return(propSign)
}



# Extract R-squared values of SWH in each global position and time using the fitted model from the function fitModel
# @param fit Model fit
# @param n.long Length number of longitude grid cells
# @param n.lat Length number of latitude grid cells
extractRSquared = function(fit,n.long,n.lat) {
  rSquared <- array(NA,dim = c(n.long, n.lat))
  for(j in 1:n.long) {
    if(length(fit[[j]])!=0){
      n.lat2 <- length(fit[[j]])
    
      for(k in 1:n.lat2) {
        if(!is.null(fit[[j]][[k]]))
        rSquared[j,k] <- summary(fit[[j]][[k]])$r.squared
      }
    }
  }
  return(rSquared)
}

# Extract estimated coefficients, SD and p-values values of SWH in each global position and time using the fitted model from the function fitModel
# @param fit Model fit
# @param n.long Length number of longitude grid cells
# @param n.lat Length number of latitude grid cells
extractCoeffsLinFit = function(fit,n.long,n.lat) {
    amat <- array(NA,dim = c(n.long, n.lat))
    bmat <- array(NA,dim = c(n.long, n.lat))
    amatSD <- array(NA,dim = c(n.long, n.lat))
    bmatSD <- array(NA,dim = c(n.long, n.lat))
    amatP <- array(NA,dim = c(n.long, n.lat))
    bmatP <- array(NA,dim = c(n.long, n.lat))
    Pmat80 <- array(NA,dim = c(n.long, n.lat))
    Pmat95 <- array(NA,dim = c(n.long, n.lat))
    Pmat9580 <- array(NA,dim = c(n.long, n.lat))

    for(j in 1:n.long) {
        if(length(fit[[j]])!=0){
            n.lat2 <- length(fit[[j]])
            
            for(k in 1:n.lat2) {
                if(!is.null(fit[[j]][[k]])){
                    amat[j,k] <- summary(fit[[j]][[k]])$coefficients[2,1]
                    amatSD[j,k] <- summary(fit[[j]][[k]])$coefficients[2,2]
                    amatP[j,k] <- summary(fit[[j]][[k]])$coefficients[2,4]
                    bmat[j,k] <- summary(fit[[j]][[k]])$coefficients[3,1]
                    bmatSD[j,k] <- summary(fit[[j]][[k]])$coefficients[3,2]
                    bmatP[j,k] <- summary(fit[[j]][[k]])$coefficients[3,4]
                }
            }
        }
    }
    ind80 <- which(((1-amatP)>0.80) & ((1-bmatP)>0.80))
    ind95 <- which(((1-amatP)>0.95) & ((1-bmatP)>0.95))
    Pmat80[ind80] <- 1
    Pmat95[ind90] <- 1
    Pmat9580[ind80] <- 1
    Pmat9580[ind95] <- 2
    return(list(amat = amat,amatSD = amatSD,amatP = amatP,bmat = bmat,bmatSD = bmatSD,bmatP = bmatP,Pmat80 = Pmat80,Pmat90 = Pmat90,Pmat9580 = Pmat9580))
}

# Predicts values of SWH in each global position and time using the fitted model from the function fitModel
# @param covariates0 List of covariates over the same global positions (SLP, gradient SLP, ...)
# @param n.long Length number of longitude grid cells
# @param n.lat Length number of latitude grid cells
# @param n.times0 Number of observations in time
# @param gradMethod The method used to compute the gradients
# @param model The (regression) model used to fit SWH and covariates
predictSWH = function(fit, covariates0, n.long, n.lat, n.times0, model) {
    predictions = array(dim = c(n.long, n.lat, n.times0))
    if(model == "WangSwail2006") {
        SLP = covariates0$SLP
        SLP.grad = covariates0$SLP.grad
        for(j in 1:n.long) {
            for(k in 1:n.lat) {
                if( (sum(is.na(SWH[j,k,])) + sum(is.na(SLP[j,k,])) + sum(is.na(SLP.grad[j,k,]))) == 0) {
                    predictions[j,k,] = predict(object = fit[[j]][[k]], newdata = data.frame(SLPjk = SLP[j,k,], gradjk = SLP.grad[j,k,])) # use the prediction function for lm in R.
                }
            }
        }
    } else if(model == "Wangetal2012") {
        SLP = covariates0$SLP
        SLP.grad = covariates0$SLP.grad
        n.seasons = 4 # JFM, AMJ, JAS, OND
        seasons = list(c(1:3), c(4:6), c(7:9), c(10:12))
        months = rep(1:12, round(n.times0/12))
        for(j in 1:n.long) {
            for(k in 1:n.lat) {
                if( (sum(is.na(SWH[j,k,])) + sum(is.na(SLP[j,k,])) + sum(is.na(SLP.grad[j,k,]))) == 0) {
                    for(i in 1:n.seasons) {
                        idx.months = months %in% seasons[[i]]
                        predictions[j,k,idx.months] = predict(object = fit[[j]][[k]][[i]], newdata = data.frame(SLPjk = SLP[j,k,idx.months], gradjk = BoxCox(SLP.grad[j,k,idx.months]))) # use the prediction function for lm in R.
                    }
                }
            }
        }    
    } else {
        ...    
    }
    return(predictions)
}


# Latitude part of prediction of Wang06
predictWang06Lat <- function(j) {
    predictionsj = list()
    for(k in 1:n.lat) {
        #cat("k =", k, "\n")
        if( (sum(is.na(SLP.anom.test[j,k,])) + sum(is.na(SLP.grad.test[j,k,]))) + sum(is.na(fit[[j]][[k]])) == 0) {
            #print(fit[[j]][[k]])
            #print(SLP.anom.test[j,k,])
            #print(SLP.grad.test[j,k,])
            predictionsj[[k]] = predict(object = fit[[j]][[k]], newdata = data.frame(SLPjk = SLP.anom.test[j,k,], gradjk = SLP.grad.test[j,k,])) # use the prediction function for lm in R.
        } else {
            predictionsj[[k]] = rep(NA, n.times.test)
        }
    }
    return(predictionsj)
}


# Latitude part of prediction of Wang12
predictWang12Lat <- function(j) {
    predictionsj = list()
    for(k in 1:n.lat) {
        predictionsjseas = rep(NA, n.times.test)
        if( (sum(is.na(SLP.anom.test[j,k,])) + sum(is.na(SLP.grad.sq.BC.test[j,k,]))) + sum(is.na(fit[[j]][[k]])) == 0) {
            for(i in 1:n.seasons.wang12) {
                idx.months = months.wang12 %in% seasons.wang12[[i]]
                predictionsjseas[idx.months] = predict(object = fit[[j]][[k]][[i]], newdata = data.frame(SLPjk = SLP.anom.test[j,k,idx.months], gradjk = SLP.grad.sq.BC.test[j,k,idx.months])) # use the prediction function for lm in R.
            }
            predictionsj[[k]] = predictionsjseas
        } else {
            predictionsj[[k]] = rep(NA, n.times.test)
        }
    }    
    return(predictionsj)
}

# Parallel version: Predicts values of SWH in each global position and time using the fitted model from the function fitModel
# @param covariates0 List of covariates over the same global positions (SLP, gradient SLP, ...)
# @param n.long Length number of longitude grid cells
# @param n.lat Length number of latitude grid cells
# @param n.times0 Number of observations in time
# @param gradMethod The method used to compute the gradients
# @param model The (regression) model used to fit SWH and covariates
predictSWHpar = function(model, mc.cores = 1) {
    if(model == "WangSwail2006") {
        j = 1:n.long
        predictions <- mclapply(j, predictWang06Lat, mc.cores = mc.cores)
    } else if(model == "Wangetal2012") {
        j = 1:n.long
        predictions <- mclapply(j, predictWang12Lat, mc.cores = mc.cores)
    } else {
        ...    
    }
    return(predictions)
}


# Computes the rank of a model based on the fit (e.g. R^2) or the predictions from the model
rankModel = function(response, fit, prediction, model, rank.method, control = NULL) {
    if(rank.method == "RMSE") {
        rank = sqrt(mean((response - prediction)^2))
    } else if(rank.method == "relRMSE") {
        rank <- sqrt(mean((response - prediction)^2))/control$baselineSWH
    } else if(rank.method == "MAE") {
        if(control$distr == "boxcox") {
            rank <- abs(qboxcox(0.5, mean = prediction, sd = control$se, lambda = control$lambda) - response)
        }
    } else if(rank.method == "Rsquared") {
        rank <- summary(fit)$r.squared
    } else if(rank.method == "LOGS") {
        if(control$distr == "normal") {
            rank <- logs(y = response, family = "normal", mean = prediction, sd = control$se)
        } else if(control$distr == "boxcox") {
            rank <- dboxcox(x = response, mean = prediction, sd = control$se, lambda = control$lambda, log = TRUE)
        } else {
            ...
        }
    } else if(rank.method == "CRPS") {
        rank <- mean(crps(y = response, family = "normal", mean = prediction, sd = control$se))
    } else if(rank.method == "PSS") { # See the appendix in Lin and Wang (2011) 
        qPSS = control$qPSS # The probability related to the quantile (Q). Compute Q from qPSS based on the real observations.
        #cat("qPSS =", qPSS, "\n")
        Q = quantile(response, probs = qPSS, na.rm = TRUE)
        # Compute contingency table 
        a = sum( (response>=Q) & (prediction>=Q) )
        b = sum( (response<Q) & (prediction>=Q) )
        c = sum( (response>=Q) & (prediction<Q) )
        d = sum( (response<Q) & (prediction<Q) )
        rank = a/(a+c) - b/(b+d) # The PSS score
    } else if(rank.method == "FBI") { # See the appendix in Lin and Wang (2011) 
        qFBI = control$qFBI # The probability related to the quantile (Q). Compute Q from qFBI based on the real observations.
        Q = quantile(response, probs = qFBI, na.rm = TRUE)
        # Compute contingency table based on the quantile Q
        a = sum( (response>=Q) & (prediction>=Q) )
        b = sum( (response<Q) & (prediction>=Q) )
        c = sum( (response>=Q) & (prediction<Q) )
        d = sum( (response<Q) & (prediction<Q) )
        rank = (a + b)/(a + c) # The FBI score
    } else {
        ...
    }
    return(rank)
}

rankModelLat <- function(j, rank.method, control = NULL) {
    rankj = rep(NA, n.lat)
    if(rank.method == "RMSE") {
        for(k in 1:n.lat) {
            #cat("j =", j, "k =", k, "\n")
            if( (sum(is.na(SWH.test[j,k,])) + sum(is.na(predictions[[j]][[k]]))) == 0) {
                rankj[k] = sqrt(mean((SWH.test[j,k,] - predictions[[j]][[k]])^2))
            }
        }
    } else if(rank.method == "Rsquared") {
        for(k in 1:n.lat) {
            if(model == "WangSwail2006") {
                if( (sum(is.na(SWH.test[j,k,])) + sum(is.na(predictions[[j]][[k]]))) == 0) {
                    rankj[k] <- summary(fit[[j]][[k]])$r.squared
                } 
            } else if(model == "Wangetal2012") {
                if( (sum(is.na(SWH.test[j,k,])) + sum(is.na(predictions[[j]][[k]]))) == 0) {
                    rank[k] = 0
                    for(i in 1:n.seasons.wang12) {
                        rankj[k] <- summary(fit[[j]][[k]][[i]])$r.squared/n.seasons.wang12
                    }
                }                
            }
        }
    } else if(rank.method == "PSS") { # See the appendix in Lin and Wang (2011) 
        qPSS = control$qPSS # The probability related to the quantile (Q). Compute Q from qPSS based on the real observations.
        for(k in 1:n.lat) {
            if( (sum(is.na(SWH.test[j,k,])) + sum(is.na(predictions[[j]][[k]]))) == 0) {
                Q = quantile(SWH.test[j,k,], probs = qPSS)
                # Compute contingency table 
                a = sum( (SWH.test[j,k,]>=Q) & (predictions[[j]][[k]]>=Q) )
                b = sum( (SWH.test[j,k,]<Q) & (predictions[[j]][[k]]>=Q) )
                c = sum( (SWH.test[j,k,]>=Q) & (predictions[[j]][[k]]<Q) )
                d = sum( (SWH.test[j,k,]<Q) & (predictions[[j]][[k]]<Q) )
                rankj[k] = a/(a+c) - b/(b+d) # The PSS score
            }
        }
    } else if(rank.method == "FBI") { # See the appendix in Lin and Wang (2011) 
        qFBI = control$qFBI # The probability related to the quantile (Q). Compute Q from qFBI based on the real observations.
        for(k in 1:n.lat) {
            if( (sum(is.na(SWH.test[j,k,])) + sum(is.na(predictions[[j]][[k]]))) == 0) {
                Q = quantile(SWH.test[j,k,], probs = qFBI)
                # Compute contingency table based on the quantile Q
                a = sum( (SWH.test[j,k,]>=Q) & (predictions[[j]][[k]]>=Q) )
                b = sum( (SWH.test[j,k,]<Q) & (predictions[[j]][[k]]>=Q) )
                c = sum( (SWH.test[j,k,]>=Q) & (predictions[[j]][[k]]<Q) )
                d = sum( (SWH.test[j,k,]<Q) & (predictions[[j]][[k]]<Q) )
                rankj[k] = (a + b)/(a + c) # The FBI score
            }
        }
    } else {
        ...
    }
    return(rankj)
}

# Paralell version: Computes the rank of a model based on the fit (e.g. R^2) or the predictions from the model
# @param SWH SWH in every global position over a given time span (training set). For Wang et al. (2012) it is the BoxCox transform of SWH.
# @param predictions The predictions from function predictSWH
# @param fit The fitted model from function fitModel
# @param n.long Length number of longitude grid cells
# @param n.lat Length number of latitude grid cells
# @param rank.method Method to compute difference to the observations
# @control List of other variables to compute the ranks
rankModelPar <- function(rank.method, control = NULL, mc.cores = 1) {
    j = 1:n.long
    rank <- mclapply(j, rankModelLat, rank.method = rank.method, control = control, mc.cores = mc.cores)
    rank <- matrix(unlist(rank), nrow = n.long, byrow = TRUE)
    return(rank)
}


runParallel <- function(j, model, part) {
    #cat("j =", j, "\n")
    #flush.console()
    rankj = array(dim = c(n.lat, length(rank.methods)))
    for(k in 1:n.lat) {
        if( (sum(is.na(SWH[j,k,])) + sum(is.na(SLP[j,k,])) + sum(is.na(SLP.grad[j,k,]))) == 0) {
            if(model == "WangSwail2006") { # Use anomalies
                # fit model
                SLP.anom.pred <- SLP.anom[j,k,idx.training]  # Renaming to be able to call lm.predict. Any better way to do this?
                SLP.grad.pred <- SLP.grad[j,k,idx.training]                    
                fit <- lm(SWH.anom[j,k,idx.training] ~ SLP.anom.pred + SLP.grad.pred)
                
                # predict
                prediction <- predict(object = fit, newdata = data.frame(SLP.anom.pred = SLP.anom[j,k,idx.test], SLP.grad.pred = SLP.grad[j,k,idx.test]))
            
                # rank
                for(l in 1:length(rank.methods)) {
                    rankj[k,l] = rankModel(response = SWH.anom[j,k,idx.test], fit, prediction, model, rank.methods[l], control)
                }
            } else if(model == "Wangetal2012") { # Do not use anomalies (???)
                rankj[k,] = 0                
                for(i in 1:n.seasons.wang12) {
                    # fit model
                    SWH.BC.training <- BoxCoxLambda(SWH[j,k, idx.training & idx.months[[i]] ])
                    SWH.BC.lambda <- SWH.BC.training$lambda
                    SWH.BC.training <- SWH.BC.training$data
                    
                    SLP.anom.pred <- SLP.anom[j,k, idx.training & idx.months[[i]] ]  # Renaming to be able to call lm.predict. Any better way to do this?

                    SLP.grad.sq.BC.training <- BoxCoxLambda(SLP.grad[j,k, idx.training & idx.months[[i]] ]^2)
                    SLP.grad.sq.lambda <- SLP.grad.sq.BC.training$lambda
                    SLP.grad.sq.BC.training <- SLP.grad.sq.BC.training$data

                    # Add covariates related to PC

                    fit <- lm(SWH.BC.training ~ SLP.anom.pred + SLP.grad.sq.BC.training)
                
                    # predict
                    SLP.anom.pred <- SLP.anom[j,k, idx.test & idx.months[[i]] ]  # Renaming to be able to call lm.predict. Any better way to do this?
                    
                    SLP.grad.sq.BC.test <- BoxCoxLambda(SLP.grad[j,k, idx.test & idx.months[[i]] ]^2)
                    SLP.grad.sq.lambda <- SLP.grad.sq.BC.test$lambda
                    SLP.grad.sq.BC.test <- SLP.grad.sq.BC.test$data
                    # Rescale the test data
                    SLP.grad.sq.BC.test <- sd(SLP.grad.sq.BC.training)/sd(SLP.grad.sq.BC.test)*(SLP.grad.sq.BC.test - mean(SLP.grad.sq.BC.test)) + mean(SLP.grad.sq.BC.training)
                    
                    prediction <- predict(object = fit, newdata = data.frame(SLP.anom.pred = SLP.anom.pred, SLP.grad.sq.BC.training = SLP.grad.sq.BC.test))
                
                    # rank
                    SWH.BC.test <- BoxCoxLambda(SWH[j,k, idx.test & idx.months[[i]] ])
                    SWH.BC.lambda <- SWH.BC.test$lambda
                    SWH.BC.test <- SWH.BC.test$data
                    # Rescale the test data
                    SWH.BC.test <- sd(SWH.BC.training)/sd(SWH.BC.test)*(SWH.BC.test - mean(SWH.BC.test)) + mean(SWH.BC.training)
                   
                    for(l in 1:length(rank.methods)) {
                        rankj[k,l] = rankj[k,l] + rankModel(response = SWH.BC.test, fit, prediction, model, rank.methods[l], control)
                    }
                }
            }
        }
    }
    return(rankj)
}

runParallelERAINTERIM <- function(j, model, part) {
    #cat("j =", j, "\n")
    flush.console()
    rankj = array(dim = c(n.latSWH, length(rank.methods)))
    for(k in 1:n.latSWH) {
        #cat("k = ", k, "\n")
        if(model == "WangSwail2006") { # !OBS! The Wang06 code below has not been updated for ERA INTERIM data
            # fit model
            SLP.anom.pred <- SLP.anom[j,k,idx.training]  # Renaming to be able to call lm.predict. Any better way to do this?
            SLP.grad.pred <- SLP.grad[j,k,idx.training]                    
            fit <- lm(SWH.anom[j,k,idx.training] ~ SLP.anom.pred + SLP.grad.pred)
                
            # predict
            prediction <- predict(object = fit, newdata = data.frame(SLP.anom.pred = SLP.anom[j,k,idx.test], SLP.grad.pred = SLP.grad[j,k,idx.test]))
            
            # rank
            for(l in 1:length(rank.methods)) {
                rankj[k,l] = rankModel(response = SWH.anom[j,k,idx.test], fit, prediction, model, rank.methods[l], control)
            }
          } else if(model == "Wangetal2012") { # Do not use anomalies (???)
            rankj[k,] = 0                
            
            idx.longSLP = (1:n.longSLP)[ longitudeSLP == longitudeSWH[j] ]
            idx.latSLP = (1:n.latSLP)[ latitudeSLP == latitudeSWH[k] ]
            # fit model

            idx.ar = idx.training-1
            idx.ar[idx.ar < 1] = 1
            predictors.training = as.data.frame(cbind(SLP.anom.scaled[idx.longSLP, idx.latSLP, idx.training],
              SLP.grad.bc.anom.scaled[idx.longSLP, idx.latSLP, idx.training], SLP.PC.scaled.training, SLP.grad.PC.scaled.training,
              SWH.bc.scaled[j, k, idx.ar]))
            idx.ar = idx.test-1
            idx.ar[idx.ar < 1] = 1
            predictors.test = as.data.frame(cbind( SLP.anom.scaled[idx.longSLP, idx.latSLP, idx.test],
              SLP.grad.bc.anom.scaled[idx.longSLP, idx.latSLP, idx.test], SLP.PC.scaled.test, SLP.grad.PC.scaled.test,
                                                  SWH.bc.scaled[j, k, idx.ar]))

            
            if( sum(is.na(SWH.bc.scaled[j, k, ])) + sum(is.na(predictors.test)) <= 10 & sum(is.na(predictors.training)) <= 5) {
                #fit <- lm(SWH.bc.scaled[j, k, idx.training] ~ ., data = predictors.training) # Use arima(...) later
                fit <- arima(SWH.bc.scaled[j, k, idx.training], order = c(1,0,0), xreg = predictors.training) # AR(1) noise
                
                # predict
                #SWH.bc.scaled.pred <- predict(object = fit, newdata = predictors.test) #lm
                SWH.bc.scaled.pred <- predict(fit, newxreg = predictors.test)$pred #arima
            
                # Descale and detransform before computing different raknkings
                SWH.bc.pred <- SWH.bc.scaled.pred*SWH.bc.sd[j, k] # descale
                SWH.pred = InvBoxCox(SWH.bc.pred, SWH.bc.lambdas[j, k]) #detransform
                    
                SWH.bc.test = SWH.bc.scaled[j, k, idx.test]*SWH.bc.sd[j, k] #descale
                SWH.test = InvBoxCox(SWH.bc.test, SWH.bc.lambdas[j, k]) #detransfrom
                
                #for(l in 1:length(rank.methods)) {
                #    rankj[k,l] = rankj[k,l] + rankModel(response = SWH.test, fit, SWH.pred, model, rank.methods[l], control)
                #}
                
                # Handling NA's. Occurs rarely.
                not.NA <- !is.na(SWH.test) & !is.na(SWH.pred)
                if(sum(!not.NA ) > 0) {
                    cat("j =", j, "k =", k, "sum(is.na(SWH.test)) =", sum(is.na(SWH.test)), "sum(is.na(SWH.pred)) =", sum(is.na(SWH.pred)), "\n")
                }
                
                rankj[k,1] = rankj[k,1] + rankModel(response = SWH.test[not.NA], fit, SWH.pred[not.NA], model, "RMSE")                
                rankj[k,2] = rankj[k,2] + rankModel(response = SWH.test[not.NA], fit, SWH.pred[not.NA], model, "relRMSE", list(baselineSWH = SWH.baseline[j,k]))                
                rankj[k,3] = rankj[k,3] + rankModel(response = SWH.test[not.NA], fit, SWH.pred[not.NA], model, "PSS", list(qPSS = 0.5))
                rankj[k,4] = rankj[k,4] + rankModel(response = SWH.test[not.NA], fit, SWH.pred[not.NA], model, "PSS", list(qPSS = 0.95))
                rankj[k,5] = rankj[k,5] + rankModel(response = SWH.test[not.NA], fit, SWH.pred[not.NA], model, "PSS", list(qPSS = 0.975))
                rankj[k,6] = rankj[k,6] + rankModel(response = SWH.test[not.NA], fit, SWH.pred[not.NA], model, "FBI", list(qFBI = 0.5))
                rankj[k,7] = rankj[k,7] + rankModel(response = SWH.test[not.NA], fit, SWH.pred[not.NA], model, "FBI", list(qFBI = 0.95))
                rankj[k,8] = rankj[k,8] + rankModel(response = SWH.test[not.NA], fit, SWH.pred[not.NA], model, "FBI", list(qFBI = 0.975))
                rankj[k,9] = rankj[k,9] + fit$coef[1] #AR1
                rankj[k,10] = rankj[k,10] + fit$coef[2] #intercept
                rankj[k,11] = rankj[k,11] + fit$coef[3] #SLP
                rankj[k,12] = rankj[k,12] + fit$coef[4]

            } else {
                rankj[k,] = NA
                #cat("j = ", j, "k =", k, "sum(is.na(SWH.bc.scaled[j, k, ])) =", sum(is.na(SWH.bc.scaled[j, k, ])), "sum(is.na(predictors.test))", sum(is.na(predictors.test)), "sum(is.na(predictors.training))", sum(is.na(predictors.training)), "\n")
                #print(SWH.bc.scaled[j, k, ])
                #print(predictors.test)
                #print(predictors.training)
                #readline()
            }
        }
    }
    #print(rankj)
    #cat("\n")
    return(rankj)
}

studyCoeffParallelERAINTERIM <- function(j) {
    coeffj = array(dim = c(n.latSWH, n.coeff, 4)) # both estimates, standard errors t-val and p-val
    for(k in 1:n.latSWH) {
        #cat("k = ", k, "\n")
        
        idx.longSLP = (1:n.longSLP)[ longitudeSLP == longitudeSWH[j] ]
        idx.latSLP = (1:n.latSLP)[ latitudeSLP == latitudeSWH[k] ]

        # fit model
        idx.ar = idx-1
        idx.ar[idx.ar < 1] = 1
        predictors = as.data.frame(cbind(SLP.anom.scaled[idx.longSLP, idx.latSLP, idx], SLP.grad.bc.anom.scaled[idx.longSLP, idx.latSLP, idx], SWH.bc.scaled[j, k, idx.ar]))
        #predictors = as.data.frame(cbind(SLP.anom.scaled[idx.longSLP, idx.latSLP], SLP.grad.bc.anom.scaled[idx.longSLP, idx.latSLP, idx], SLP.PC.scaled.training, SLP.grad.PC.scaled.training, SWH.bc.scaled[j, k, idx.ar]))

        if( sum(is.na(SWH.bc.scaled[j, k, ])) + sum(is.na(predictors)) <= 10) {
            fit <- summary(lm(SWH.bc.scaled[j, k, idx] ~ ., data = predictors)) # Use arima(...) later)
            coeffj[k,,] = fit$coefficients
            #fit <- arima(SWH.bc.scaled[j, k, idx.training], order = c(1,0,0), xreg = predictors) # AR(1) noise
        }
    }
    return(coeffj)
}

fitmodelERAINTERIM <- function(j, model, part) {
    
    flush.console()
    rankj = array(dim = c(n.latSWH, length(rank.methods)))
    for(k in 1:n.latSWH) {
        rankj[k,] = 0                
        
        idx.longSLP = (1:n.longSLP)[ longitudeSLP == longitudeSWH[j] ]
        idx.latSLP = (1:n.latSLP)[ latitudeSLP == latitudeSWH[k] ]

        ## fit model

        idx.ar = idx.training-1
        idx.ar[idx.ar < 1] = 1
        predictors.training = as.data.frame(cbind(SLP.anom.scaled[idx.longSLP, idx.latSLP, idx.training],
                                                  SLP.grad.bc.anom.scaled[idx.longSLP, idx.latSLP, idx.training], SLP.PC.scaled.training, SLP.grad.PC.scaled.training,
                                                  SWH.bc.scaled[j, k, idx.ar]))
        idx.ar = idx.test-1
        idx.ar[idx.ar < 1] = 1
        predictors.test = as.data.frame(cbind( SLP.anom.scaled[idx.longSLP, idx.latSLP, idx.test],
                                              SLP.grad.bc.anom.scaled[idx.longSLP, idx.latSLP, idx.test], SLP.PC.scaled.test, SLP.grad.PC.scaled.test,
                                              SWH.bc.scaled[j, k, idx.ar]))

        
        if( sum(is.na(SWH.bc.scaled[j, k, ])) + sum(is.na(predictors.test)) <= 10 & sum(is.na(predictors.training)) <= 5) {
                                        #fit <- lm(SWH.bc.scaled[j, k, idx.training] ~ ., data = predictors.training) # Use arima(...) later
            if(k == 1) {
                init = NULL
            } else {
                init = fit$coef
            }
            fit <- arima(SWH.bc.scaled[j, k, idx.training], order = c(1,0,0), xreg = predictors.training, init = init) # AR(1) noise
            
            ## predict
            ##SWH.bc.scaled.pred <- predict(object = fit, newdata = predictors.test) #lm
            SWH.bc.scaled.pred <- predict(fit, newxreg = predictors.test)$pred #arima
            
            ## Descale and detransform before computing different raknkings
            SWH.bc.pred <- SWH.bc.scaled.pred*SWH.bc.sd[j, k] # descale
            SWH.pred = InvBoxCox(SWH.bc.pred, SWH.bc.lambdas[j, k]) #detransform
            
            SWH.bc.test = SWH.bc.scaled[j, k, idx.test]*SWH.bc.sd[j, k] #descale
            SWH.test = InvBoxCox(SWH.bc.test, SWH.bc.lambdas[j, k]) #detransfrom
            
            ##for(l in 1:length(rank.methods)) {
            ##    rankj[k,l] = rankj[k,l] + rankModel(response = SWH.test, fit, SWH.pred, model, rank.methods[l], control)
            ##}
            
            ## Handling NA's. Occurs rarely.
            not.NA <- !is.na(SWH.test) & !is.na(SWH.pred)
            if(sum(!not.NA ) > 0) {
                cat("j =", j, "k =", k, "sum(is.na(SWH.test)) =", sum(is.na(SWH.test.jfm)), "sum(is.na(SWH.pred)) =", sum(is.na(SWH.pred)), "\n")
            }
            
            rankj[k,1] = rankj[k,1] + rankModel(response = SWH.test[not.NA], fit, SWH.pred[not.NA], model, "RMSE")                
            rankj[k,2] = rankj[k,2] + rankModel(response = SWH.test[not.NA], fit, SWH.pred[not.NA], model, "relRMSE", list(baselineSWH = SWH.baseline[j,k]))                
            rankj[k,3] = rankj[k,3] + rankModel(response = SWH.test[not.NA], fit, SWH.pred[not.NA], model, "PSS", list(qPSS = 0.5))
            rankj[k,4] = rankj[k,4] + rankModel(response = SWH.test[not.NA], fit, SWH.pred[not.NA], model, "PSS", list(qPSS = 0.95))
            rankj[k,5] = rankj[k,5] + rankModel(response = SWH.test[not.NA], fit, SWH.pred[not.NA], model, "PSS", list(qPSS = 0.975))
            rankj[k,6] = rankj[k,6] + rankModel(response = SWH.test[not.NA], fit, SWH.pred[not.NA], model, "FBI", list(qFBI = 0.5))
            rankj[k,7] = rankj[k,7] + rankModel(response = SWH.test[not.NA], fit, SWH.pred[not.NA], model, "FBI", list(qFBI = 0.95))
            rankj[k,8] = rankj[k,8] + rankModel(response = SWH.test[not.NA], fit, SWH.pred[not.NA], model, "FBI", list(qFBI = 0.975))
            rankj[k,9] = rankj[k,9] + fit$coef[1] #AR1
            rankj[k,10] = rankj[k,10] + fit$coef[2] #intercept
            rankj[k,11] = rankj[k,11] + fit$coef[3] #SLP
            rankj[k,12] = rankj[k,12] + fit$coef[4]

        } else {
            rankj[k,] = NA
        }
    }
    
    return(rankj)
}








runParallelWangVanemWalker.4seasons <- function(j) {
    rankj = array(dim = c(length(area.latSWH), n.models.4seasons, n.training.test, n.rank.methods) )
    rankj.all = array(dim = c(length(area.latSWH), n.models.4seasons, n.training.test, 4, n.time.max) )
        
    for(p in 1:length(training.test)) {

        ## Divide in training and test set
        idx.training = training.test[[p]][[1]]
        idx.test = training.test[[p]][[2]]
    
        for(k in 1:length(area.latSWH)) {            
            idx.longSLP = which(longitudeSLP == longitudeSWH[j])
            idx.latSLP = which(latitudeSLP == latitudeSWH[area.latSWH[k]])

            if(sum(is.na(SWH.bc.scaled[j, area.latSWH[k],])) < na.thresh & sum(is.na( SLP.anom.scaled[idx.longSLP, idx.latSLP, ])) < na.thresh & sum(is.na( SLP.grad.bc.anom.scaled[idx.longSLP, idx.latSLP, ])) < na.thresh & !is.na(SWH.bc.sd[j,area.latSWH[k]]) & !is.na(SLP.grad.bc.baseline[j,area.latSWH[k]]) & !is.na(SLP.grad.bc.anom.sd[j,area.latSWH[k]])) {
                
                not.NA = which(!is.na(SWH.bc.scaled[j, area.latSWH[k],]) & !is.na( SLP.anom.scaled[idx.longSLP, idx.latSLP, ]) & !is.na( SLP.grad.bc.anom.scaled[idx.longSLP, idx.latSLP, ]))
                idx.training = idx.training[idx.training %in% not.NA]
                idx.test = idx.test[idx.test %in% not.NA]

                SLP.PC.scaled.training <- SLP.PC.scaled[idx.training, ]
                SLP.PC.scaled.test <- SLP.PC.scaled[idx.test, ]

                SLP.grad.PC.scaled.training <- SLP.grad.PC.scaled[idx.training, ]
                SLP.grad.PC.scaled.test <- SLP.grad.PC.scaled[idx.test, ]
                
                if(change.lambda) {                    
                    # Descale and InvBoxCox to get original SWH
                    SWH.org.training = InvBoxCox(SWH.bc.scaled[j, area.latSWH[k], idx.training]*SWH.bc.sd[j,area.latSWH[k]], SWH.bc.lambdas[j,area.latSWH[k]])
                    SWH.org.test = InvBoxCox(SWH.bc.scaled[j, area.latSWH[k], idx.test]*SWH.bc.sd[j,area.latSWH[k]], SWH.bc.lambdas[j,area.latSWH[k]])
                    # Perform BoxCox with another lambda value
                    SWH.bc.dummy.training = BoxCoxLambda(SWH.org.training)
                    SWH.bc.lambda.training = SWH.bc.dummy.training$lambda
                    SWH.bc.training = SWH.bc.dummy.training$data
                    SWH.bc.test = BoxCoxLambdaKnown2(SWH.org.test, SWH.bc.lambda.training)
                    # Upgrade values
                    SWH.bc.scaled[j, area.latSWH[k], idx.training] = SWH.bc.training/SWH.bc.sd[j,area.latSWH[k]]
                    SWH.bc.scaled[j, area.latSWH[k], idx.test] = SWH.bc.test/SWH.bc.sd[j,area.latSWH[k]]
                    SWH.bc.lambdas[j,area.latSWH[k]] = SWH.bc.lambda.training
                    
                    # Descale and InvBoxCox to get original SLP.grad
                    SLP.grad.org.training = InvBoxCox(SLP.grad.bc.anom.scaled[j,area.latSWH[k],idx.training]*SLP.grad.bc.anom.sd[j,area.latSWH[k]] + SLP.grad.bc.baseline[j,area.latSWH[k]], SLP.grad.bc.lambdas[j,area.latSWH[k]])
                    SLP.grad.org.test = InvBoxCox(SLP.grad.bc.anom.scaled[j,area.latSWH[k],idx.test]*SLP.grad.bc.anom.sd[j,area.latSWH[k]] + SLP.grad.bc.baseline[j,area.latSWH[k]], SLP.grad.bc.lambdas[j,area.latSWH[k]])
                    # Perform BoxCox with another lambda value
                    SLP.grad.bc.dummy.training = BoxCoxLambda(SLP.grad.org.training)
                    SLP.grad.bc.lambda.training = SLP.grad.bc.dummy.training$lambda
                    SLP.grad.bc.training = SLP.grad.bc.dummy.training$data
                    SLP.grad.bc.test = BoxCoxLambdaKnown2(SLP.grad.org.test, SLP.grad.bc.lambda.training)
                    # Upgrade values
                    SLP.grad.bc.anom.scaled[j,area.latSWH[k],idx.training] = (SLP.grad.bc.training - SLP.grad.bc.baseline[j,area.latSWH[k]])/SLP.grad.bc.anom.sd[j,area.latSWH[k]]
                    SLP.grad.bc.anom.scaled[j,area.latSWH[k],idx.test] = (SLP.grad.bc.test - SLP.grad.bc.baseline[j,area.latSWH[k]])/SLP.grad.bc.anom.sd[j,area.latSWH[k]]
                    SLP.grad.bc.lambdas[j,area.latSWH[k]] = SLP.grad.bc.lambda.training
                }
                
                ##### Wang et al. 2012 LASSO#####
                cat("4seasons, j =", j, "k =", area.latSWH[k], "\n")
                # fit model            
                idx.ar1 = c(rep(idx.training[1],1), idx.training[1:(length(idx.training)-1)])
                idx.ar2 = c(rep(idx.training[1],2), idx.training[1:(length(idx.training)-2)])
                predictors.training = as.data.frame(cbind(SLP.anom.scaled[idx.longSLP, idx.latSLP, idx.training], SLP.grad.bc.anom.scaled[idx.longSLP, idx.latSLP, idx.training], SLP.PC.scaled.training, SLP.grad.PC.scaled.training, SWH.bc.scaled[j, area.latSWH[k], idx.ar1], SWH.bc.scaled[j, area.latSWH[k], idx.ar2]))
                idx.ar1 = c(rep(idx.test[1],1), idx.test[1:(length(idx.test)-1)])
                idx.ar2 = c(rep(idx.test[1],2), idx.test[1:(length(idx.test)-2)])                
                predictors.test = as.data.frame(cbind( SLP.anom.scaled[idx.longSLP, idx.latSLP, idx.test], SLP.grad.bc.anom.scaled[idx.longSLP, idx.latSLP, idx.test], SLP.PC.scaled.test, SLP.grad.PC.scaled.test, SWH.bc.scaled[j, area.latSWH[k], idx.ar1], SWH.bc.scaled[j, area.latSWH[k], idx.ar2]))

                # Start with LASSO selection
                cv = cv.glmnet(as.matrix(predictors.training), SWH.bc.scaled[j, area.latSWH[k], idx.training], family = "gaussian", alpha = 1, nfold = 10)
                lasso = glmnet(as.matrix(predictors.training), SWH.bc.scaled[j, area.latSWH[k], idx.training], alpha = 1)
                minindex = which.min(abs(lasso$lambda - cv$lambda.min))    
                beta = lasso$beta[,minindex]
                predictors.training = predictors.training[, which( abs(beta) > 1e-6 )]
                predictors.test = predictors.test[, which( abs(beta) > 1e-6 )]
                
                fit <- lm(SWH.bc.scaled[j, area.latSWH[k], idx.training] ~ ., data = predictors.training) # Use arima(...) later
                fits = summary(fit)
                #fit <- arima(SWH.bc.scaled[j, area.latSWH[k], idx.training], order = c(1,0,0), xreg = predictors.training, method = "ML") # AR(1) noise
                
                # predict
                SWH.bc.scaled.pred <- predict(object = fit, newdata = predictors.test) #lm
                SWH.bc.scaled.pred.se <- fits$sigma
                #SWH.bc.scaled.dummy <- predict(fit, newxreg = predictors.test, se.fit = TRUE) #arima                
                #SWH.bc.scaled.pred <- SWH.bc.scaled.dummy$pred
                #SWH.bc.scaled.pred.se <- SWH.bc.scaled.dummy$se
                
                # Descale and detransform before computing different raknkings
                SWH.bc.pred <- SWH.bc.scaled.pred*SWH.bc.sd[j, area.latSWH[k]] # descale
                SWH.bc.pred.se <- SWH.bc.scaled.pred.se*SWH.bc.sd[j, area.latSWH[k]] # descale
                SWH.pred = InvBoxCox(SWH.bc.pred, SWH.bc.lambdas[j, area.latSWH[k]]) #detransform
                    
                SWH.bc.test = SWH.bc.scaled[j, area.latSWH[k], idx.test]*SWH.bc.sd[j, area.latSWH[k]] #descale
                SWH.test = InvBoxCox(SWH.bc.test, SWH.bc.lambdas[j, area.latSWH[k]]) #detransfrom

                # Handling NA's. Occurs rarely.
                not.NA <- !is.na(SWH.test) & !is.na(SWH.pred)
                if(sum(!not.NA ) > 0) {
                    cat("j =", j, "k =", k, "sum(is.na(SWH.test)) =", sum(is.na(SWH.test)), "sum(is.na(SWH.pred)) =", sum(is.na(SWH.pred)), "\n")
                }
                
                # To compute the Reliability index 
                #bc.sample = rnorm(n.bc.samp, mean = SWH.bc.pred[not.NA], sd = as.numeric(SWH.bc.pred.se))
                #sample = InvBoxCox(inv.bc.sample, SWH.bc.lambdas[j, area.latSWH[k]])
                #f.SWH = approxfun(density(sample)) # The distribution of the inverse of the normally distributed data
                #F.SWH = ecdf(sample) # The cum. distribution of the inverse of the normally distributed data
                
                rankj[k,1,p,1] = mean(rankModel(response = SWH.test[not.NA], fit, SWH.bc.pred[not.NA], model, "MAE", list(distr = "boxcox", lambda = SWH.bc.lambdas[j, area.latSWH[k]], se = SWH.bc.pred.se)), na.rm = TRUE)
                rankj[k,1,p,2] = mean(rankModel(response = SWH.test[not.NA], fit, SWH.bc.pred[not.NA], model, "LOGS", list(distr = "boxcox", lambda = SWH.bc.lambdas[j, area.latSWH[k]], se = SWH.bc.pred.se)), na.rm = TRUE)
                rankj[k,1,p,3] = reliabilityIndex(ppredbc(SWH.bc.test[not.NA], mean = SWH.bc.pred[not.NA], sd = as.numeric(SWH.bc.pred.se), lambda = SWH.bc.lambdas[j, area.latSWH[k]]), n.bins)
                rankj[k,1,p,4] = reliabilityIndexSquare(ppredbc(SWH.bc.test[not.NA], mean = SWH.bc.pred[not.NA], sd = as.numeric(SWH.bc.pred.se), lambda = SWH.bc.lambdas[j, area.latSWH[k]]), n.bins)
                rankj[k,1,p,5] = RMSE(SWH.test[not.NA], SWH.bc.pred[not.NA], as.numeric(SWH.bc.pred.se), SWH.bc.lambdas[j, area.latSWH[k]], n.bc.samp)
                rankj[k,1,p,6] = mean(CRPS(SWH.test[not.NA], SWH.bc.pred[not.NA], as.numeric(SWH.bc.pred.se), SWH.bc.lambdas[j, area.latSWH[k]], n.bc.samp), na.rm = TRUE)
                rankj[k,1,p,7] = mean(pSupport(mean = SWH.bc.pred[not.NA], sd = as.numeric(SWH.bc.pred.se), lambda = SWH.bc.lambdas[j, area.latSWH[k]]))

                #rankj.all[k,1,p,1,1:length(SWH.test[not.NA])] = rankModel(response = SWH.test[not.NA], fit, SWH.bc.pred[not.NA], model, "MAE", list(distr = "boxcox", lambda = SWH.bc.lambdas[j, area.latSWH[k]], se = SWH.bc.pred.se))
                #rankj.all[k,1,p,2,1:length(SWH.test[not.NA])] = rankModel(response = SWH.test[not.NA], fit, SWH.bc.pred[not.NA], model, "LOGS", list(distr = "boxcox", lambda = SWH.bc.lambdas[j, area.latSWH[k]], se = SWH.bc.pred.se))
                rankj.all[k,1,p,1,1:length(SWH.test[not.NA])] = ppredbc(SWH.bc.test[not.NA], mean = SWH.bc.pred[not.NA], sd = as.numeric(SWH.bc.pred.se), lambda = SWH.bc.lambdas[j, area.latSWH[k]])
                rankj.all[k,1,p,2,1:length(SWH.test[not.NA])] = NA
                #}                

                #U = ppredbc(SWH.bc.test[not.NA], mean = SWH.bc.pred[not.NA], sd = as.numeric(SWH.bc.pred.se), lambda = SWH.bc.lambdas[j, area.latSWH[k]])
                #save(U, file = paste("OutputComputations/Wang_season_", s, "j_", j, "k_", k, ".Rdata", sep = ""))
                #png(filename = "Figures/PIT_s1.png", height = 17, width = 17, units = "cm", res = 300)
                #hist(U, col = "light blue", breaks = 10, probability = TRUE, xlab = "PIT", main = paste("Longitude =", longitudeSWH[j], "Latitude =", latitudeSWH[k]))
                #dev.off()
                
                
                ##### Wang et al. 2012 No PC/EOF#####
                #cat("Wang et al. 2012 No PC/EOF, j =", j, "k =", area.latSWH[k], "\n")                
                idx.ar1 = c(rep(idx.training[1],1), idx.training[1:(length(idx.training)-1)])
                idx.ar2 = c(rep(idx.training[1],2), idx.training[1:(length(idx.training)-2)])
                predictors.training = as.data.frame(cbind(SLP.anom.scaled[idx.longSLP, idx.latSLP, idx.training], SLP.grad.bc.anom.scaled[idx.longSLP, idx.latSLP, idx.training], SWH.bc.scaled[j, area.latSWH[k], idx.ar1], SWH.bc.scaled[j, area.latSWH[k], idx.ar2]))
                
                idx.ar1 = c(rep(idx.test[1],1), idx.test[1:(length(idx.test)-1)])
                idx.ar2 = c(rep(idx.test[1],2), idx.test[1:(length(idx.test)-2)])                
                predictors.test = as.data.frame(cbind( SLP.anom.scaled[idx.longSLP, idx.latSLP, idx.test], SLP.grad.bc.anom.scaled[idx.longSLP, idx.latSLP, idx.test], SWH.bc.scaled[j, area.latSWH[k], idx.ar1], SWH.bc.scaled[j, area.latSWH[k], idx.ar2]))
                
                #if( sum(is.na(SWH.bc.scaled[j, area.latSWH[k], ])) + sum(is.na(predictors.test)) <= 10 & sum(is.na(predictors.training)) <= 5) {
                # Start with LASSO selection
                cv = cv.glmnet(as.matrix(predictors.training), SWH.bc.scaled[j, area.latSWH[k], idx.training], family = "gaussian", alpha = 1, nfold = 10)
                lasso = glmnet(as.matrix(predictors.training), SWH.bc.scaled[j, area.latSWH[k], idx.training], alpha = 1)
                minindex = which.min(abs(lasso$lambda - cv$lambda.min))    
                beta = lasso$beta[,minindex]
                predictors.training = predictors.training[, which( abs(beta) > 1e-6 )]
                predictors.test = predictors.test[, which( abs(beta) > 1e-6 )]
                
                fit <- lm(SWH.bc.scaled[j, area.latSWH[k], idx.training] ~ ., data = predictors.training) # Use arima(...) later
                fits = summary(fit)
                #fit <- arima(SWH.bc.scaled[j, area.latSWH[k], idx.training], order = c(1,0,0), xreg = predictors.training, method = "ML") # AR(1) noise
                
                # predict
                SWH.bc.scaled.pred <- predict(object = fit, newdata = predictors.test) #lm
                SWH.bc.scaled.pred.se <- fits$sigma
                #SWH.bc.scaled.dummy <- predict(fit, newxreg = predictors.test, se.fit = TRUE) #arima                
                #SWH.bc.scaled.pred <- SWH.bc.scaled.dummy$pred
                #SWH.bc.scaled.pred.se <- SWH.bc.scaled.dummy$se
                
                # Descale and detransform before computing different raknkings
                SWH.bc.pred <- SWH.bc.scaled.pred*SWH.bc.sd[j, area.latSWH[k]] # descale
                SWH.bc.pred.se <- SWH.bc.scaled.pred.se*SWH.bc.sd[j, area.latSWH[k]] # descale
                SWH.pred = InvBoxCox(SWH.bc.pred, SWH.bc.lambdas[j, area.latSWH[k]]) #detransform
                    
                SWH.bc.test = SWH.bc.scaled[j, area.latSWH[k], idx.test]*SWH.bc.sd[j, area.latSWH[k]] #descale
                SWH.test = InvBoxCox(SWH.bc.test, SWH.bc.lambdas[j, area.latSWH[k]]) #detransfrom

                # Handling NA's. Occurs rarely.
                not.NA <- !is.na(SWH.test) & !is.na(SWH.pred)
                if(sum(!not.NA ) > 0) {
                    cat("j =", j, "k =", k, "sum(is.na(SWH.test)) =", sum(is.na(SWH.test)), "sum(is.na(SWH.pred)) =", sum(is.na(SWH.pred)), "\n")
                }
                
                # To compute the Reliability index 
                #bc.sample = rnorm(n.bc.samp, mean = SWH.bc.pred[not.NA], sd = as.numeric(SWH.bc.pred.se))
                #sample = InvBoxCox(inv.bc.sample, SWH.bc.lambdas[j, area.latSWH[k]])
                #f.SWH = approxfun(density(sample)) # The distribution of the inverse of the normally distributed data
                #F.SWH = ecdf(sample) # The cum. distribution of the inverse of the normally distributed data

                rankj[k,2,p,1] = mean(rankModel(response = SWH.test[not.NA], fit, SWH.bc.pred[not.NA], model, "MAE", list(distr = "boxcox", lambda = SWH.bc.lambdas[j, area.latSWH[k]], se = SWH.bc.pred.se)), na.rm = TRUE)
                rankj[k,2,p,2] = mean(rankModel(response = SWH.test[not.NA], fit, SWH.bc.pred[not.NA], model, "LOGS", list(distr = "boxcox", lambda = SWH.bc.lambdas[j, area.latSWH[k]], se = SWH.bc.pred.se)), na.rm = TRUE)
                rankj[k,2,p,3] = reliabilityIndex(ppredbc(SWH.bc.test[not.NA], mean = SWH.bc.pred[not.NA], sd = as.numeric(SWH.bc.pred.se), lambda = SWH.bc.lambdas[j, area.latSWH[k]]), n.bins)
                rankj[k,2,p,4] = reliabilityIndexSquare(ppredbc(SWH.bc.test[not.NA], mean = SWH.bc.pred[not.NA], sd = as.numeric(SWH.bc.pred.se), lambda = SWH.bc.lambdas[j, area.latSWH[k]]), n.bins)
                rankj[k,2,p,5] = RMSE(SWH.test[not.NA], SWH.bc.pred[not.NA], as.numeric(SWH.bc.pred.se), SWH.bc.lambdas[j, area.latSWH[k]], n.bc.samp)
                rankj[k,2,p,6] = mean(CRPS(SWH.test[not.NA], SWH.bc.pred[not.NA], as.numeric(SWH.bc.pred.se), SWH.bc.lambdas[j, area.latSWH[k]], n.bc.samp), na.rm = TRUE)
                rankj[k,2,p,7] = mean(pSupport(mean = SWH.bc.pred[not.NA], sd = as.numeric(SWH.bc.pred.se), lambda = SWH.bc.lambdas[j, area.latSWH[k]]))
                
                #rankj.all[k,2,p,1,1:length(SWH.test[not.NA])] = rankModel(response = SWH.test[not.NA], fit, SWH.bc.pred[not.NA], model, "MAE", list(distr = "boxcox", lambda = SWH.bc.lambdas[j, area.latSWH[k]], se = SWH.bc.pred.se))
                #rankj.all[k,2,p,2,1:length(SWH.test[not.NA])] = rankModel(response = SWH.test[not.NA], fit, SWH.bc.pred[not.NA], model, "LOGS", list(distr = "boxcox", lambda = SWH.bc.lambdas[j, area.latSWH[k]], se = SWH.bc.pred.se))
                #rankj.all[k,2,p,3,1:length(SWH.test[not.NA])] = pnorm(SWH.bc.test[not.NA], mean = SWH.bc.pred[not.NA], sd = as.numeric(SWH.bc.pred.se))
                #rankj.all[k,2,p,4,1:length(SWH.test[not.NA])] = CRPS(SWH.test[not.NA], SWH.bc.pred[not.NA], as.numeric(SWH.bc.pred.se), SWH.bc.lambdas[j, area.latSWH[k]], n.bc.samp)                
                rankj.all[k,2,p,1,1:length(SWH.test[not.NA])] = ppredbc(SWH.bc.test[not.NA], mean = SWH.bc.pred[not.NA], sd = as.numeric(SWH.bc.pred.se), lambda = SWH.bc.lambdas[j, area.latSWH[k]])
                rankj.all[k,2,p,2,1:length(SWH.test[not.NA])] = NA
                #} 
            }
        }
    }
    return(list(rankj, rankj.all))
}


runParallelWangVanemWalker.fourier <- function(j) {
    rankj = array(dim = c(length(area.latSWH), n.models.fourier, n.training.test, n.rank.methods) )
    rankj.all =array(dim = c(length(area.latSWH), n.models.fourier, n.training.test, n.rank.methods.all, length(training.test[[1]][[2]])) )
    #rankj.all =array(dim = c(length(area.latSWH), n.models.fourier, n.training.test, n.rank.methods.all, length(training.test[[1]][[1]])) ) # When computing residuals
     
    for(p in 1:length(training.test)) {
        
        for(k in 1:length(area.latSWH)) {            
            ## Divide in training and test set
            idx.training = training.test[[p]][[1]]
            idx.test = training.test[[p]][[2]]
            
            idx.rank = 0
            idx.longSLP = which(longitudeSLP == longitudeSWH[j])
            idx.latSLP = which(latitudeSLP == latitudeSWH[area.latSWH[k]])

            if(sum(is.na(SWH[j, area.latSWH[k],])) < na.thresh & sum(is.na( SLP[idx.longSLP, idx.latSLP, ])) < na.thresh & sum(is.na( SLP.grad[idx.longSLP, idx.latSLP, ])) < na.thresh) {
                cat("Fourier, j =", j, "k =", area.latSWH[k], "\n")
                
                not.NA = which(!is.na(SWH[j, area.latSWH[k],]) & !is.na( SLP[idx.longSLP, idx.latSLP, ]) & !is.na( SLP.grad[idx.longSLP, idx.latSLP, ]))
                idx.training = idx.training[idx.training %in% not.NA]
                n.training = length(idx.training)
                idx.test = idx.test[idx.test %in% not.NA]
                n.test = length(idx.test)
            
                SWH.bc.dummy.training = BoxCoxLambda(SWH[j, area.latSWH[k],idx.training])
                SWH.bc.lambda.training = SWH.bc.dummy.training$lambda
                SWH.bc.training = SWH.bc.dummy.training$data
                SWH.bc.test = BoxCoxLambdaKnown2(SWH[j, area.latSWH[k],idx.test], SWH.bc.lambda.training)
                
                SLP.grad.bc.dummy.training = BoxCoxLambda(SLP.grad[idx.longSLP, idx.latSLP,idx.training])
                SLP.grad.bc.lambda.training = SLP.grad.bc.dummy.training$lambda
                SLP.grad.bc.training = SLP.grad.bc.dummy.training$data
                SLP.grad.bc.test = BoxCoxLambdaKnown2(SLP.grad[idx.longSLP, idx.latSLP,idx.test], SLP.grad.bc.lambda.training)
                
                # Standarizing the variables
                SWH.bc.mean.training = mean(SWH.bc.training)
                SWH.bc.sd.training = sd(SWH.bc.training)
                SWH.bc.standard.training = (SWH.bc.training - SWH.bc.mean.training)/SWH.bc.sd.training
                SWH.bc.standard.test = (SWH.bc.test - SWH.bc.mean.training)/SWH.bc.sd.training

                SLP.mean.training = mean( SLP[idx.longSLP, idx.latSLP, idx.training] )
                SLP.sd.training = sd( SLP[idx.longSLP, idx.latSLP, idx.training] )
                SLP.standard.training = (SLP[idx.longSLP, idx.latSLP, idx.training] - SLP.mean.training)/SLP.sd.training
                SLP.standard.test = (SLP[idx.longSLP, idx.latSLP, idx.test] - SLP.mean.training)/SLP.sd.training
                
                SLP.grad.bc.mean.training = mean(SLP.grad.bc.training)
                SLP.grad.bc.sd.training = sd(SLP.grad.bc.training)
                SLP.grad.bc.standard.training = (SLP.grad.bc.training - SLP.grad.bc.mean.training)/SLP.grad.bc.sd.training
                SLP.grad.bc.standard.test = (SLP.grad.bc.test - SLP.grad.bc.mean.training)/SLP.grad.bc.sd.training
                            
                # Compute Fourier periodics of the explanatory variables
                fourier.training = intercept.fourier[idx.training,]
                fourier.test = intercept.fourier[idx.test,]
                
                SLP.fourier.training = fourier.training*SLP.standard.training
                SLP.fourier.test = fourier.test*SLP.standard.test

                SLP.grad.bc.fourier.training = fourier.training*SLP.grad.bc.standard.training
                SLP.grad.bc.fourier.test = fourier.test*SLP.grad.bc.standard.test

                SWH.bc.fourier.training = fourier.training*SWH.bc.standard.training
                SWH.bc.fourier.test = fourier.test*SWH.bc.standard.test

                source("make_ar.r", local = TRUE)
                
                ##### Vanem&Walker spatial model LASSO #####
                for(neig in spatial.neighborhoods) {
                    #cat("Vanem&Walker spatial model LASSO, j =", j, "k =", area.latSWH[k], "neig =", neig, "\n")                
                    # Compute mean, max og mean of the neighborhood of the current point
                    SLP.spatmax = apply(SLP[ max(idx.longSLP-neig,1):min(idx.longSLP+neig, dim(SLP)[1]),max(idx.latSLP-neig,1):min(idx.latSLP+neig, dim(SLP)[2]),], 3, max, na.rm = T)
                    SLP.spatmin = apply(SLP[ max(idx.longSLP-neig,1):min(idx.longSLP+neig, dim(SLP)[1]),max(idx.latSLP-neig,1):min(idx.latSLP+neig, dim(SLP)[2]),], 3, min, na.rm = T)
                    SLP.spatmean = apply(SLP[ max(idx.longSLP-neig,1):min(idx.longSLP+neig, dim(SLP)[1]),max(idx.latSLP-neig,1):min(idx.latSLP+neig, dim(SLP)[2]),], 3, mean, na.rm = T)

                    SLP.grad.spatmax = apply(SLP.grad[ max(idx.longSLP-neig,1):min(idx.longSLP+neig, dim(SLP)[1]),max(idx.latSLP-neig,1):min(idx.latSLP+neig, dim(SLP)[2]),], 3, max, na.rm = T)
                    SLP.grad.spatmin = apply(SLP.grad[ max(idx.longSLP-neig,1):min(idx.longSLP+neig, dim(SLP)[1]),max(idx.latSLP-neig,1):min(idx.latSLP+neig, dim(SLP)[2]),], 3, min, na.rm = T)
                    SLP.grad.spatmean = apply(SLP.grad[ max(idx.longSLP-neig,1):min(idx.longSLP+neig, dim(SLP)[1]),max(idx.latSLP-neig,1):min(idx.latSLP+neig, dim(SLP)[2]),], 3, mean, na.rm = T)

                    # Split in test and training
                    SLP.spatmax.training = SLP.spatmax[idx.training]
                    SLP.spatmin.training = SLP.spatmin[idx.training]
                    SLP.spatmean.training = SLP.spatmean[idx.training]

                    SLP.grad.spatmax.training = SLP.grad.spatmax[idx.training]
                    SLP.grad.spatmin.training = SLP.grad.spatmin[idx.training]
                    SLP.grad.spatmean.training = SLP.grad.spatmean[idx.training]

                    SLP.spatmax.test = SLP.spatmax[idx.test]
                    SLP.spatmin.test = SLP.spatmin[idx.test]
                    SLP.spatmean.test = SLP.spatmean[idx.test]

                    SLP.grad.spatmax.test = SLP.grad.spatmax[idx.test]
                    SLP.grad.spatmin.test = SLP.grad.spatmin[idx.test]
                    SLP.grad.spatmean.test = SLP.grad.spatmean[idx.test]

                    # Standardize SLP
                    SLP.spatmax.mean.training = mean(SLP.spatmax.training)
                    SLP.spatmax.sd.training = sd(SLP.spatmax.training)
                    SLP.spatmax.standard.training = (SLP.spatmax.training - SLP.spatmax.mean.training)/SLP.spatmax.sd.training
                    SLP.spatmax.standard.test = (SLP.spatmax.test - SLP.spatmax.mean.training)/SLP.spatmax.sd.training
                    
                    SLP.spatmin.mean.training = mean(SLP.spatmin.training)
                    SLP.spatmin.sd.training = sd(SLP.spatmin.training)
                    SLP.spatmin.standard.training = (SLP.spatmin.training - SLP.spatmin.mean.training)/SLP.spatmin.sd.training
                    SLP.spatmin.standard.test = (SLP.spatmin.test - SLP.spatmin.mean.training)/SLP.spatmin.sd.training

                    SLP.spatmean.mean.training = mean(SLP.spatmean.training)
                    SLP.spatmean.sd.training = sd(SLP.spatmean.training)
                    SLP.spatmean.standard.training = (SLP.spatmean.training - SLP.spatmean.mean.training)/SLP.spatmean.sd.training
                    SLP.spatmean.standard.test = (SLP.spatmean.test - SLP.spatmean.mean.training)/SLP.spatmean.sd.training
                    
                    # Box Cox and standardize SLP.grad
                    SLP.grad.spatmax.bc.dummy = BoxCoxLambda(SLP.grad.spatmax.training)
                    SLP.grad.spatmax.bc.lambda.training = SLP.grad.spatmax.bc.dummy$lambda
                    SLP.grad.spatmax.bc.training = SLP.grad.spatmax.bc.dummy$data
                    SLP.grad.spatmax.bc.test = BoxCoxLambdaKnown2(SLP.grad.spatmax.test, SLP.grad.spatmax.bc.lambda.training)
                    
                    SLP.grad.spatmax.bc.mean.training = mean(SLP.grad.spatmax.bc.training)
                    SLP.grad.spatmax.bc.sd.training = sd(SLP.grad.spatmax.bc.training)
                    SLP.grad.spatmax.bc.standard.training = (SLP.grad.spatmax.bc.training - SLP.grad.spatmax.bc.mean.training)/SLP.grad.spatmax.bc.sd.training
                    SLP.grad.spatmax.bc.standard.test = (SLP.grad.spatmax.bc.test - SLP.grad.spatmax.bc.mean.training)/SLP.grad.spatmax.bc.sd.training

                    SLP.grad.spatmin.bc.dummy = BoxCoxLambda(SLP.grad.spatmin.training)
                    SLP.grad.spatmin.bc.lambda.training = SLP.grad.spatmin.bc.dummy$lambda
                    SLP.grad.spatmin.bc.training = SLP.grad.spatmin.bc.dummy$data
                    SLP.grad.spatmin.bc.test = BoxCoxLambdaKnown2(SLP.grad.spatmin.test, SLP.grad.spatmin.bc.lambda.training)
                    
                    SLP.grad.spatmin.bc.mean.training = mean(SLP.grad.spatmin.bc.training)
                    SLP.grad.spatmin.bc.sd.training = sd(SLP.grad.spatmin.bc.training)
                    SLP.grad.spatmin.bc.standard.training = (SLP.grad.spatmin.bc.training - SLP.grad.spatmin.bc.mean.training)/SLP.grad.spatmin.bc.sd.training
                    SLP.grad.spatmin.bc.standard.test = (SLP.grad.spatmin.bc.test - SLP.grad.spatmin.bc.mean.training)/SLP.grad.spatmin.bc.sd.training

                    SLP.grad.spatmean.bc.dummy = BoxCoxLambda(SLP.grad.spatmean.training)
                    SLP.grad.spatmean.bc.lambda.training = SLP.grad.spatmean.bc.dummy$lambda
                    SLP.grad.spatmean.bc.training = SLP.grad.spatmean.bc.dummy$data
                    SLP.grad.spatmean.bc.test = BoxCoxLambdaKnown2(SLP.grad.spatmean.test, SLP.grad.spatmean.bc.lambda.training)
                    
                    SLP.grad.spatmean.bc.mean.training = mean(SLP.grad.spatmean.bc.training)
                    SLP.grad.spatmean.bc.sd.training = sd(SLP.grad.spatmean.bc.training)
                    SLP.grad.spatmean.bc.standard.training = (SLP.grad.spatmean.bc.training - SLP.grad.spatmean.bc.mean.training)/SLP.grad.spatmean.bc.sd.training
                    SLP.grad.spatmean.bc.standard.test = (SLP.grad.spatmean.bc.test - SLP.grad.spatmean.bc.mean.training)/SLP.grad.spatmean.bc.sd.training

                    # Compute Fourier covariates
                    SLP.spatmax.fourier.training = fourier.training*SLP.spatmax.standard.training
                    SLP.spatmax.fourier.test = fourier.test*SLP.spatmax.standard.test

                    SLP.spatmin.fourier.training = fourier.training*SLP.spatmin.standard.training
                    SLP.spatmin.fourier.test = fourier.test*SLP.spatmin.standard.test

                    SLP.spatmean.fourier.training = fourier.training*SLP.spatmean.standard.training
                    SLP.spatmean.fourier.test = fourier.test*SLP.spatmean.standard.test

                    SLP.grad.spatmax.bc.fourier.training = fourier.training*SLP.grad.spatmax.bc.standard.training
                    SLP.grad.spatmax.bc.fourier.test = fourier.test*SLP.grad.spatmax.bc.standard.test
                    
                    SLP.grad.spatmin.bc.fourier.training = fourier.training*SLP.grad.spatmin.bc.standard.training
                    SLP.grad.spatmin.bc.fourier.test = fourier.test*SLP.grad.spatmin.bc.standard.test

                    SLP.grad.spatmean.bc.fourier.training = fourier.training*SLP.grad.spatmean.bc.standard.training
                    SLP.grad.spatmean.bc.fourier.test = fourier.test*SLP.grad.spatmean.bc.standard.test

                    predictors.training.spatial = as.data.frame(cbind(SLP.spatmax.standard.training, SLP.spatmin.standard.training, SLP.spatmean.standard.training, SLP.grad.spatmax.bc.standard.training, SLP.grad.spatmin.bc.standard.training, SLP.grad.spatmean.bc.standard.training, SLP.spatmax.fourier.training, SLP.spatmin.fourier.training, SLP.spatmean.fourier.training, SLP.grad.spatmax.bc.fourier.training, SLP.grad.spatmin.bc.fourier.training, SLP.grad.spatmean.bc.fourier.training))
                    
                    predictors.test.spatial = as.data.frame(cbind(SLP.spatmax.standard.test, SLP.spatmin.standard.test, SLP.spatmean.standard.test, SLP.grad.spatmax.bc.standard.test, SLP.grad.spatmin.bc.standard.test, SLP.grad.spatmean.bc.standard.test, SLP.spatmax.fourier.test, SLP.spatmin.fourier.test, SLP.spatmean.fourier.test, SLP.grad.spatmax.bc.fourier.test, SLP.grad.spatmin.bc.fourier.test, SLP.grad.spatmean.bc.fourier.test))
                    
                    predictors.training = as.data.frame( cbind(predictors.training.spatial, intercept.fourier[idx.training,], SLP.standard.training, SLP.grad.bc.standard.training, SLP.fourier.training, SLP.grad.bc.fourier.training, SWH.bc.ar.training.m1, SWH.bc.fourier.ar.training.m1, SWH.bc.ar.training.m2, SWH.bc.fourier.ar.training.m2) )

                    predictors.test = as.data.frame( cbind(predictors.test.spatial, intercept.fourier[idx.test,], SLP.standard.test, SLP.grad.bc.standard.test, SLP.fourier.test, SLP.grad.bc.fourier.test, SWH.bc.ar.test.m1, SWH.bc.fourier.ar.test.m1, SWH.bc.ar.test.m2, SWH.bc.fourier.ar.test.m2) )

                    source("LASSO_fit_predict_rank.bcfirst.r", local = TRUE)

                    predictors.training = as.data.frame( cbind(predictors.training.spatial, intercept.fourier[idx.training,], SLP.standard.training, SLP.grad.bc.standard.training, SLP.fourier.training, SLP.grad.bc.fourier.training, SWH.bc.ar.training.m1, SWH.bc.fourier.ar.training.m1, SWH.bc.ar.training.m2, SWH.bc.fourier.ar.training.m2, SWH.bc.ar.training.m3, SWH.bc.fourier.ar.training.m3, SWH.bc.ar.training.m4, SWH.bc.fourier.ar.training.m4, SWH.bc.ar.training.m5, SWH.bc.fourier.ar.training.m5) )
                    
                    predictors.test = as.data.frame( cbind(predictors.test.spatial, intercept.fourier[idx.test,], SLP.standard.test, SLP.grad.bc.standard.test, SLP.fourier.test, SLP.grad.bc.fourier.test, SWH.bc.ar.test.m1, SWH.bc.fourier.ar.test.m1, SWH.bc.ar.test.m2, SWH.bc.fourier.ar.test.m2, SWH.bc.ar.test.m3, SWH.bc.fourier.ar.test.m3, SWH.bc.ar.test.m4, SWH.bc.fourier.ar.test.m4, SWH.bc.ar.test.m5, SWH.bc.fourier.ar.test.m5) )

                    source("LASSO_fit_predict_rank.bcfirst.r", local = TRUE)

                    predictors.training = as.data.frame( cbind(predictors.training.spatial, intercept.fourier[idx.training,], SLP.standard.training, SLP.grad.bc.standard.training, SLP.fourier.training, SLP.grad.bc.fourier.training, SWH.bc.ar.training.m1, SWH.bc.fourier.ar.training.m1, SWH.bc.ar.training.m2, SWH.bc.fourier.ar.training.m2, SWH.bc.ar.training.m3, SWH.bc.fourier.ar.training.m3, SWH.bc.ar.training.m4, SWH.bc.fourier.ar.training.m4, SWH.bc.ar.training.m5, SWH.bc.fourier.ar.training.m5, SWH.bc.ar.training.m6, SWH.bc.fourier.ar.training.m6, SWH.bc.ar.training.m7, SWH.bc.fourier.ar.training.m7, SWH.bc.ar.training.m8, SWH.bc.fourier.ar.training.m8, SWH.bc.ar.training.m9, SWH.bc.fourier.ar.training.m9, SWH.bc.ar.training.m10, SWH.bc.fourier.ar.training.m10) )

                    predictors.test = as.data.frame( cbind(predictors.test.spatial, intercept.fourier[idx.test,], SLP.standard.test, SLP.grad.bc.standard.test, SLP.fourier.test, SLP.grad.bc.fourier.test, SWH.bc.ar.test.m1, SWH.bc.fourier.ar.test.m1, SWH.bc.ar.test.m2, SWH.bc.fourier.ar.test.m2, SWH.bc.ar.test.m3, SWH.bc.fourier.ar.test.m3, SWH.bc.ar.test.m4, SWH.bc.fourier.ar.test.m4, SWH.bc.ar.test.m5, SWH.bc.fourier.ar.test.m5, SWH.bc.ar.test.m6, SWH.bc.fourier.ar.test.m6, SWH.bc.ar.test.m7, SWH.bc.fourier.ar.test.m7, SWH.bc.ar.test.m8, SWH.bc.fourier.ar.test.m8, SWH.bc.ar.test.m9, SWH.bc.fourier.ar.test.m9, SWH.bc.ar.test.m10, SWH.bc.fourier.ar.test.m10) )

                    source("LASSO_fit_predict_rank.bcfirst.r", local = TRUE)

                }
            }
        }
    }
    return(list(rankj, rankj.all))
}


get.preddistr <- function(j) {
    preddistr.meanj = array(dim = c(length(area.latSWH), length(training.test[[1]][[2]])))
    preddistr.sdj = rep(NA, length(area.latSWH))
    preddistr.lambdaj = rep(NA, length(area.latSWH))
    
    for(p in 1:length(training.test)) {
        
        for(k in 1:length(area.latSWH)) {            
            ## Divide in training and test set
            idx.training = training.test[[p]][[1]]
            idx.test = training.test[[p]][[2]]
            
            idx.rank = 0
            idx.longSLP = which(longitudeSLP == longitudeSWH[j])
            idx.latSLP = which(latitudeSLP == latitudeSWH[area.latSWH[k]])

            if(sum(is.na(SWH[j, area.latSWH[k],])) < na.thresh & sum(is.na( SLP[idx.longSLP, idx.latSLP, ])) < na.thresh & sum(is.na( SLP.grad[idx.longSLP, idx.latSLP, ])) < na.thresh) {
                cat("Fourier, j =", j, "k =", area.latSWH[k], "\n")
                
                not.NA = which(!is.na(SWH[j, area.latSWH[k],]) & !is.na( SLP[idx.longSLP, idx.latSLP, ]) & !is.na( SLP.grad[idx.longSLP, idx.latSLP, ]))
                idx.training = idx.training[idx.training %in% not.NA]
                n.training = length(idx.training)
                idx.test = idx.test[idx.test %in% not.NA]
                n.test = length(idx.test)
            
                SWH.bc.dummy.training = BoxCoxLambda(SWH[j, area.latSWH[k],idx.training])
                SWH.bc.lambda.training = SWH.bc.dummy.training$lambda
                SWH.bc.training = SWH.bc.dummy.training$data
                SWH.bc.test = BoxCoxLambdaKnown2(SWH[j, area.latSWH[k],idx.test], SWH.bc.lambda.training)
                
                SLP.grad.bc.dummy.training = BoxCoxLambda(SLP.grad[idx.longSLP, idx.latSLP,idx.training])
                SLP.grad.bc.lambda.training = SLP.grad.bc.dummy.training$lambda
                SLP.grad.bc.training = SLP.grad.bc.dummy.training$data
                SLP.grad.bc.test = BoxCoxLambdaKnown2(SLP.grad[idx.longSLP, idx.latSLP,idx.test], SLP.grad.bc.lambda.training)
                
                # Standarizing the variables
                SWH.bc.mean.training = mean(SWH.bc.training)
                SWH.bc.sd.training = sd(SWH.bc.training)
                SWH.bc.standard.training = (SWH.bc.training - SWH.bc.mean.training)/SWH.bc.sd.training
                SWH.bc.standard.test = (SWH.bc.test - SWH.bc.mean.training)/SWH.bc.sd.training

                SLP.mean.training = mean( SLP[idx.longSLP, idx.latSLP, idx.training] )
                SLP.sd.training = sd( SLP[idx.longSLP, idx.latSLP, idx.training] )
                SLP.standard.training = (SLP[idx.longSLP, idx.latSLP, idx.training] - SLP.mean.training)/SLP.sd.training
                SLP.standard.test = (SLP[idx.longSLP, idx.latSLP, idx.test] - SLP.mean.training)/SLP.sd.training
                
                SLP.grad.bc.mean.training = mean(SLP.grad.bc.training)
                SLP.grad.bc.sd.training = sd(SLP.grad.bc.training)
                SLP.grad.bc.standard.training = (SLP.grad.bc.training - SLP.grad.bc.mean.training)/SLP.grad.bc.sd.training
                SLP.grad.bc.standard.test = (SLP.grad.bc.test - SLP.grad.bc.mean.training)/SLP.grad.bc.sd.training
                            
                # Compute Fourier periodics of the explanatory variables
                fourier.training = intercept.fourier[idx.training,]
                fourier.test = intercept.fourier[idx.test,]
                
                SLP.fourier.training = fourier.training*SLP.standard.training
                SLP.fourier.test = fourier.test*SLP.standard.test

                SLP.grad.bc.fourier.training = fourier.training*SLP.grad.bc.standard.training
                SLP.grad.bc.fourier.test = fourier.test*SLP.grad.bc.standard.test

                SWH.bc.fourier.training = fourier.training*SWH.bc.standard.training
                SWH.bc.fourier.test = fourier.test*SWH.bc.standard.test

                source("make_ar.r", local = TRUE)
                
                ##### Vanem&Walker spatial model LASSO #####
                for(neig in spatial.neighborhoods) {
                    #cat("Vanem&Walker spatial model LASSO, j =", j, "k =", area.latSWH[k], "neig =", neig, "\n")                
                    # Compute mean, max og mean of the neighborhood of the current point
                    SLP.spatmax = apply(SLP[ max(idx.longSLP-neig,1):min(idx.longSLP+neig, dim(SLP)[1]),max(idx.latSLP-neig,1):min(idx.latSLP+neig, dim(SLP)[2]),], 3, max, na.rm = T)
                    SLP.spatmin = apply(SLP[ max(idx.longSLP-neig,1):min(idx.longSLP+neig, dim(SLP)[1]),max(idx.latSLP-neig,1):min(idx.latSLP+neig, dim(SLP)[2]),], 3, min, na.rm = T)
                    SLP.spatmean = apply(SLP[ max(idx.longSLP-neig,1):min(idx.longSLP+neig, dim(SLP)[1]),max(idx.latSLP-neig,1):min(idx.latSLP+neig, dim(SLP)[2]),], 3, mean, na.rm = T)

                    SLP.grad.spatmax = apply(SLP.grad[ max(idx.longSLP-neig,1):min(idx.longSLP+neig, dim(SLP)[1]),max(idx.latSLP-neig,1):min(idx.latSLP+neig, dim(SLP)[2]),], 3, max, na.rm = T)
                    SLP.grad.spatmin = apply(SLP.grad[ max(idx.longSLP-neig,1):min(idx.longSLP+neig, dim(SLP)[1]),max(idx.latSLP-neig,1):min(idx.latSLP+neig, dim(SLP)[2]),], 3, min, na.rm = T)
                    SLP.grad.spatmean = apply(SLP.grad[ max(idx.longSLP-neig,1):min(idx.longSLP+neig, dim(SLP)[1]),max(idx.latSLP-neig,1):min(idx.latSLP+neig, dim(SLP)[2]),], 3, mean, na.rm = T)

                    # Split in test and training
                    SLP.spatmax.training = SLP.spatmax[idx.training]
                    SLP.spatmin.training = SLP.spatmin[idx.training]
                    SLP.spatmean.training = SLP.spatmean[idx.training]

                    SLP.grad.spatmax.training = SLP.grad.spatmax[idx.training]
                    SLP.grad.spatmin.training = SLP.grad.spatmin[idx.training]
                    SLP.grad.spatmean.training = SLP.grad.spatmean[idx.training]

                    SLP.spatmax.test = SLP.spatmax[idx.test]
                    SLP.spatmin.test = SLP.spatmin[idx.test]
                    SLP.spatmean.test = SLP.spatmean[idx.test]

                    SLP.grad.spatmax.test = SLP.grad.spatmax[idx.test]
                    SLP.grad.spatmin.test = SLP.grad.spatmin[idx.test]
                    SLP.grad.spatmean.test = SLP.grad.spatmean[idx.test]

                    # Standardize SLP
                    SLP.spatmax.mean.training = mean(SLP.spatmax.training)
                    SLP.spatmax.sd.training = sd(SLP.spatmax.training)
                    SLP.spatmax.standard.training = (SLP.spatmax.training - SLP.spatmax.mean.training)/SLP.spatmax.sd.training
                    SLP.spatmax.standard.test = (SLP.spatmax.test - SLP.spatmax.mean.training)/SLP.spatmax.sd.training
                    
                    SLP.spatmin.mean.training = mean(SLP.spatmin.training)
                    SLP.spatmin.sd.training = sd(SLP.spatmin.training)
                    SLP.spatmin.standard.training = (SLP.spatmin.training - SLP.spatmin.mean.training)/SLP.spatmin.sd.training
                    SLP.spatmin.standard.test = (SLP.spatmin.test - SLP.spatmin.mean.training)/SLP.spatmin.sd.training

                    SLP.spatmean.mean.training = mean(SLP.spatmean.training)
                    SLP.spatmean.sd.training = sd(SLP.spatmean.training)
                    SLP.spatmean.standard.training = (SLP.spatmean.training - SLP.spatmean.mean.training)/SLP.spatmean.sd.training
                    SLP.spatmean.standard.test = (SLP.spatmean.test - SLP.spatmean.mean.training)/SLP.spatmean.sd.training
                    
                    # Box Cox and standardize SLP.grad
                    SLP.grad.spatmax.bc.dummy = BoxCoxLambda(SLP.grad.spatmax.training)
                    SLP.grad.spatmax.bc.lambda.training = SLP.grad.spatmax.bc.dummy$lambda
                    SLP.grad.spatmax.bc.training = SLP.grad.spatmax.bc.dummy$data
                    SLP.grad.spatmax.bc.test = BoxCoxLambdaKnown2(SLP.grad.spatmax.test, SLP.grad.spatmax.bc.lambda.training)
                    
                    SLP.grad.spatmax.bc.mean.training = mean(SLP.grad.spatmax.bc.training)
                    SLP.grad.spatmax.bc.sd.training = sd(SLP.grad.spatmax.bc.training)
                    SLP.grad.spatmax.bc.standard.training = (SLP.grad.spatmax.bc.training - SLP.grad.spatmax.bc.mean.training)/SLP.grad.spatmax.bc.sd.training
                    SLP.grad.spatmax.bc.standard.test = (SLP.grad.spatmax.bc.test - SLP.grad.spatmax.bc.mean.training)/SLP.grad.spatmax.bc.sd.training

                    SLP.grad.spatmin.bc.dummy = BoxCoxLambda(SLP.grad.spatmin.training)
                    SLP.grad.spatmin.bc.lambda.training = SLP.grad.spatmin.bc.dummy$lambda
                    SLP.grad.spatmin.bc.training = SLP.grad.spatmin.bc.dummy$data
                    SLP.grad.spatmin.bc.test = BoxCoxLambdaKnown2(SLP.grad.spatmin.test, SLP.grad.spatmin.bc.lambda.training)
                    
                    SLP.grad.spatmin.bc.mean.training = mean(SLP.grad.spatmin.bc.training)
                    SLP.grad.spatmin.bc.sd.training = sd(SLP.grad.spatmin.bc.training)
                    SLP.grad.spatmin.bc.standard.training = (SLP.grad.spatmin.bc.training - SLP.grad.spatmin.bc.mean.training)/SLP.grad.spatmin.bc.sd.training
                    SLP.grad.spatmin.bc.standard.test = (SLP.grad.spatmin.bc.test - SLP.grad.spatmin.bc.mean.training)/SLP.grad.spatmin.bc.sd.training

                    SLP.grad.spatmean.bc.dummy = BoxCoxLambda(SLP.grad.spatmean.training)
                    SLP.grad.spatmean.bc.lambda.training = SLP.grad.spatmean.bc.dummy$lambda
                    SLP.grad.spatmean.bc.training = SLP.grad.spatmean.bc.dummy$data
                    SLP.grad.spatmean.bc.test = BoxCoxLambdaKnown2(SLP.grad.spatmean.test, SLP.grad.spatmean.bc.lambda.training)
                    
                    SLP.grad.spatmean.bc.mean.training = mean(SLP.grad.spatmean.bc.training)
                    SLP.grad.spatmean.bc.sd.training = sd(SLP.grad.spatmean.bc.training)
                    SLP.grad.spatmean.bc.standard.training = (SLP.grad.spatmean.bc.training - SLP.grad.spatmean.bc.mean.training)/SLP.grad.spatmean.bc.sd.training
                    SLP.grad.spatmean.bc.standard.test = (SLP.grad.spatmean.bc.test - SLP.grad.spatmean.bc.mean.training)/SLP.grad.spatmean.bc.sd.training

                    # Compute Fourier covariates
                    SLP.spatmax.fourier.training = fourier.training*SLP.spatmax.standard.training
                    SLP.spatmax.fourier.test = fourier.test*SLP.spatmax.standard.test

                    SLP.spatmin.fourier.training = fourier.training*SLP.spatmin.standard.training
                    SLP.spatmin.fourier.test = fourier.test*SLP.spatmin.standard.test

                    SLP.spatmean.fourier.training = fourier.training*SLP.spatmean.standard.training
                    SLP.spatmean.fourier.test = fourier.test*SLP.spatmean.standard.test

                    SLP.grad.spatmax.bc.fourier.training = fourier.training*SLP.grad.spatmax.bc.standard.training
                    SLP.grad.spatmax.bc.fourier.test = fourier.test*SLP.grad.spatmax.bc.standard.test
                    
                    SLP.grad.spatmin.bc.fourier.training = fourier.training*SLP.grad.spatmin.bc.standard.training
                    SLP.grad.spatmin.bc.fourier.test = fourier.test*SLP.grad.spatmin.bc.standard.test

                    SLP.grad.spatmean.bc.fourier.training = fourier.training*SLP.grad.spatmean.bc.standard.training
                    SLP.grad.spatmean.bc.fourier.test = fourier.test*SLP.grad.spatmean.bc.standard.test

                    predictors.training.spatial = as.data.frame(cbind(SLP.spatmax.standard.training, SLP.spatmin.standard.training, SLP.spatmean.standard.training, SLP.grad.spatmax.bc.standard.training, SLP.grad.spatmin.bc.standard.training, SLP.grad.spatmean.bc.standard.training, SLP.spatmax.fourier.training, SLP.spatmin.fourier.training, SLP.spatmean.fourier.training, SLP.grad.spatmax.bc.fourier.training, SLP.grad.spatmin.bc.fourier.training, SLP.grad.spatmean.bc.fourier.training))
                    
                    predictors.test.spatial = as.data.frame(cbind(SLP.spatmax.standard.test, SLP.spatmin.standard.test, SLP.spatmean.standard.test, SLP.grad.spatmax.bc.standard.test, SLP.grad.spatmin.bc.standard.test, SLP.grad.spatmean.bc.standard.test, SLP.spatmax.fourier.test, SLP.spatmin.fourier.test, SLP.spatmean.fourier.test, SLP.grad.spatmax.bc.fourier.test, SLP.grad.spatmin.bc.fourier.test, SLP.grad.spatmean.bc.fourier.test))
                    

                    predictors.training = as.data.frame( cbind(predictors.training.spatial, intercept.fourier[idx.training,], SLP.standard.training, SLP.grad.bc.standard.training, SLP.fourier.training, SLP.grad.bc.fourier.training, SWH.bc.ar.training.m1, SWH.bc.fourier.ar.training.m1, SWH.bc.ar.training.m2, SWH.bc.fourier.ar.training.m2, SWH.bc.ar.training.m3, SWH.bc.fourier.ar.training.m3, SWH.bc.ar.training.m4, SWH.bc.fourier.ar.training.m4, SWH.bc.ar.training.m5, SWH.bc.fourier.ar.training.m5) )
                    
                    predictors.test = as.data.frame( cbind(predictors.test.spatial, intercept.fourier[idx.test,], SLP.standard.test, SLP.grad.bc.standard.test, SLP.fourier.test, SLP.grad.bc.fourier.test, SWH.bc.ar.test.m1, SWH.bc.fourier.ar.test.m1, SWH.bc.ar.test.m2, SWH.bc.fourier.ar.test.m2, SWH.bc.ar.test.m3, SWH.bc.fourier.ar.test.m3, SWH.bc.ar.test.m4, SWH.bc.fourier.ar.test.m4, SWH.bc.ar.test.m5, SWH.bc.fourier.ar.test.m5) )

                    colnames(predictors.training) = paste("V", 1:dim(predictors.training)[2], sep = "")
                    colnames(predictors.test) = paste("V", 1:dim(predictors.test)[2], sep = "")

                    # Start with LASSO selection
                    cv = cv.glmnet(as.matrix(predictors.training), SWH.bc.standard.training, family = "gaussian", alpha = 1, nfold = 10)
                    lasso = glmnet(as.matrix(predictors.training), SWH.bc.standard.training, alpha = 1)
                    minindex = which.min(abs(lasso$lambda - cv$lambda.min))    
                    beta = lasso$beta[,minindex]
                    predictors.training = predictors.training[, which( abs(beta) > 1e-6 )]
                    predictors.test = predictors.test[, which( abs(beta) > 1e-6 )]

                    fit <- lm(SWH.bc.standard.training ~ ., data = predictors.training) # Use arima(...) later
                    fits = summary(fit)

                    # predict
                    SWH.bc.standard.pred <- predict(object = fit, newdata = predictors.test) #lm
                    SWH.bc.standard.pred.se <- fits$sigma

                    # Descale and detransform before computing different raknkings
                    SWH.bc.pred = SWH.bc.standard.pred*SWH.bc.sd.training + SWH.bc.mean.training
                    SWH.bc.pred.se = SWH.bc.standard.pred.se*SWH.bc.sd.training
                    SWH.pred = InvBoxCox(SWH.bc.pred, SWH.bc.lambda.training) #detransform

                    preddistr.meanj[k,idx.test - length(training.test[[1]][[1]])] = SWH.bc.pred
                    preddistr.sdj[k] = SWH.bc.pred.se
                    preddistr.lambdaj[k]= SWH.bc.lambda.training
                }
            }
        }
    }
    return(list(preddistr.meanj, preddistr.sdj, preddistr.lambdaj))
}
