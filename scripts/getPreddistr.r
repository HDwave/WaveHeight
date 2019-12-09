
# Estimates model parameters and generates predictive distributions for all time points in the test set
# @param j Longitude index
# @param k Latitude index   
# @param neig Size of spatial neighorhood for spatial covariates
# @return A list of predictive means, predictive spread and lambda for Box-Cox transformation
getPreddistr <- function(j, k, neig) {
    pred.mean = rep(NA, length(training.test[[2]]))
    pred.sd = NA
    pred.lambda = NA
    
    ## Divide in training and test set
    idx.training = training.test[[1]]
    idx.test = training.test[[2]]
    
    idx.rank = 0

    ## Make sure covariates are from the correct location
    idx.longSLP = which(longitudeSLP == longitudeSWH[j])
    idx.latSLP = which(latitudeSLP == latitudeSWH[k])

    ## Only continue with analysis if sufficient available data
    if(sum(is.na(SWH[j, k,])) < na.thresh &
       sum(is.na( SLP[idx.longSLP, idx.latSLP, ])) < na.thresh &
       sum(is.na( SLP.grad[idx.longSLP, idx.latSLP, ])) < na.thresh) {
        
        cat("Fourier: j =", j, ", k =", k, "\n")
        
        ## Remove missing data 
        not.NA = which(!is.na(SWH[j, k,]) &
                       !is.na( SLP[idx.longSLP, idx.latSLP, ]) &
                       !is.na( SLP.grad[idx.longSLP, idx.latSLP, ]))
        idx.training = idx.training[idx.training %in% not.NA]
        n.training = length(idx.training)
        idx.test = idx.test[idx.test %in% not.NA]
        n.test = length(idx.test)
        
        ## Estimate lambda and Box-Cox transform SWH in training period
        SWH.bc.dummy.training = BoxCoxLambda(SWH[j, k,idx.training])
        SWH.bc.lambda.training = SWH.bc.dummy.training$lambda
        SWH.bc.training = SWH.bc.dummy.training$data
        ## Box-Cox transform SWH in test period with same lambda
        SWH.bc.test = BoxCoxLambdaKnown(SWH[j, k,idx.test], SWH.bc.lambda.training)
        
        ## Estimate lambda and Box-Cox transform SLP.grad in training period 
        SLP.grad.bc.dummy.training = BoxCoxLambda(SLP.grad[idx.longSLP, idx.latSLP,idx.training])
        SLP.grad.bc.lambda.training = SLP.grad.bc.dummy.training$lambda
        SLP.grad.bc.training = SLP.grad.bc.dummy.training$data
        ## Box-Cox transform SLP.grad in test period with same lambda
        SLP.grad.bc.test = BoxCoxLambdaKnown(SLP.grad[idx.longSLP, idx.latSLP,idx.test], SLP.grad.bc.lambda.training)
        
        ## Standardizing the variables
        SWH.bc.mean.training = mean(SWH.bc.training)
        SWH.bc.sd.training = sd(SWH.bc.training)
        SWH.bc.standard.training = (SWH.bc.training - SWH.bc.mean.training)/
            SWH.bc.sd.training
        SWH.bc.standard.test = (SWH.bc.test - SWH.bc.mean.training)/
            SWH.bc.sd.training

        SLP.mean.training = mean( SLP[idx.longSLP, idx.latSLP, idx.training] )
        SLP.sd.training = sd( SLP[idx.longSLP, idx.latSLP, idx.training] )
        SLP.standard.training = (SLP[idx.longSLP, idx.latSLP, idx.training] - SLP.mean.training)/
            SLP.sd.training
        SLP.standard.test = (SLP[idx.longSLP, idx.latSLP, idx.test] - SLP.mean.training)/
            SLP.sd.training
        
        SLP.grad.bc.mean.training = mean(SLP.grad.bc.training)
        SLP.grad.bc.sd.training = sd(SLP.grad.bc.training)
        SLP.grad.bc.standard.training = (SLP.grad.bc.training - SLP.grad.bc.mean.training)/
            SLP.grad.bc.sd.training
        SLP.grad.bc.standard.test = (SLP.grad.bc.test - SLP.grad.bc.mean.training)/
            SLP.grad.bc.sd.training
        
        ## Compute Fourier periodics of the explanatory variables
        fourier.training = intercept.fourier[idx.training,]
        fourier.test = intercept.fourier[idx.test,]
        
        SLP.fourier.training = fourier.training*SLP.standard.training
        SLP.fourier.test = fourier.test*SLP.standard.test
        
        SLP.grad.bc.fourier.training = fourier.training*SLP.grad.bc.standard.training
        SLP.grad.bc.fourier.test = fourier.test*SLP.grad.bc.standard.test
        
        SWH.bc.fourier.training = fourier.training*SWH.bc.standard.training
        SWH.bc.fourier.test = fourier.test*SWH.bc.standard.test
        
        source("make_ar.r", local = TRUE)
                
        ## Vanem&Walker spatial model LASSO #####

        ## Compute mean, max og mean of the neighborhood of the current point
        SLP.spatmax = apply(SLP[ max(idx.longSLP-neig,1):min(idx.longSLP+neig, dim(SLP)[1]),
                                max(idx.latSLP-neig,1):min(idx.latSLP+neig, dim(SLP)[2]),], 3, max, na.rm = T)
        SLP.spatmin = apply(SLP[ max(idx.longSLP-neig,1):min(idx.longSLP+neig, dim(SLP)[1]),
                                max(idx.latSLP-neig,1):min(idx.latSLP+neig, dim(SLP)[2]),], 3, min, na.rm = T)
        SLP.spatmean = apply(SLP[ max(idx.longSLP-neig,1):min(idx.longSLP+neig, dim(SLP)[1]),
                                 max(idx.latSLP-neig,1):min(idx.latSLP+neig, dim(SLP)[2]),], 3, mean, na.rm = T)
        
        SLP.grad.spatmax = apply(SLP.grad[ max(idx.longSLP-neig,1):min(idx.longSLP+neig, dim(SLP)[1]),
                                          max(idx.latSLP-neig,1):min(idx.latSLP+neig, dim(SLP)[2]),], 3, max, na.rm = T)
        SLP.grad.spatmin = apply(SLP.grad[ max(idx.longSLP-neig,1):min(idx.longSLP+neig, dim(SLP)[1]),
                                          max(idx.latSLP-neig,1):min(idx.latSLP+neig, dim(SLP)[2]),], 3, min, na.rm = T)
        SLP.grad.spatmean = apply(SLP.grad[ max(idx.longSLP-neig,1):min(idx.longSLP+neig, dim(SLP)[1]),
                                           max(idx.latSLP-neig,1):min(idx.latSLP+neig, dim(SLP)[2]),], 3, mean, na.rm = T)
        
        ## Split in test and training
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
        
        ## Standardize SLP spatial variables
        SLP.spatmax.mean.training = mean(SLP.spatmax.training)
        SLP.spatmax.sd.training = sd(SLP.spatmax.training)
        SLP.spatmax.standard.training = (SLP.spatmax.training - SLP.spatmax.mean.training)/
            SLP.spatmax.sd.training
        SLP.spatmax.standard.test = (SLP.spatmax.test - SLP.spatmax.mean.training)/
            SLP.spatmax.sd.training
        
        SLP.spatmin.mean.training = mean(SLP.spatmin.training)
        SLP.spatmin.sd.training = sd(SLP.spatmin.training)
        SLP.spatmin.standard.training = (SLP.spatmin.training - SLP.spatmin.mean.training)/
            SLP.spatmin.sd.training
        SLP.spatmin.standard.test = (SLP.spatmin.test - SLP.spatmin.mean.training)/
            SLP.spatmin.sd.training
        
        SLP.spatmean.mean.training = mean(SLP.spatmean.training)
        SLP.spatmean.sd.training = sd(SLP.spatmean.training)
        SLP.spatmean.standard.training = (SLP.spatmean.training - SLP.spatmean.mean.training)/
            SLP.spatmean.sd.training
        SLP.spatmean.standard.test = (SLP.spatmean.test - SLP.spatmean.mean.training)/
            SLP.spatmean.sd.training
        
        ## Box-Cox transform and standardize SLP.grad
        SLP.grad.spatmax.bc.dummy = BoxCoxLambda(SLP.grad.spatmax.training)
        SLP.grad.spatmax.bc.lambda.training = SLP.grad.spatmax.bc.dummy$lambda
        SLP.grad.spatmax.bc.training = SLP.grad.spatmax.bc.dummy$data
        SLP.grad.spatmax.bc.test = BoxCoxLambdaKnown(SLP.grad.spatmax.test, SLP.grad.spatmax.bc.lambda.training)
        
        SLP.grad.spatmax.bc.mean.training = mean(SLP.grad.spatmax.bc.training)
        SLP.grad.spatmax.bc.sd.training = sd(SLP.grad.spatmax.bc.training)
        SLP.grad.spatmax.bc.standard.training = (SLP.grad.spatmax.bc.training - SLP.grad.spatmax.bc.mean.training)/
            SLP.grad.spatmax.bc.sd.training
        SLP.grad.spatmax.bc.standard.test = (SLP.grad.spatmax.bc.test - SLP.grad.spatmax.bc.mean.training)/
            SLP.grad.spatmax.bc.sd.training
        
        SLP.grad.spatmin.bc.dummy = BoxCoxLambda(SLP.grad.spatmin.training)
        SLP.grad.spatmin.bc.lambda.training = SLP.grad.spatmin.bc.dummy$lambda
        SLP.grad.spatmin.bc.training = SLP.grad.spatmin.bc.dummy$data
        SLP.grad.spatmin.bc.test = BoxCoxLambdaKnown(SLP.grad.spatmin.test, SLP.grad.spatmin.bc.lambda.training)
        
        SLP.grad.spatmin.bc.mean.training = mean(SLP.grad.spatmin.bc.training)
        SLP.grad.spatmin.bc.sd.training = sd(SLP.grad.spatmin.bc.training)
        SLP.grad.spatmin.bc.standard.training = (SLP.grad.spatmin.bc.training - SLP.grad.spatmin.bc.mean.training)/
            SLP.grad.spatmin.bc.sd.training
        SLP.grad.spatmin.bc.standard.test = (SLP.grad.spatmin.bc.test - SLP.grad.spatmin.bc.mean.training)/
            SLP.grad.spatmin.bc.sd.training
        
        SLP.grad.spatmean.bc.dummy = BoxCoxLambda(SLP.grad.spatmean.training)
        SLP.grad.spatmean.bc.lambda.training = SLP.grad.spatmean.bc.dummy$lambda
        SLP.grad.spatmean.bc.training = SLP.grad.spatmean.bc.dummy$data
        SLP.grad.spatmean.bc.test = BoxCoxLambdaKnown(SLP.grad.spatmean.test, SLP.grad.spatmean.bc.lambda.training)
        
        SLP.grad.spatmean.bc.mean.training = mean(SLP.grad.spatmean.bc.training)
        SLP.grad.spatmean.bc.sd.training = sd(SLP.grad.spatmean.bc.training)
        SLP.grad.spatmean.bc.standard.training = (SLP.grad.spatmean.bc.training - SLP.grad.spatmean.bc.mean.training)/
            SLP.grad.spatmean.bc.sd.training
        SLP.grad.spatmean.bc.standard.test = (SLP.grad.spatmean.bc.test - SLP.grad.spatmean.bc.mean.training)/
            SLP.grad.spatmean.bc.sd.training
        
        ## Compute Fourier predictors for seasonally varying coefficients
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

        ## Create full set of spatial predictors for the training and test sets 
        predictors.training.spatial = as.data.frame(cbind(SLP.spatmax.standard.training,
                                                          SLP.spatmin.standard.training,
                                                          SLP.spatmean.standard.training,
                                                          SLP.grad.spatmax.bc.standard.training,
                                                          SLP.grad.spatmin.bc.standard.training,
                                                          SLP.grad.spatmean.bc.standard.training,
                                                          SLP.spatmax.fourier.training,
                                                          SLP.spatmin.fourier.training,
                                                          SLP.spatmean.fourier.training,
                                                          SLP.grad.spatmax.bc.fourier.training,
                                                          SLP.grad.spatmin.bc.fourier.training,
                                                          SLP.grad.spatmean.bc.fourier.training))
                    
        predictors.test.spatial = as.data.frame(cbind(SLP.spatmax.standard.test,
                                                      SLP.spatmin.standard.test,
                                                      SLP.spatmean.standard.test,
                                                      SLP.grad.spatmax.bc.standard.test,
                                                      SLP.grad.spatmin.bc.standard.test,
                                                      SLP.grad.spatmean.bc.standard.test,
                                                      SLP.spatmax.fourier.test,
                                                      SLP.spatmin.fourier.test,
                                                      SLP.spatmean.fourier.test,
                                                      SLP.grad.spatmax.bc.fourier.test,
                                                      SLP.grad.spatmin.bc.fourier.test,
                                                      SLP.grad.spatmean.bc.fourier.test))
        
        ## Create set of local predictors with AR order 5
        predictors.training = as.data.frame( cbind(predictors.training.spatial,
                                                   intercept.fourier[idx.training,],
                                                   SLP.standard.training,
                                                   SLP.grad.bc.standard.training,
                                                   SLP.fourier.training,
                                                   SLP.grad.bc.fourier.training,
                                                   SWH.bc.ar.training.m1,
                                                   SWH.bc.fourier.ar.training.m1,
                                                   SWH.bc.ar.training.m2,
                                                   SWH.bc.fourier.ar.training.m2,
                                                   SWH.bc.ar.training.m3,
                                                   SWH.bc.fourier.ar.training.m3,
                                                   SWH.bc.ar.training.m4,
                                                   SWH.bc.fourier.ar.training.m4,
                                                   SWH.bc.ar.training.m5,
                                                   SWH.bc.fourier.ar.training.m5) )
        
        predictors.test = as.data.frame( cbind(predictors.test.spatial,
                                               intercept.fourier[idx.test,],
                                               SLP.standard.test,
                                               SLP.grad.bc.standard.test,
                                               SLP.fourier.test,
                                               SLP.grad.bc.fourier.test,
                                               SWH.bc.ar.test.m1,
                                               SWH.bc.fourier.ar.test.m1,
                                               SWH.bc.ar.test.m2,
                                               SWH.bc.fourier.ar.test.m2,
                                               SWH.bc.ar.test.m3,
                                               SWH.bc.fourier.ar.test.m3,
                                               SWH.bc.ar.test.m4,
                                               SWH.bc.fourier.ar.test.m4,
                                               SWH.bc.ar.test.m5,
                                               SWH.bc.fourier.ar.test.m5) )
        
        colnames(predictors.training) = paste("V", 1:dim(predictors.training)[2], sep = "")
        colnames(predictors.test) = paste("V", 1:dim(predictors.test)[2], sep = "")
        cat("Total number of potential predictors:",dim(predictors.training)[2], "\n")
        
        ## Start with LASSO selection
        cv = cv.glmnet(as.matrix(predictors.training), SWH.bc.standard.training, family = "gaussian", alpha = 1, nfold = 10)
        lasso = glmnet(as.matrix(predictors.training), SWH.bc.standard.training, alpha = 1)
        minindex = which.min(abs(lasso$lambda - cv$lambda.min))    
        beta = lasso$beta[,minindex]
        predictors.training = predictors.training[, which( abs(beta) > 1e-6 )]
        predictors.test = predictors.test[, which( abs(beta) > 1e-6 )]
        cat("Number of predictors selected by LASSO:",dim(predictors.training)[2], "\n")
      
        fit <- lm(SWH.bc.standard.training ~ ., data = predictors.training) 
        fits = summary(fit)
        
        ## predict
        SWH.bc.standard.pred <- predict(object = fit, newdata = predictors.test) #lm
        SWH.bc.standard.pred.se <- fits$sigma
        
        ## Descale and detransform before computing different raknkings
        SWH.bc.pred = SWH.bc.standard.pred*SWH.bc.sd.training + SWH.bc.mean.training
        SWH.bc.pred.se = SWH.bc.standard.pred.se*SWH.bc.sd.training
        SWH.pred = InvBoxCox(SWH.bc.pred, SWH.bc.lambda.training) #detransform
        
        pred.mean[idx.test - length(training.test[[1]])] = SWH.bc.pred
        pred.sd = SWH.bc.pred.se
        pred.lambda = SWH.bc.lambda.training
    }
    return(list(pred.SWH=SWH.pred, pred.mean=pred.mean, pred.sd=pred.sd, pred.lambda=pred.lambda, fits=fits))
}
