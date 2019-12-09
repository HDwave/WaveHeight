# Compute quantiles of the Box Cox distribution
qBoxCox = function(p, mean, sd, lambda) {
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
