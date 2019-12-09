# Compute pit values of the Box Cox distribution
pBoxCox = function(x, mean, sd, lambda) {
    g.x  <- BoxCoxLambdaKnown(x, lambda)
    if(round(lambda,2) < 0) {
        p = pnorm((g.x-mean)/sd) / pnorm((-1/lambda - mean)/sd)
    } else if(round(lambda,2) == 0) {
        p = pnorm((g.x-mean)/sd)
    } else { # lambda > 0
        p = (pnorm((g.x-mean)/sd) - pnorm((-1/lambda - mean)/sd)) / pnorm((1/lambda + mean)/sd)
    }
    return(p)
}
