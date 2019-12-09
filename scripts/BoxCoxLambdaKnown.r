# Performs Box-Cox transformation with a known lambda parameter 
# @param obs The observations to be transformed
# @param lambda The lambda parameter to be used in the transformation
# @return The observations transformed
BoxCoxLambdaKnown = function(obs, lambda) {
    if(round(lambda,3) == 0) {
        obs = log(obs)
    } else {
        obs = (obs^lambda - 1)/lambda
    }
    return(obs)
}
