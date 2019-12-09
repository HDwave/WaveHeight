# Calculates the Fourier terms for modeling seasonality
# @param x Coefficients
# @param K Number of Fourier terms   
# @param m Number of observations per period
# @return A matrix with 2xK Fourier terms 
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
