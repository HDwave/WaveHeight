lag = 1
SWH.bc.ar.training.m1 = c(rep(SWH.bc.standard.training[1],lag), SWH.bc.standard.training[1:(n.training-lag)])
SWH.bc.fourier.ar.training.m1 = rbind(matrix(rep(SWH.bc.fourier.training[1,],lag), nrow = lag, byrow = TRUE), SWH.bc.fourier.training[1:(n.training-lag),])

lag = 2
SWH.bc.ar.training.m2 = c(rep(SWH.bc.standard.training[1],lag), SWH.bc.standard.training[1:(n.training-lag)])
SWH.bc.fourier.ar.training.m2 = rbind(matrix(rep(SWH.bc.fourier.training[1,],lag), nrow = lag, byrow = TRUE), SWH.bc.fourier.training[1:(n.training-lag),])

lag = 3
SWH.bc.ar.training.m3 = c(rep(SWH.bc.standard.training[1],lag), SWH.bc.standard.training[1:(n.training-lag)])
SWH.bc.fourier.ar.training.m3 = rbind(matrix(rep(SWH.bc.fourier.training[1,],lag), nrow = lag, byrow = TRUE), SWH.bc.fourier.training[1:(n.training-lag),])

lag = 4
SWH.bc.ar.training.m4 = c(rep(SWH.bc.standard.training[1],lag), SWH.bc.standard.training[1:(n.training-lag)])
SWH.bc.fourier.ar.training.m4 = rbind(matrix(rep(SWH.bc.fourier.training[1,],lag), nrow = lag, byrow = TRUE), SWH.bc.fourier.training[1:(n.training-lag),])

lag = 5
SWH.bc.ar.training.m5 = c(rep(SWH.bc.standard.training[1],lag), SWH.bc.standard.training[1:(n.training-lag)])
SWH.bc.fourier.ar.training.m5 = rbind(matrix(rep(SWH.bc.fourier.training[1,],lag), nrow = lag, byrow = TRUE), SWH.bc.fourier.training[1:(n.training-lag),])

lag = 6
SWH.bc.ar.training.m6 = c(rep(SWH.bc.standard.training[1],lag), SWH.bc.standard.training[1:(n.training-lag)])
SWH.bc.fourier.ar.training.m6 = rbind(matrix(rep(SWH.bc.fourier.training[1,],lag), nrow = lag, byrow = TRUE), SWH.bc.fourier.training[1:(n.training-lag),])

lag = 7
SWH.bc.ar.training.m7 = c(rep(SWH.bc.standard.training[1],lag), SWH.bc.standard.training[1:(n.training-lag)])
SWH.bc.fourier.ar.training.m7 = rbind(matrix(rep(SWH.bc.fourier.training[1,],lag), nrow = lag, byrow = TRUE), SWH.bc.fourier.training[1:(n.training-lag),])

lag = 8
SWH.bc.ar.training.m8 = c(rep(SWH.bc.standard.training[1],lag), SWH.bc.standard.training[1:(n.training-lag)])
SWH.bc.fourier.ar.training.m8 = rbind(matrix(rep(SWH.bc.fourier.training[1,],lag), nrow = lag, byrow = TRUE), SWH.bc.fourier.training[1:(n.training-lag),])

lag = 9
SWH.bc.ar.training.m9 = c(rep(SWH.bc.standard.training[1],lag), SWH.bc.standard.training[1:(n.training-lag)])
SWH.bc.fourier.ar.training.m9 = rbind(matrix(rep(SWH.bc.fourier.training[1,],lag), nrow = lag, byrow = TRUE), SWH.bc.fourier.training[1:(n.training-lag),])

lag = 10
SWH.bc.ar.training.m10 = c(rep(SWH.bc.standard.training[1],lag), SWH.bc.standard.training[1:(n.training-lag)])
SWH.bc.fourier.ar.training.m10 = rbind(matrix(rep(SWH.bc.fourier.training[1,],lag), nrow = lag, byrow = TRUE), SWH.bc.fourier.training[1:(n.training-lag),])


lag = 1
SWH.bc.ar.test.m1 = c(rep(SWH.bc.standard.test[1],lag), SWH.bc.standard.test[1:(n.test-lag)])
SWH.bc.fourier.ar.test.m1 = rbind(matrix(rep(SWH.bc.fourier.test[1,],lag), nrow = lag, byrow = TRUE), SWH.bc.fourier.test[1:(n.test-lag),])

lag = 2
SWH.bc.ar.test.m2 = c(rep(SWH.bc.standard.test[1],lag), SWH.bc.standard.test[1:(n.test-lag)])
SWH.bc.fourier.ar.test.m2 = rbind(matrix(rep(SWH.bc.fourier.test[1,],lag), nrow = lag, byrow = TRUE), SWH.bc.fourier.test[1:(n.test-lag),])

lag = 3
SWH.bc.ar.test.m3 = c(rep(SWH.bc.standard.test[1],lag), SWH.bc.standard.test[1:(n.test-lag)])
SWH.bc.fourier.ar.test.m3 = rbind(matrix(rep(SWH.bc.fourier.test[1,],lag), nrow = lag, byrow = TRUE), SWH.bc.fourier.test[1:(n.test-lag),])

lag = 4
SWH.bc.ar.test.m4 = c(rep(SWH.bc.standard.test[1],lag), SWH.bc.standard.test[1:(n.test-lag)])
SWH.bc.fourier.ar.test.m4 = rbind(matrix(rep(SWH.bc.fourier.test[1,],lag), nrow = lag, byrow = TRUE), SWH.bc.fourier.test[1:(n.test-lag),])

lag = 5
SWH.bc.ar.test.m5 = c(rep(SWH.bc.standard.test[1],lag), SWH.bc.standard.test[1:(n.test-lag)])
SWH.bc.fourier.ar.test.m5 = rbind(matrix(rep(SWH.bc.fourier.test[1,],lag), nrow = lag, byrow = TRUE), SWH.bc.fourier.test[1:(n.test-lag),])

lag = 6
SWH.bc.ar.test.m6 = c(rep(SWH.bc.standard.test[1],lag), SWH.bc.standard.test[1:(n.test-lag)])
SWH.bc.fourier.ar.test.m6 = rbind(matrix(rep(SWH.bc.fourier.test[1,],lag), nrow = lag, byrow = TRUE), SWH.bc.fourier.test[1:(n.test-lag),])

lag = 7
SWH.bc.ar.test.m7 = c(rep(SWH.bc.standard.test[1],lag), SWH.bc.standard.test[1:(n.test-lag)])
SWH.bc.fourier.ar.test.m7 = rbind(matrix(rep(SWH.bc.fourier.test[1,],lag), nrow = lag, byrow = TRUE), SWH.bc.fourier.test[1:(n.test-lag),])

lag = 8
SWH.bc.ar.test.m8 = c(rep(SWH.bc.standard.test[1],lag), SWH.bc.standard.test[1:(n.test-lag)])
SWH.bc.fourier.ar.test.m8 = rbind(matrix(rep(SWH.bc.fourier.test[1,],lag), nrow = lag, byrow = TRUE), SWH.bc.fourier.test[1:(n.test-lag),])

lag = 9
SWH.bc.ar.test.m9 = c(rep(SWH.bc.standard.test[1],lag), SWH.bc.standard.test[1:(n.test-lag)])
SWH.bc.fourier.ar.test.m9 = rbind(matrix(rep(SWH.bc.fourier.test[1,],lag), nrow = lag, byrow = TRUE), SWH.bc.fourier.test[1:(n.test-lag),])

lag = 10
SWH.bc.ar.test.m10 = c(rep(SWH.bc.standard.test[1],lag), SWH.bc.standard.test[1:(n.test-lag)])
SWH.bc.fourier.ar.test.m10 = rbind(matrix(rep(SWH.bc.fourier.test[1,],lag), nrow = lag, byrow = TRUE), SWH.bc.fourier.test[1:(n.test-lag),])
