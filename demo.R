source('GeX_Meth_Discretiser_Library.R')
load('signal.rdata')

hist(signal,100)
hist(EL_statistic(signal))

print(detect_bimodal_signal((signal)))