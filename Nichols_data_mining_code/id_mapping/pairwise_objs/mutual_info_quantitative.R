#pairwise MI for quantitative data

setwd('../..')
source("Nichols_preload.R")

install.packages('mpmi')


start.time = Sys.time()
mi_all=mpmi::cminjk(t(All_Data_NAimputed))
end.time = Sys.time()
end.time - start.time  #Time difference of 27.39689 mins

mi_all=mi_all/log(2) #convert nats to bits
class(mi_all)
dim(mi_all)

save(mi_all,file="Data/mi_all.RData")