#execute by Rscript .... in terminal
##This script takes around 1.5hr on my PC


setwd('../..')
source("Nichols_preload.R")

load("Data/uniqueChemIndex.RData")
TF=( names(uniqueChemIndex) %in% c("NH4Cl (MOPS)","Iron excess-FeSO4","Iron starvation-FeSO4","Acetate (M9)",
                                   "Glucosamine (M9)","Glucose (M9)","Glycerol (M9)","Maltose (M9)","N-acetyl Glucosamine","Succinate (M9)") )
used_cond=names(uniqueChemIndex)[!TF]
cond_indices=unlist(uniqueChemIndex[used_cond])

dat=All_Data_NAimputed[,cond_indices]


#Create a data frame that has strain1 -  strain2 - PCC -Spearmen...(All meaningful distances)


#Quantitative:
##PCC
pcc= (cor(t(dat))) %>% as.dist %>% melt_dist

##Spearman
spearman=cor(t(dat),method="spearman") %>% as.dist %>% melt_dist

#Note: Above pcc and spearman based distance are not real distance. It's just for my own convenience => real distance should be 1 - corr (no absolute value)


##Euclidean
#euclidean=get_dist(dat,method="euclidean") %>% melt_dist

##Mutual information
start.time = Sys.time()
mi_all=mpmi::cminjk(t(dat))
end.time = Sys.time()
end.time - start.time  #Time difference of 27.39689 mins
save(mi_all,file="Data/mi_all_noMinimalMedia.RData")
load("Data/mi_all_noMinimalMedia.RData")
mi_all=mi_all/log(2) #convert nats to bits
mi=mi_all %>% as.dist %>% melt_dist





#Qualitative:

dat=Ternary_Data_324cutff_NAremoved[,cond_indices]

##Euclidean
#euclidean_qualitative=get_dist(dat,method="euclidean") %>% melt_dist

#Not sure this is needed for the paper
#cond_indices=?
#euclidean_collapsedCond=get_dist(Ternary_Data_324cutff_condCollapsed[,cond_indices],method="euclidean") %>% melt_dist


##Mutual Information (Sho)
##Ref: https://cran.r-project.org/web/packages/infotheo/infotheo.pdf
start.time = Sys.time()
mi_ternary=(infotheo::mutinformation(t(dat) %>% as.data.frame) %>% infotheo::natstobits() ) %>% melt_dist 
end.time = Sys.time()
end.time-start.time #Time difference of 1.059961 hours
save(mi_ternary,file="Data/mi_ternary_noMinimalMedia.RData")
load("Data/mi_ternary_noMinimalMedia.RData")


##MI for collapsed conditions. Not sure this is needed for the paper

#dat=Ternary_Data_324cutff_condCollapsed[,]
#cond_indices=?
#start.time = Sys.time()
#mi_ternary_collapsedCond=( 1- mutinformation(t(dat) %>% as.data.frame) %>% natstobits ) %>% melt_dist 
#end.time = Sys.time()
#end.time-start.time #Time difference of 22.47311 mins
#save(mi_ternary_collapsedCond,file="Data/mi_ternary_collapsedCond_noMinimalMedia.RData")
#load("Data/mi_ternary_collapsedCond_noMinimalMedia.RData")





#If the first 2 columns are the same across all distance dataframes here, I can simply use cbind()
identical(pcc[,1:2],spearman[,1:2])
#identical(spearman[,1:2],euclidean[,1:2])
#identical(euclidean[,1:2],mi[,1:2])
identical(spearman[,1:2],mi[,1:2])
#identical(euclidean[,1:2],euclidean_qualitative[,1:2])
#identical(euclidean_qualitative[,1:2],mi_ternary[,1:2])
identical(mi[,1:2],mi_ternary[,1:2])


strain1strain2_allSimilarities=cbind(pcc,
                                     spearman=spearman$value,
                                     #euclidean=euclidean$value,
                                     mi=mi$value,
                                     #euclidean_qualitative=euclidean_qualitative$value,
                                     mi_ternary=mi_ternary$value,
                                     stringsAsFactors=F)

colnames(strain1strain2_allSimilarities)[1:3]=c("strain1","strain2","pcc")

strain1strain2_allSimilarities$strain1=as.integer(strain1strain2_allSimilarities$strain1)
strain1strain2_allSimilarities$strain2=as.integer(strain1strain2_allSimilarities$strain2)

for(i in 3:dim(strain1strain2_allSimilarities)[2]){
  strain1strain2_allSimilarities[,i]=as.numeric(strain1strain2_allSimilarities[,i])
}

rownames(strain1strain2_allSimilarities)=as.character(1:7914231) #correct the rownames (don't know why it started with 2)

strain1strain2_allSimilarities_noMinimalMedia=strain1strain2_allSimilarities
save(strain1strain2_allSimilarities_noMinimalMedia,file="Data/strain1strain2_allSimilarities_noMinimalMedia.RData")








