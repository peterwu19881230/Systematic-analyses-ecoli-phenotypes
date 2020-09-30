#Create the obj: strain1 - strain 2 - (all annotations) - (all distances)

setwd('../..')
source("Nichols_preload.R")

load("Data/strain1_strain2.RData")
load("Data/strain1strain2_allAnnotations.RData")
load("Data/strain1strain2_allSimilarities.RData")

#Check if strain1 - strain2 are all the same for strain1strain2_allAnnotations, strain1strain2_allSimilarities 
##identical(strain1strain2_allAnnotations[,1:2],strain1strain2_allSimilarities[,1:2])
##I don't know why the above doesn't work (may be due to some hidden properties of the 2 objs)
##(In both obj: 1. strain1, strain2 are both integer 2. Rownames are the same)

identical(as.integer(strain1strain2_allSimilarities$strain1),strain1strain2_allAnnotations$strain1)
identical(as.integer(strain1strain2_allSimilarities$strain2),strain1strain2_allAnnotations$strain2)


strain1strain2_allAnnotations_allSimilarities=cbind(strain1strain2_allAnnotations,strain1strain2_allSimilarities[,-(1:2)],stringsAsFactors=F)

#Create additional T,F column where T stands for: 2 strains both have at least 1 significant phenotype
dim(sigQuantitative_Data_324cutoff_condCollapsed)

TF=apply(sigQuantitative_Data_324cutoff_condCollapsed,1,FUN=function(row){
  sum(row!=0)>0
})

index=rownames(sigQuantitative_Data_324cutoff_condCollapsed)[TF] %>% as.integer

sigOrNot=( (strain1_strain2$strain1 %in% index) & (strain1_strain2$strain2 %in% index) ) #2438736

strain1strain2_allAnnotations_allSimilarities$sigPhenoInBoth=sigOrNot




#delete so the next time these 2 lazy loading objs can be regenerated
file1="Data/sourced/strain1strain2_allAnnotations_allSimilarities.RData.rdb"
file2="Data/sourced/strain1strain2_allAnnotations_allSimilarities.RData.rdx"
if(file.exists(file1) & file.exists(file2)){
  file.remove(file1)
  file.remove(file2)
} 

#change some colnames
names(strain1strain2_allAnnotations_allSimilarities)[names(strain1strain2_allAnnotations_allSimilarities)=='regulator']='regulon'
names(strain1strain2_allAnnotations_allSimilarities)[names(strain1strain2_allAnnotations_allSimilarities)=='Pwy']='pwy'


save(strain1strain2_allAnnotations_allSimilarities,file="Data/sourced/strain1strain2_allAnnotations_allSimilarities.RData")









