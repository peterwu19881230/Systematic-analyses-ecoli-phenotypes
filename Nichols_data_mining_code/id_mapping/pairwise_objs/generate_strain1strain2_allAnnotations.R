#(run this using terminal: Rscript script.R)
#last time I ran on my 13' mac took 2.680446 mins
setwd("/Users/peterwu/Dropbox/Nichols_Data_mining/")
source("Nichols_preload.R")

#Construct strain1strain2_allAnnotations 

#dependency: id_allAttributes


#1st and 2nd column for strain pairs
strain1strain2_allAnnotations=as.data.frame(t(combn(1:3979,2)))
names(strain1strain2_allAnnotations)=c("strain1","strain2")


#Add columns for various annotations: Pwy, pcomplex, regulator_name, operon, kegg_modules


attrs=c("Pwy", "pcomplex", "regulator_name", "operon", "kegg_modules")


#function to construct X for each annotation set (contains dummy variables )
get_dummies=function(id_annot,nrow=3979){
    X=matrix(,ncol=0,nrow=nrow)
    
    #parallel lapply using multicores
    unique_annots=unique(id_annot[[2]])
    unique_annots=unique_annots[!is.na(unique_annots)]
    annot_dummy_vecs=mclapply(unique_annots,FUN=function(annot){
      annot_dummy=rep(FALSE,nrow)
      annot_true_index=as.numeric(id_annot[[1]][id_annot[[2]] %in% annot])
      annot_dummy[annot_true_index]=TRUE
      return(annot_dummy)
    },mc.cores=7)    
    
    
    X=Reduce(cbind,annot_dummy_vecs)
    return(X)
}


start.time = Sys.time()



for(i in 1:length(attrs)){
  
  attr=attrs[i]
  
  #construct X for each annotation set (contains dummy variables ) 
  id_annot=unique(id_allAttributes[,c("ids",attr)])
  X=get_dummies(id_annot)
  
  #X*t(X) to get coannotation similarity matrix
  id_sim_matrix=ifelse(X%*%t(X)>=1,1,0)
  
  #melt the coannotation similarity matrix
  result=melt_dist(as.dist(id_sim_matrix))
  
  cat(attr," is done\n")
  save(result,file=paste("Data/strain1strain2_annotCol_",attr,".RData",sep=""))
}

end.time = Sys.time()
end.time - start.time 


load("Data/strain1strain2_annotCol_Pwy.RData")
strain1strain2_allAnnotations=cbind(strain1strain2_allAnnotations,Pwy=result[[3]])
load("Data/strain1strain2_annotCol_pcomplex.RData")
strain1strain2_allAnnotations=cbind(strain1strain2_allAnnotations,pcomplex=result[[3]])
load("Data/strain1strain2_annotCol_operon.RData")
strain1strain2_allAnnotations=cbind(strain1strain2_allAnnotations,operon=result[[3]])
load("Data/strain1strain2_annotCol_regulator_name.RData")
strain1strain2_allAnnotations=cbind(strain1strain2_allAnnotations,regulator=result[[3]])
load("Data/strain1strain2_annotCol_kegg_modules.RData")
strain1strain2_allAnnotations=cbind(strain1strain2_allAnnotations,kegg_modules=result[[3]])




#type check
class(strain1strain2_allAnnotations)
apply(strain1strain2_allAnnotations,2,class)

save(strain1strain2_allAnnotations,file="Data/strain1strain2_allAnnotations.RData")