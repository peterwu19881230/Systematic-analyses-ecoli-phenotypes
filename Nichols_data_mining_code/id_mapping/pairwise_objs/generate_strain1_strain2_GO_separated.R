#Pairwise semantic similarity of Nichols strains (by sets of GO annotations (BP, MF, CC are applied separately) )




id_allGO=id_allAttributes[,c("ids","GO")] %>% unique
id_allGO=id_allGO[!is.na(id_allGO$GO),]

id_BP=id_allGO[AnnotationDbi::Ontology(id_allGO$GO)=="BP" & !is.na(AnnotationDbi::Ontology(id_allGO$GO)),]
id_MF=id_allGO[AnnotationDbi::Ontology(id_allGO$GO)=="MF" & !is.na(AnnotationDbi::Ontology(id_allGO$GO)),]
id_CC=id_allGO[AnnotationDbi::Ontology(id_allGO$GO)=="CC" & !is.na(AnnotationDbi::Ontology(id_allGO$GO)),]



pairwise_GOannotation=function(id_GO){
  strain1strain2_transposed=as.data.frame((combn(1:3979,2))) 
  
  library(parallel)
  no_cores <- detectCores()-1
  cl <- makeCluster(no_cores)
  
  
  clusterExport(cl=cl,"id_GO",envir=environment()) 
  #envir=environment() allows clusterExport to search the environment within this function
  #ref: https://stackoverflow.com/questions/12023403/using-parlapply-and-clusterexport-inside-a-function/12024448
  
  
  GOannotation=parLapply(cl=cl,strain1strain2_transposed,fun=function(strain1strain2){ 
    ##Although input is not a list, it will be coerced by as.list() according to ?lapply
    
    strain1=as.character(strain1strain2[1])
    strain2=as.character(strain1strain2[2])
    
    annotation1=id_GO$GO[id_GO$ids==strain1]
    annotation2=id_GO$GO[id_GO$ids==strain2]
    
    return(list(annotation1,annotation2))
    
  })
  
  
  stopCluster(cl)
  return(GOannotation)
}




start.time=Sys.time()

pairwise_id_BP=pairwise_GOannotation(id_BP)
pairwise_id_MF=pairwise_GOannotation(id_MF)
pairwise_id_CC= pairwise_GOannotation(id_CC)

#Note: id without GO annotations have character(0) in the list
save(pairwise_id_BP,pairwise_id_MF,pairwise_id_CC,file="Data/pairwise_GOAnnotation_separated.RData")


end.time=Sys.time()
end.time-start.time #Time difference of 16.3797 mins





##parallelized code to get distances from Wang method 
#BiocManager::install("AnnotationHub", version = "3.8")
library(AnnotationHub)
library(GOSemSim)

library(parallel)
no_cores <- detectCores()-1
cl <- makeCluster(no_cores)


#Export packages to each cluster (https://stackoverflow.com/questions/33761123/r-how-can-i-export-methods-provided-by-a-package-to-a-psock-cluster)
clusterEvalQ(cl,c(library(AnnotationHub),library(GOSemSim),library(parallel))) # I tried to pass more than 1 library and it worked


ECK_GO_MF=godata('org.EcK12.eg.db',ont="MF",computeIC = F)
ECK_GO_BP=godata('org.EcK12.eg.db',ont="BP",computeIC = F)
ECK_GO_CC=godata('org.EcK12.eg.db',ont="CC",computeIC = F)

#Export required variables to each cluster (http://gforge.se/2015/02/how-to-go-parallel-in-r-basics-tips/)
clusterExport(cl=cl, c("ECK_GO_MF","ECK_GO_BP","ECK_GO_CC")) # I tried to pass more than 1 variable and it worked

start.time=Sys.time()

Wang_pairwise_similarity_BP=parLapply(cl=cl,X=pairwise_id_BP,fun=function(annots){
  mgoSim(annots[[1]],annots[[2]], semData=ECK_GO_BP, measure="Wang",combine="BMA")
}
)

Wang_pairwise_similarity_MF=parLapply(cl=cl,X=pairwise_id_MF,fun=function(annots){
  mgoSim(annots[[1]],annots[[2]], semData=ECK_GO_MF, measure="Wang",combine="BMA")
}
) 

Wang_pairwise_similarity_CC=parLapply(cl=cl,X=pairwise_id_CC,fun=function(annots){
  mgoSim(annots[[1]],annots[[2]], semData=ECK_GO_CC, measure="Wang",combine="BMA")
}
) 


save(Wang_pairwise_similarity_BP,Wang_pairwise_similarity_MF,Wang_pairwise_similarity_CC,file="Data/Wang_pairwise_similarity_separate.RData") #Time difference of 8.254868 mins

end.time=Sys.time()
end.time-start.time 
stopCluster(cl)



load("Data/strain1_strain2.RData")
strain1_strain2_WangBP_WangMF_WangCC=cbind(strain1_strain2,
                                           BP=unlist(Wang_pairwise_similarity_BP),
                                           MF=unlist(Wang_pairwise_similarity_MF),
                                           CC=unlist(Wang_pairwise_similarity_CC))


##test if they are correct:
##GO1=id_allAttributes[id_allAttributes$ids=="1","GO"] %>% unique
##GO1_BP=GO1[AnnotationDbi::Ontology(GO1)=="BP" & !is.na(AnnotationDbi::Ontology(GO1)=="BP")]

##GO3=id_allAttributes[id_allAttributes$ids=="3","GO"] %>% unique
##GO3_BP=GO3[AnnotationDbi::Ontology(GO3)=="BP" & !is.na(AnnotationDbi::Ontology(GO3)=="BP")]

##mgoSim(GO1_BP,GO3_BP,semData = ECK_GO_BP, measure="Wang",combine="BMA") #0.298


save(strain1_strain2_WangBP_WangMF_WangCC,file="Data/strain1_strain2_WangBP_WangMF_WangCC.RData")











