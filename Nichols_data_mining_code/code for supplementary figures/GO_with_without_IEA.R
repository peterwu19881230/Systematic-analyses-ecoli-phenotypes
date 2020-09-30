#plot the GO violin plots after removing GO BP annotations with IEA code



##get id_BPnoIEA
#=========================================================================================================
##read the .gaf file
ECgene_association.ecocyc=read.csv("Data/2017_05_ECgene_association.ecocyc.csv",stringsAsFactors=F)

##Give the column names according to the teamwiki page: https://hexamer.tamu.edu/team/wiki/index.php/Making_GAF
colnames(ECgene_association.ecocyc)<-c("DB","DB Object ID","DB Object Symbol","Qualifier","GO ID","DB:Reference (IDB:Reference)","Evidence Code","With (or) From","Aspect","DB Object Name","DB Object Synonym (ISynonym)","DB Object Type","Taxon(ITaxon)","Date","Assigned By","Annotation Extension","Gene Product Form ID")

#remove IEAs
ECgene_association.ecocyc_noIEA=ECgene_association.ecocyc[ECgene_association.ecocyc[["Evidence Code"]]!="IEA",]

 
##map based on DB Object ID (EcoCyc symbol)
load("Data/ECK_1st_table.RData")
tab1=ECK_1st_table[,c("ids","EcoCycID")]
tab2=ECgene_association.ecocyc_noIEA[,c("DB Object ID","GO ID")]
temp=inner_join(tab1,tab2,by=c("EcoCycID"="DB Object ID"))

names(temp)[3]="Data"
temp$'Data Type'="GO ID"
temp=temp[,c("ids","Data","Data Type")]
temp=temp[order(temp$ids %>% as.numeric),]
id.GOnoIEA=temp[,1:2]
#=========================================================================================================



id_allGO=id.GOnoIEA %>% unique
names(id_allGO)=c("ids","GO")
id_allGO=id_allGO[!is.na(id_allGO$GO),]


id_BP=id_allGO[AnnotationDbi::Ontology(id_allGO$GO)=="BP" & !is.na(AnnotationDbi::Ontology(id_allGO$GO)),]

##calculate no. of annotations without IEA
#ECK_id_BP=left_join(id_BP,ECK_1st_table[,c("ECK","ids")] %>% unique,by="ids")
#ECK_BP=ECK_id_BP[,c("ECK","GO")]
#dim(ECK_BP)[1] # 4574 BP annotations (no IEA)


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

end.time=Sys.time()
end.time-start.time #Time difference of 17.57923 mins





##parallelized code to get distances from Wang method 
#BiocManager::install("AnnotationHub", version = "3.8")
library(AnnotationHub)
library(GOSemSim)

library(parallel)
no_cores <- detectCores()-1
cl <- makeCluster(no_cores)


#Export packages to each cluster (https://stackoverflow.com/questions/33761123/r-how-can-i-export-methods-provided-by-a-package-to-a-psock-cluster)
clusterEvalQ(cl,c(library(AnnotationHub),library(GOSemSim),library(parallel))) # I tried to pass more than 1 library and it worked


ECK_GO_BP=godata('org.EcK12.eg.db',ont="BP",computeIC = F)


#Export required variables to each cluster (http://gforge.se/2015/02/how-to-go-parallel-in-r-basics-tips/)
clusterExport(cl=cl, c("ECK_GO_BP")) # I tried to pass more than 1 variable and it worked

start.time=Sys.time()

Wang_pairwise_similarity_BP=parLapply(cl=cl,X=pairwise_id_BP,fun=function(annots){
  mgoSim(annots[[1]],annots[[2]], semData=ECK_GO_BP, measure="Wang",combine="BMA")
}
)


end.time=Sys.time()
end.time-start.time #Time difference of 10.89509 mins
stopCluster(cl)



load("Data/strain1_strain2.RData")
strain1_strain2_WangBPnoIEA=cbind(strain1_strain2,
                                  BP=unlist(Wang_pairwise_similarity_BP)
                                  )





save(strain1_strain2_WangBPnoIEA,file="Data/strain1_strain2_WangBPnoIEA.RData")
