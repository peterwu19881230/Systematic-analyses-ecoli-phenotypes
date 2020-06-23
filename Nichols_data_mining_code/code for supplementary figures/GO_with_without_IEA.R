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





save(strain1_strain2_WangBPnoIEA,file="strain1_strain2_WangBPnoIEA.RData")




###if the above .RData exists, start from this line:


#GO violin using the same similarity cutoffs

load("Data/strain1_strain2_WangBPnoIEA.RData")

#remove pairs that don't have semantic similarity; change class of strains from chr to numeric
strain1_strain2_WangBPnoIEA=strain1_strain2_WangBPnoIEA[!is.na(strain1_strain2_WangBPnoIEA$BP),]
strain1_strain2_WangBPnoIEA$strain1=as.numeric(strain1_strain2_WangBPnoIEA$strain1)
strain1_strain2_WangBPnoIEA$strain2=as.numeric(strain1_strain2_WangBPnoIEA$strain2)

strain1strain2_allAnnotations_allDistances_WangBPnoIEA=left_join(strain1_strain2_WangBPnoIEA,strain1strain2_allAnnotations_allDistances,by=c("strain1","strain2"))
names(strain1strain2_allAnnotations_allDistances_WangBPnoIEA)[3]="BPnoNA"


#Violin using GO annotations


options(repr.plot.width = 15, repr.plot.height = 15/2)

#get p5+p6
xlabs=c("All","|PCC|>0.75","Mutual information>0.15","Mutual information>0.32")

##cutoff index for pairs that have pcc>0.75
cuoff_index_pcc=sum((1-strain1strain2_allAnnotations_allDistances_WangBPnoIEA$pcc)>0.75)
cuoff_index_mi_1=sum((1-strain1strain2_allAnnotations_allDistances_WangBPnoIEA$mi_ternary)>0.15)
cuoff_index_mi_2=sum((1-strain1strain2_allAnnotations_allDistances_WangBPnoIEA$mi_ternary_collapsedCond)>0.32)
##note: no need to sort because I used sum()


df=strain1strain2_allAnnotations_allDistances_WangBPnoIEA


##the first box
#------------------------------------------------
box1=df$BPnoNA[!is.na(df$BPnoNA)]
#------------------------------------------------

##the second box
#------------------------------------------------
similarity="pcc"
box2=df$BPnoNA[order(df[[similarity]])][1:cuoff_index_pcc]
box2=box2[!is.na(box2)]
#------------------------------------------------

##the third box
#------------------------------------------------
similarity="mi_ternary"
box3=df$BPnoNA[order(df[[similarity]])][1:cuoff_index_mi_1]
box3=box3[!is.na(box3)]
#------------------------------------------------

##the forth box
#------------------------------------------------
similarity="mi_ternary_collapsedCond"
box4=df$BPnoNA[order(df[[similarity]])][1:cuoff_index_mi_2]
box4=box4[!is.na(box4)]
#------------------------------------------------


df_all=data.frame(All=box1)
df1=data.frame(aboveCutoff1=box2)
df2=data.frame(aboveCutoff2=box3)
df3=data.frame(aboveCutoff3=box4)


p7=ggplot() +
  geom_violin(data = df_all,aes(xlabs[1],All)) +
  geom_boxplot(data = df_all,aes(xlabs[1],All),width=0.1,outlier.shape = NA)+ #outlier.shape decides the shape of outliers. Here I don't let them show
  geom_violin(data = df1,aes(xlabs[2],aboveCutoff1)) +
  geom_boxplot(data = df1,aes(xlabs[2],aboveCutoff1),width=0.1,outlier.shape = NA)+
  geom_violin(data = df2,aes(xlabs[3],aboveCutoff2)) +
  geom_boxplot(data = df2,aes(xlabs[3],aboveCutoff2),width=0.1,outlier.shape = NA)+
  geom_violin(data = df3,aes(xlabs[4],aboveCutoff3)) +
  geom_boxplot(data = df3,aes(xlabs[4],aboveCutoff3),width=0.1,outlier.shape = NA)+
  scale_x_discrete("",limits=xlabs)+ 
  #I want the x axis to be empty. And if I don't use this, the order is not right
  scale_y_continuous("Semantic similarity")+
  theme(text=element_text(size=25), 
        axis.text.y=element_text(size=25),
        axis.title=element_text(size=25))

p7+
  annotate(geom="text", x=2.2, y=1, label="***",color="black",size=10)+
  annotate(geom="text", x=3.2, y=1, label="***",color="black",size=10)+
  annotate(geom="text", x=4.2, y=1, label="***",color="black",size=10)



#significance test
wilcox.test(box2,box1,alternative="greater") #p-value < 2.2e-16
wilcox.test(box3,box1,alternative="greater") #p-value < 2.2e-16
wilcox.test(box4,box1,alternative="greater") #p-value < 2.2e-16