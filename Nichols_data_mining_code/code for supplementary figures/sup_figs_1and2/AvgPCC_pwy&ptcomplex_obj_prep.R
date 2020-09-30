#Object preparation for AvgPCC_pwy&ptcomplex.R


#For Violin plot & distribution
#=========================================================================================================
#Part I: The following gathers all avg from the same genes in pwy

id_pwy=id_allAttributes[,c("ids","Pwy")] %>% unique
id_pwy=id_pwy[!is.na(id_pwy$Pwy),]
pathway_ids=attr_list(id_pwy[,c("Pwy","ids")])

samePwy_df=data.frame(id1=NULL,id2=NULL,pcc=NULL)
for(i in 1:length(pathway_ids)){
  if(length(pathway_ids[[i]])==1){ #There are 7 pwy that has only 1 gene
    print("The following only has 1 gene:")
    print(pathway_ids[i])
    next
  }
  comb=t(combn(pathway_ids[[i]],2))
  print(dim(comb)[1])
  
  for(j in 1:dim(comb)[1]){
    pair=data.frame(id1=comb[j,1],id2=comb[j,2],pcc=cor_strains[comb[j,1],comb[j,2]])
    samePwy_df=rbind(samePwy_df,pair)
  }
}

#Clean pairs that have been co-annotated in more than 1 pathway. I only want to count those pairs once.
sum(duplicated(samePwy_df))
samePwy_df=samePwy_df[!duplicated(samePwy_df),]


#Part II: Ptcomplexes
##I used inptcomplex_strain1strain2.pcc.samePTcomplex_Annot defined in another script
id_pcomplex=id_allAttributes[,c("ids","pcomplex")] %>% unique
id_pcomplex=id_pcomplex[!is.na(id_pcomplex$pcomplex),]


pcomplex_ids=attr_list(id_pcomplex[,c("pcomplex","ids")])

#The following gathers all avg from genes in the same protein complex
samePcomplex_df=data.frame(id1=NULL,id2=NULL,pcc=NULL)
for(i in 1:length(pcomplex_ids)){
  if(grepl(pattern = "MONOMER",x = names(pcomplex_ids)[i]) | length(pcomplex_ids[[i]])==1){ 
    ##If there is only 1 gene in Nichols' that's involved in the protein complex, remove it
    
    print("The following only has 1 gene:")
    print(pcomplex_ids[i])
    next
  }
  comb=t(combn(pcomplex_ids[[i]],2))
  print(dim(comb)[1])
  
  for(j in 1:dim(comb)[1]){
    pair=data.frame(id1=comb[j,1],id2=comb[j,2],pcc=cor_strains[comb[j,1],comb[j,2]])
    samePcomplex_df=rbind(samePcomplex_df,pair)
  }
}

#Clean pairs that have been co-annotated in more than 1 pathway. I only want to count those pairs once.
sum(duplicated(samePcomplex_df))
samePcomplex_df=samePcomplex_df[!duplicated(samePcomplex_df),]


#Pwy and/or Pcomplex 
samePwyPcomplex_df=rbind(samePwy_df,samePcomplex_df)
samePwyPcomplex_df=samePwyPcomplex_df[duplicated(samePwyPcomplex_df),]

samePwyORPcomplex_df=rbind(samePwy_df,samePcomplex_df)
samePwyORPcomplex_df=samePwyORPcomplex_df[!duplicated(samePwyORPcomplex_df),]




#The following gathers all avg from genes in pwy
cor_in_pwy=c()
for(i in 1:length(pathway_ids)){
  if(length(pathway_ids[[i]])==1){ #There are 7 pwy that has only 1 gene
    print("The following only has 1 gene:")
    print(pathway_ids[i])
    next
  }
  comb=t(combn(pathway_ids[[i]],2))
  print(dim(comb)[1])
  
  for(j in 1:dim(comb)[1]){
    cor_in_pwy=c(cor_in_pwy,cor_strains[comb[j,1],comb[j,2]])
  }
}


#The following gathers all avg from genes in ptcomplex
cor_in_pcomplex=c()
for(i in 1:length(pcomplex_ids)){
  if(length(pcomplex_ids[[i]])==1){ #There are 7 pwy that has only 1 gene
    print("The following only has 1 gene:")
    print(pcomplex_ids[i])
    next
  }
  comb=t(combn(pcomplex_ids[[i]],2))
  print(dim(comb)[1])
  
  for(j in 1:dim(comb)[1]){
    cor_in_pcomplex=c(cor_in_pcomplex,cor_strains[comb[j,1],comb[j,2]])
  }
}




#=========================================================================================================





#For boxplots
#=========================================================================================================
#Pwy (Curtis has done part of this)

##(!)In this experiment the pathway annotations shouldn't be limited to what was used in Nichols' fig.S1, so I use all of them
#Function to get tables for the following boxplots
#Dependency: my self-written function attr_list()

#annotationSet is any annotation set from id_allAttributes, "absolute" decides whther to use absolute value of pcc or not
get_boxplot_table=function(annotationSet,absolute=F){
  
  id_annotation=id_allAttributes[,c("ids",annotationSet)] %>% unique
  id_annotation=id_annotation[!is.na(id_annotation[[annotationSet]]),]
  annotation_ids=attr_list(cbind(id_annotation[[annotationSet]],id_annotation$ids)) 
  
  cors_for_annotationSet_all=list()
  cors_for_annotationSet=data.frame()
  count=1
  for(i in 1:length(annotation_ids)){
    if(length(annotation_ids[[i]]) ==1){ #Before I exclude 2 gene pwy/ptcomplex because Nichols' use >=3 gene pathways
      print("The following only has 1 gene found:")
      print(annotation_ids[i])
      next
    }
    comb=t(combn(annotation_ids[[i]],2))
    
    print(dim(comb)[1])
    
    pcc=c()
    for(j in 1:dim(comb)[1]){
      if(absolute==F){
        pcc[j]=cor_strains[comb[j,1],comb[j,2]] #extract pcc from cor_strains
      }else(
        pcc[j]=cor_strains[comb[j,1],comb[j,2]] %>% abs #extract pcc from cor_strains
      )
      
    }
    cors_for_annotationSet=rbind(cors_for_annotationSet,data.frame(annotation=names(annotation_ids)[i],
                                                                   avg_pcc=mean(pcc),
                                                                   median_pcc=median(pcc),
                                                                   std_pcc=sd(pcc),
                                                                   no_gene_used=length(annotation_ids[[i]])),
                                 stringsAsFactors=F)
    
    
    cors_for_annotationSet_all[[count]]=pcc
    names(cors_for_annotationSet_all)[count]=names(annotation_ids[i])
    count=count+1
  }
  
  names(cors_for_annotationSet)[1]=annotationSet #make the name of the annotation column = "annotationSet" argument (eg. "Pwy")
  cors_for_annotationSet[[annotationSet]]=as.character(cors_for_annotationSet[[annotationSet]]) #I don't know why this got converted into factor in rbind()
  
  
  
  
  
  #Return a list that has: 1. cors_for_annotationSet (dataframe) 2. cors_for_annotationSet_all (list)
  return(list(cors_for_annotationSet,cors_for_annotationSet_all)) 
  
}



##Pwy
temp=get_boxplot_table("Pwy") ##Do "non-absolute pcc"
cors_for_pwys=temp[[1]] 
cors_for_pwys_all=temp[[2]]

temp=get_boxplot_table("Pwy",absolute=T) ##Do absolute pcc
abs_cors_for_pwys=temp[[1]] 
abs_cors_for_pwys_all=temp[[2]] 



##Pcomplex
temp=get_boxplot_table("pcomplex") ##Do "non-absolute pcc"
cors_for_pcomplexes=temp[[1]] 
cors_for_pcomplexes_all=temp[[2]]

temp=get_boxplot_table("pcomplex",absolute=T) ##Do absolute pcc
abs_cors_for_pcomplexes=temp[[1]] 
abs_cors_for_pcomplexes_all=temp[[2]] 


rm(temp)

#=========================================================================================================


save(cors_for_pwys,cors_for_pwys_all,
     abs_cors_for_pwys,abs_cors_for_pwys_all,
     cors_for_pcomplexes,cors_for_pcomplexes_all,
     abs_cors_for_pcomplexes,abs_cors_for_pcomplexes_all,
     cor_in_pwy,cor_in_pcomplex,
     samePwy_df,samePcomplex_df,samePwyORPcomplex_df,samePwyPcomplex_df,file="Data/cors_for_Pwy_Pcomplex.RData")





