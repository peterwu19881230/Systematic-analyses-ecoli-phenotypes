#Compute different types of Nichols' phenotype data


#Primary data preparation
All_Data=read.csv("Data/allData.csv")[,-1]
Ternary_Data_NAnotimputed=TernaryConvert(matrix=All_Data,thresh=3.463)
#save(All_Data,Ternary_Data_NAnotimputed,file="Data/sourced/PrimaryDataSets.RData")



## Impute all NA with mean of all scores
avg=mean(as.numeric(as.matrix(All_Data)),na.rm=T)
All_Data_NAimputed=All_Data
indexForNA=which(is.na(All_Data_NAimputed),arr.ind=T)
for(i in 1:dim(indexForNA)[1]){
  row=indexForNA[i,1]
  col=indexForNA[i,2]
  All_Data_NAimputed[row,col]=avg
}
which(is.na(All_Data_NAimputed),arr.ind=T)##Verify that there is no NA anymore
#save(All_Data_NAimputed,file="Data/sourced/All_Data_NAimputed.RData")


##Convert All_Data_NAimputed to Ternary
matrix<-matrix(,3979,324)
rownames(matrix)<-rownames(All_Data_NAimputed)
colnames(matrix)<-colnames(All_Data_NAimputed)
matrix[All_Data_NAimputed>-3.463&All_Data_NAimputed<3.463]<-0
matrix[All_Data_NAimputed<=-3.463]<--1
matrix[All_Data_NAimputed>=3.463]<-1
Ternary_Data=matrix
#save(Ternary_Data,file="Data/sourced/Ternary_Data")

#Ternary data using 324 cutoffs for each condition
dat=All_Data
Ternary_Data_324cutff=matrix(numeric(3979*324),nrow=3979)
for(i in 1:324){
  Ternary_Data_324cutff[,i]=TernaryConvert(dat[,i] %>% as.matrix,thresh =min_forEachCond[i] )
}


rownames(Ternary_Data_324cutff)=1:3979
#save(Ternary_Data_324cutff,file="Data/sourced/Ternary_Data_324cutff.RData")


#Impute the NA values with 0. Ref: https://gist.github.com/Jfortin1/d4888d68359a36fbda60
impute_matrix=function(matrix){
  missing=which(is.na(matrix), arr.ind=T)
  if (length(missing)!=0){
    for (j in 1:nrow(missing)){
      matrix[missing[j,1],missing[j,2]]=0
    }
  }
  matrix
}

Ternary_Data_324cutff_NAremoved=impute_matrix(Ternary_Data_324cutff)
rownames(Ternary_Data_324cutff_NAremoved)=1:3979
#save(Ternary_Data_324cutff_NAremoved,file="Data/sourced/Ternary_Data_324cutff_NAremoved.RData")


##Add various ternary data based on differnt FDR: 1%, 10%, 15%, 20%
##========================================================================

convertToTernary_byFDR=function(dat=All_Data,FDR_cond=FDR_cond,FDR=0.05){ #FDR_cond is defined in phenotypeCutoff.R
  
  #determine the cutoff for each condition
  min_forEachCond=sapply(FDR_cond,FUN=function(q_fdr){
    q_fdr$score[q_fdr$FDR<=FDR] %>% abs %>% min   
  })
  
  #generate the ternary data by using the cutoffs to filter
  Ternary_Data_324cutff=matrix(numeric(3979*324),nrow=3979)
  for(i in 1:324){
    Ternary_Data_324cutff[,i]=TernaryConvert(dat[,i] %>% as.matrix,thresh =min_forEachCond[i] )
  }
  
  rownames(Ternary_Data_324cutff)=1:3979
  
  #impute NAs with 0. Ref: https://gist.github.com/Jfortin1/d4888d68359a36fbda60
  matrix=Ternary_Data_324cutff
  
  missing=which(is.na(matrix), arr.ind=T)
  if (length(missing)!=0){
    for (j in 1:nrow(missing)){
      matrix[missing[j,1],missing[j,2]]=0
    }
  }
  
  Ternary_Data_324cutff_NAremoved=matrix
  
  return(Ternary_Data_324cutff_NAremoved)
}

Ternary_Data_324cutff_NAremoved_1_FDR=convertToTernary_byFDR(dat=All_Data,FDR_cond=FDR_cond,FDR=0.01)
Ternary_Data_324cutff_NAremoved_10_FDR=convertToTernary_byFDR(dat=All_Data,FDR_cond=FDR_cond,FDR=0.1)
Ternary_Data_324cutff_NAremoved_15_FDR=convertToTernary_byFDR(dat=All_Data,FDR_cond=FDR_cond,FDR=0.15)
Ternary_Data_324cutff_NAremoved_20_FDR=convertToTernary_byFDR(dat=All_Data,FDR_cond=FDR_cond,FDR=0.20)

##save(Ternary_Data_324cutff_NAremoved_1_FDR,
##     Ternary_Data_324cutff_NAremoved_10_FDR,
##     Ternary_Data_324cutff_NAremoved_15_FDR,
##     Ternary_Data_324cutff_NAremoved_20_FDR,
##     file="Data/ternary_data_various_FDR.RData")

##========================================================================




#Collapse 324 Conditions in Nichols' by the 114 treatments, and create necessary objs for further CorrVSAnnot analysis

#Criteria: if there is either 1 or -1, use that. If both are found or there are only 0s + NAs , switch to 0

dat=Ternary_Data_NAnotimputed

Ternary_Data_condCollapsed=matrix(,nrow=3979,ncol=length(uniqueChemIndex))
rownames(Ternary_Data_condCollapsed)=1:3979

for(i in 1:length(uniqueChemIndex)){
  
  conds=dat[,uniqueChemIndex[[i]]] %>% as.data.frame #as.data.frame is to prevent treatment with only 1 condition from crashing
  
  treatment=apply(conds,1,FUN = function(row){
    if( (1 %in% row) && !(-1 %in% row) ) return(1) 
    if( (-1 %in% row) && !(1 %in% row) ) return(-1) 
    return(0) # If neither of the above critera were met, use 0 
    
  })
  
  Ternary_Data_condCollapsed[,i]=treatment
  
}


#Verify the authenticity
complete(Ternary_Data_condCollapsed) # Note: this is a self-defined function
#save(Ternary_Data_condCollapsed,file="Data/sourced/Ternary_Data_condCollapsed.RData")




##Collapse 324 Conditions in Nichols' by the 114 treatments


##Criteria: if there is either 1 or -1, use that. If both are found or there are only 0s + NAs , switch to 0

dat=Ternary_Data_324cutff

Ternary_Data_324cutff_condCollapsed=matrix(,nrow=3979,ncol=length(uniqueChemIndex))
rownames(Ternary_Data_324cutff_condCollapsed)=1:3979

for(i in 1:length(uniqueChemIndex)){
  
  conds=dat[,uniqueChemIndex[[i]]] %>% as.data.frame #as.data.frame is to prevent treatment with only 1 condition from crashing
  
  treatment=apply(conds,1,FUN = function(row){
    if( (1 %in% row) && !(-1 %in% row) ) return(1) 
    if( (-1 %in% row) && !(1 %in% row) ) return(-1) 
    return(0) # If neither of the above critera were met, use 0 
    
  })
  
  Ternary_Data_324cutff_condCollapsed[,i]=treatment
  
}

rownames()

#Verify the authenticity
complete(Ternary_Data_324cutff_condCollapsed) # Note: this is a self-defined function


#save(Ternary_Data_324cutff_condCollapsed,file="Data/sourced/Ternary_Data_324cutff_condCollapsed.RData")



#Remove strains that have no phenotypes (for quantitative and qualitative)
TF=apply(Ternary_Data_324cutff_NAremoved,1,FUN = function(row){
  sum(row!=0)>0  
})



#Quantitative phenotypic profiles after removing non-significant phenotypes
significantORnot_mat=ifelse(Ternary_Data_324cutff_NAremoved!=0,1,0) #Note: ifelse can take a matrix or a dataframe and preserve the structure (ref: ifelse_matrixORdataframe.R)
sigQuantitative_Data_324cutoff= All_Data_NAimputed * significantORnot_mat #pairwise multiplication of a matrix and a numeric dataframe (I have confirmed that this works)
#save(sigQuantitative_Data_324cutoff,file="Data/sourced/sigQuantitative_Data_324cutoff.RData") 


#Quantitative phenotypic profiles after collapsing conditions (most significant phenotypes are used among different concentrations of conditions and if there's a conflict,there will be no phenotype (fitness = 0) )
sigQuantitative_Data_324cutoff_condCollapsed=as.data.frame(matrix(0,nrow=3979,ncol=length(uniqueChemIndex)))
for(i in seq_along(uniqueChemIndex)){
  
  sameChemical_mat=sigQuantitative_Data_324cutoff[,uniqueChemIndex[[i]]] %>% as.matrix
  
  sigQuantitative_Data_324cutoff_condCollapsed[,i]=apply(sameChemical_mat,1,FUN = function(Row){
    
    if(max(Row)>0 & min(Row)<0){ #If there are conflict phenotypes just set it to fitness=0
      fitness=0
    }else{
      
      index=ifelse(abs(Row)==max(abs(Row)),T,F)
      
      fitness=Row[index]
      
      if(length(fitness)!=1){ #if there are 2 identical fitness values (very unlikely), I will just take the first 1
        fitness=fitness[1]
      } 
      
    }
  
    return(fitness)  
  } 
  )
  
}  

names(sigQuantitative_Data_324cutoff_condCollapsed)=names(uniqueChemIndex)

#save(sigQuantitative_Data_324cutoff_condCollapsed,file="Data/sourced/sigQuantitative_Data_324cutoff_condCollapsed.RData")   











