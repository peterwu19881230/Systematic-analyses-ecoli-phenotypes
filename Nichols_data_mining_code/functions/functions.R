#functions:
#This file is by default sourced


#Convert original data to ternary. Input is the original matrix and the threshold (cutoff) for scores. Input can be a matrix or a data frame. thresh should be a positive real value
TernaryConvert=function(matrix,thresh){
  ternary<-matrix(,nrow(matrix),ncol(matrix)) #This creates an empty matrix (Precisely, matrix with all NAs)
  rownames(ternary)<-rownames(matrix)
  colnames(ternary)<-colnames(matrix)
  ternary[matrix>-1*thresh&matrix<thresh]=0
  ternary[matrix<=-1*thresh]=-1
  ternary[matrix>=thresh]=1
  return(ternary) 
}


#Some narrower functions are defined as below (Previous ones are for more general use):

##Distance matrix based on pearson correlation coefficient. Correlations for rows will be used. This is what I thought Nasos has used: use="pairwise.complete.obs"(He hasn't replied for my last question)
##I use abs() because anti-correlation is better than not having any relationships
pcc_dist=function(matrix){
  ###pairwise distance would cause bias, but this is the best I can do now
  dissimilarity_pearson_pairwise<- 1 - abs(cor(t(matrix),use="pairwise.complete.obs", method="pearson"))
  distance_pcc <- as.dist(dissimilarity_pearson_pairwise)
  return(distance_pcc)
}

##Distance matrix based on Spearman correlation coefficient. Correlations for rows will be used. 
##Note: I used absolute value while in factoextra:get_dist the spearman based method doesn't
##Note: The reason why I used absolute value is the same reason as when using pcc
spearman_dist=function(matrix){
  dissimilarity_spearman_pairwise<- 1 - abs(cor(t(matrix),use="pairwise.complete.obs", method="spearman"))
  distance_spearman <- as.dist(dissimilarity_spearman_pairwise)
  return(distance_spearman)
}

##Distance matrix based on Euclidean distance
euclidean_dist=function(matrix){
  get_dist(matrix,method="euclidean")
}


#Distance matrix based on Mutual Information
my_mutual_info_dist=function(matrix){
  mi=infotheo::mutinformation(t(matrix) %>% as.data.frame) %>% natstobits #I have to convert the matrix to a dataframe. Otherwise, infotheo::mutinformation() will fail
  mi_distObj=as.dist(mi)
  mi_distance=1-mi_distObj
  return(mi_distance)
}


#This function takes a distance object as input and output a melted dataframe
#!(NA will be kicked out)
melt_dist=function(distance){
  m1=as.matrix(distance)
  m1[upper.tri(m1)]=NA #I suppose that no distance can be NA, so I can use this to do filtering
  diag(m1)=NA
  library(reshape2)
  m2=melt(m1) #I suppose that no distance can be NA, so I can use this to do filtering
  if(class(m2$Var1)!="character") m2$Var1=as.character(m2$Var1) #This is to prevent numeric names (the class should still be "character") being converted to "numeric" by melt()
  if(class(m2$Var2)!="character") m2$Var2=as.character(m2$Var2) #This is to prevent numeric names (the class should still be "character") being converted to "numeric" by melt()
  
  m2=m2[!is.na(m2[,3]) | is.nan(m2[,3]),] # "| is.nan(m2[,3])" is used because I want to keep the NaN values
  
  m2=m2[,c(2,1,3)] ##reorder the columns
  names(m2)=c("object_1","object_2","value") ##name the columns
  return(m2)
}


#This function takes a distance object as input and output a melted and sorted dataframe
#!(NA will be kicked out)
meltANDsort_dist=function(distance,decreasing=F){
  m2=melt_dist(distance) #melt_dist is a self-defined function
  m2=m2[order(m2[,3],decreasing=decreasing),] ##reorder by distance
  return(m2)
}


## Compute the attribute list by inputting relational data
##name is the name vecotr, attribute is the vector of corresponding attributes
attr_list=function(name,attribute){
  uniqueName=unique(name)
  
  attribute_list=list()
  for(i in 1:length(uniqueName)){
    attribute_list[[i]]=attribute[name==uniqueName[i]]
  }
  
  names(attribute_list)=uniqueName
  
  return(attribute_list)
}

#function to obtain confusion matrix and some metrics
#========================================================================================================================
confusionMatrix_metrics=function(cumSums,rankings=seq_along(cumSums),total=length(cumSums),seed=102){ #cumSums would be a numeric vector
  
  
  df=data.frame(
    TP=cumSums,
    FP=rankings-cumSums,
    TN=total-rankings-(max(cumSums)-cumSums),
    FN=max(cumSums)-cumSums,
    
    sensitivity=cumSums/max(cumSums),  # = TP/(TP+FN) #Note: this can also be viewed as a kind of coverage
    specificity = (total-rankings-(max(cumSums)-cumSums)) / (total-max(cumSums)), # = TN / (TN+FP)
    precision= cumSums / rankings, # = TP/(TP+FP)
    accuracy =(cumSums+total-rankings-max(cumSums)+cumSums)/total # = (TP + TN)/total     ##Equation looks complicated. Have to double check
    
  )
  
  randomCoAnnotation=rep(0,total); set.seed(seed); randomCoAnnotation[sample(total,max(cumSums))]=1
  randomCumSums=cumsum(randomCoAnnotation)
  
  random_df=data.frame(
    random_TP=randomCumSums,
    random_FP=rankings-randomCumSums,
    random_TN=total-rankings-(max(randomCumSums)-randomCumSums),
    random_FN=max(randomCumSums)-randomCumSums,
    
    random_sensitivity=randomCumSums/max(randomCumSums),  # = TP/(TP+FN)
    random_specificity = (total-rankings-(max(randomCumSums)-randomCumSums)) / (total-max(randomCumSums)), # = TN / (TN+FP)
    random_precision= randomCumSums / rankings, # = TP/(TP+FP)
    random_accuracy =(randomCumSums+total-rankings-max(randomCumSums)+randomCumSums)/total # = (TP + TN)/total     ##Equation looks complicated. Have to double check
    
  )
  
  
  return(cbind(df,random_df))
  
}
#========================================================================================================================

