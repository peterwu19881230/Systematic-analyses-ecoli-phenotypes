source("Nichols_preload.R")

auxotroph_1=read_excel("Data/No_growth_on_MOPS_Glu_0.4%.xlsx",sheet=1)
auxotroph_2=read_excel("Data/No_growth_on_MOPS_Glu_0.4%.xlsx",sheet=2)
auxotroph=dplyr::full_join(auxotroph_1,auxotroph_2,by=c("Gene name","EcoCyc ID","ECK number"))

names(auxotroph)=c("associated_gene_names","EcoCycID","bNumber","ECK")
used_strain_table=left_join(auxotroph,id_allAttributes[,
                                                       c("ids","associated_gene_names","EcoCycID","bNumber","ECK")] %>% unique,
                            by=names(auxotroph))


dim(used_strain_table) #hisA was duplicated (Fine. In Nichols there are duplicated strains for hisA)
not_used_id=used_strain_table$ids %>% unique
not_used_id=not_used_id[!is.na(not_used_id)]
not_used_id=as.numeric(not_used_id)

cor_matrix=cor(t(All_Data_NAimputed[setdiff(1:3979,not_used_id),])) 
cor_matrix[is.na(cor_matrix)]=-9 #change where cor=NA to -9 
dist_=as.dist(cor_matrix)
new_cor_table=meltANDsort_dist(dist_)
new_cor_table=new_cor_table[!(new_cor_table[,3]==-9),] 
new_cor_table$object_1=as.numeric(new_cor_table$object_1)
new_cor_table$object_2=as.numeric(new_cor_table$object_2)
new_cor_table[,3]=1-abs(new_cor_table[,3]) #change pcc to pcc based distance
names(new_cor_table)=c("strain1","strain2","pcc")

##calculate cumsum based on pcc

table_1=strain1strain2_allAnnotations_allDistances
table_2=left_join(new_cor_table,strain1strain2_allAnnotations_allDistances[,1:7],by=c("strain1","strain2"))


#function to get confusion matrix based on annot and similarity


get_confusionMatrix=function(df,annot,similarity,seed=9){   
  
  if(length(similarity)!=dim(df)[1]){ #if new similarity is in the the table, just retreive from the table
    similarity=df[,similarity]
  }
  df=cbind(df[,annot],similarity)
  
  if(length(annot)>=2){
    coannotation=(rowSums(df[,annot])>=1)
    #coannotation=(rowSums(df[,annot])>=length(annot)) this is A&B&...
    cumsum_=cumsum(coannotation[order(df[,dim(df)[2]])])
    return(confusionMatrix_metrics(cumsum_,seed=seed))
  } 
  
  cumsum_=cumsum(df[,1][order(df[,2])])
  return(confusionMatrix_metrics(cumsum_,seed=seed))
}



#function to plot results
graph_corr_annot=function(metric,similarity,samples,subset,cols,ylim,xlim,lwd,annot_list,annot_list_name){
  
  options(repr.plot.width = 10, repr.plot.height = 7)    
  
  for(similarity_ in similarity){
    
    x_lab=paste0("high similarity -- ranked pairs -- low similarity (",similarity_,")")
    for(metric_ in metric){
      
      random_metric=paste0("random_",metric_)
      for(i in seq(annot_list)){
        
        #precalculate and subset to prevent memory problem
        con1=get_confusionMatrix(df=table_1,annot_list[[i]],similarity_)[1:subset,]
        con2=get_confusionMatrix(df=table_2,annot_list[[i]],similarity_)[1:subset,]
        
        
        exp_list=list(con1,con2)
        exp_1_name=paste0("Same ",paste(annot_list_name[i],collapse=" "))
        exp_2_name=paste0("Same ",paste(annot_list_name[i],collapse=" ")," without auxotrophs")
        names(exp_list)=c(exp_1_name,exp_2_name)
        
        
        
        for(i in seq(exp_list)){
          if(i==1){
            plot(samples,exp_list[[i]][[metric_]][samples],xlab=x_lab,ylab=metric_,type='l',col=cols[i],
                 ylim=ylim,xlim=xlim,lwd = lwd, cex.lab=1.5) 
            grid(lty='solid')
            
            
          }else{
            lines(samples,exp_list[[i]][[metric_]][samples],col=cols[i],lwd = lwd) 
          }
          
          
        }
        
        legend(50, 1, legend=names(exp_list),
               col=cols, lty=1,lwd = lwd, cex=1.5, box.lty=0,bg="transparent") 
        
        #add the negative control
        lines(samples,exp_list[[1]][[random_metric]][samples],col='black',lty = 'dashed',lwd=2.5)
        
      }
    }
    
    
  } 
}    


metric="precision" #this is a highly imbalanced dataset so I will not use accuracy or specificity
similarity="pcc"
samples=1:500
subset=4000
cols=c("#56B4E9","#F3518A")
y=metric
ylim=c(0,1)
xlim=c(1,max(samples))
lwd=2
annot_list=list("pcomplex")
annot_list_name=c("protein complex(es)")


graph_corr_annot(metric,similarity,samples,subset,cols,ylim,xlim,lwd,annot_list,annot_list_name)
