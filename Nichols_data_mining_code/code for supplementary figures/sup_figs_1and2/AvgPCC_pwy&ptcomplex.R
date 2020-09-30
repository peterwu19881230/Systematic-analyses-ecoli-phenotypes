setwd("~/Dropbox/OMP\ shared/peter_wu/1st_paper_nichols_reanalysis/github/Systematic-analyses-ecoli-phenotypes/Nichols_data_mining_code")
source("Nichols_preload.R")

#Random expectation
random_expectation=abs(strain1strain2_allAnnotations_allSimilarities$pcc) #avg |PCC| of all pairs


#Function to generate df for the boxplot
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pre_process_for_ggplot=function(abs_cors_for_annotations,abs_cors_for_annotations_all){
  #have to turn no_gene_used into a factor in order to get proper color filling
  abs_cors_for_annotations$no_gene_used=as.factor(abs_cors_for_annotations$no_gene_used)
  abs_cors_for_annotations[[1]]=as.factor(abs_cors_for_annotations[[1]]) #suppose the annotation column is the 1st column
  
  #Lost the way I plot using different colors for genes in n gene pathways, but I think it's not important anymore
  df.abs_cors_for_annotations_all=data.frame(annotation=rep(names(abs_cors_for_annotations_all),sapply(abs_cors_for_annotations_all,length)),abs_pcc=unlist(abs_cors_for_annotations_all))
  names(df.abs_cors_for_annotations_all)[1]=names(abs_cors_for_annotations)[1]
  
  tab=left_join(df.abs_cors_for_annotations_all,abs_cors_for_annotations[,c(names(abs_cors_for_annotations)[1],"median_pcc","no_gene_used")],by=names(abs_cors_for_annotations)[1])
  tab=arrange(tab,as.numeric(no_gene_used),-median_pcc,-abs_pcc) #tab[order(a,-b,-c),] should also work
  ##Ref about using order(a,-1): https://stackoverflow.com/questions/7793295/how-to-order-a-data-frame-by-one-descending-and-one-ascending-column
  ##Note: I feel that the part about rev() is BSing. rev() doesn't really help (I wonder if the person really had run it)
  
  
  tab[[names(abs_cors_for_annotations)[1]]]=factor(tab[[names(abs_cors_for_annotations)[1]]],levels=unique(tab[[names(abs_cors_for_annotations)[1]]]))#Do this to prevent automatic x axis reordering
  
  ##ref: https://stackoverflow.com/questions/43877663/order-multiple-variables-in-ggplot2
  ##ref: http://www.sthda.com/english/wiki/ggplot2-facet-split-a-plot-into-a-matrix-of-panels
  ##Also, I am separating them into 2 lines
  
  return(tab)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



load("Data/cors_for_Pwy_Pcomplex.RData") #load the required data
text_size=10 #text size for x labels (pwy/pcomplex names)


random_line=geom_hline(yintercept = random_expectation,colour="red",linetype="dashed",size=0.5)

##~~~~ add this block to make the random line thicker~~~~
deviation=0.003
random_line_2=geom_hline(yintercept = random_expectation+deviation,colour="red",linetype="dashed",size=0.5)
random_line_3=geom_hline(yintercept = random_expectation-deviation,colour="red",linetype="dashed",size=0.5)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


load("Data/fdr_based_p_val_list.RData")

#get annotation data (pwy, complex)
load("Data/cors_for_Pwy_Pcomplex.RData") 
##There are multiple data in it. Only these will be used: abs_cors_for_pwys, abs_cors_for_pcomplexes



#Pwy

#no. of genes in pathways= c(3:22,27,28,31,41,47,48)
#I am currently using this: no matter how many genes are used, 
# pathways of same No. of genes defined by EcoCyc should be put into 1 box 
#=========================================================================================================

#Note: up to this point there are 282 pathways being used

pwy_tab=pre_process_for_ggplot(abs_cors_for_pwys,abs_cors_for_pwys_all)

##store a table to connect names to the numbers and replace names with number  
pwy_name_table=data.frame(no=as.factor(as.numeric(pwy_tab$Pwy)) ,pwy=pwy_tab$Pwy)
pwy_tab=cbind(no=pwy_name_table$no,pwy_tab)

tab=pwy_tab[,c(1,3,4,5)]


no_Pwy=pwy_tab[,c("no","Pwy")] %>% unique
pwy_fdr=fdr_based_p_val_list[[1]][match(as.character(no_Pwy$Pwy),abs_cors_for_pwys$Pwy)]

significance_numeric=pwy_fdr



p_Pwy_list=list()
count=1
for(no_of_gene in list(2,c(3,4),c(5,6),c(7,8,9,10),c(11,12,13,14,15,16,17,18,20,21,22,26,27,28,31,41,48))){
  
  
  no_no_gene_used=tab[tab$no_gene_used %in% no_of_gene, c("no","no_gene_used")] %>% unique
  
  
  one_star=which(significance_numeric<=0.05 & significance_numeric>0.01)
  one_star=one_star[one_star %in% as.numeric(no_no_gene_used$no)]
  
  two_star=which(significance_numeric<=0.01 & significance_numeric>0.001)
  two_star=two_star[two_star %in% as.numeric(no_no_gene_used$no)]
  
  three_star=which(significance_numeric<=0.001)
  three_star=three_star[three_star %in% as.numeric(no_no_gene_used$no)]
  
  star_list=list(one_star,two_star,three_star)
  label=c("*","**","***")
  
  
  annot_list=list()
  for(i in 1:3){
    
    if(identical(star_list[[i]],integer(0))){ #exception handler (sometimes there won't be any significance found)    
      annot_list[[i]]=NULL
      next
    } 
    
    
    label_df = data.frame(no = factor(star_list[[i]],levels=as.character(1:366)),
                          no_gene_used=no_no_gene_used$no_gene_used[which(as.numeric(no_no_gene_used$no) %in% star_list[[i]])],
                          abs_pcc = 1) 
    
    annot_list[[i]]=geom_text(data=label_df,label = label[i], angle = 90)
    
    
  }
  
  
  p_Pwy_list[[count]]=ggplot(tab[tab$no_gene_used %in% no_of_gene,], aes(x=no, y=abs_pcc)) +  
    theme_classic()+
    theme(plot.title= element_text(size = 20, hjust = 0.5),
          axis.text.x = element_text(angle = 90, hjust = 1,size = text_size,face="bold"),legend.position="none",
          axis.text.y=element_text(size=10),
          axis.title=element_text(size=20))+
    facet_grid(. ~ no_gene_used, scales="free_x", space="free_x")+
    geom_boxplot()+ylab("")+xlab("")+
    annot_list[[1]]+
    annot_list[[2]]+
    annot_list[[3]]+
    random_line+random_line_2+random_line_3+
    ylim(0,1)
  
  
  
  count=count+1
}





p_Pwy_correction_list=list()
count=1
for(no_of_gene in list(2,c(3,4),c(5,6),c(7,8,9,10),c(11,12,13,14,15,16,17,18,20,21,22,26,27,28,31,41,48))){
  
  
  no_no_gene_used=tab[tab$no_gene_used %in% no_of_gene, c("no","no_gene_used")] %>% unique
  
  
  one_star=which(significance_numeric<=0.05 & significance_numeric>0.01)
  one_star=one_star[one_star %in% as.numeric(no_no_gene_used$no)]
  
  two_star=which(significance_numeric<=0.01 & significance_numeric>0.001)
  two_star=two_star[two_star %in% as.numeric(no_no_gene_used$no)]
  
  three_star=which(significance_numeric<=0.001)
  three_star=three_star[three_star %in% as.numeric(no_no_gene_used$no)]
  
  star_list=list(one_star,two_star,three_star)
  label=c("*","**","***")
  
  
  annot_list=list()
  for(i in 1:3){
    
    if(identical(star_list[[i]],integer(0))){ #exception handler (sometimes there won't be any significance found)    
      annot_list[[i]]=NULL
      next
    } 
    
    
    label_df = data.frame(no = factor(star_list[[i]],levels=as.character(1:366)),
                          no_gene_used=no_no_gene_used$no_gene_used[which(as.numeric(no_no_gene_used$no) %in% star_list[[i]])],
                          abs_pcc = 1) 
    
    annot_list[[i]]=geom_text(data=label_df,label = label[i], angle = 90)
    
    
  }
  
  
  p_Pwy_correction_list[[count]]=ggplot(tab[tab$no_gene_used %in% no_of_gene,], aes(x=no, y=abs_pcc)) +  
    theme_classic()+
    theme(plot.title= element_text(size = 20, hjust = 0.5),
          axis.text.x = element_text(angle = 90, hjust = 1,size = text_size,face="bold"),legend.position="none",
          axis.text.y=element_text(size=10),
          axis.title=element_text(size=20))+
    facet_grid(. ~ no_gene_used, scales="free_x", space="free_x")+
    geom_point(color="black",alpha=0.5)+ylab("")+xlab("")+
    annot_list[[1]]+
    annot_list[[2]]+
    annot_list[[3]]+
    random_line+random_line_2+random_line_3+
    ylim(0,1) 
  
  count=count+1
}




#Pcomplex
#=========================================================================================================


pcomplex_tab=pre_process_for_ggplot(abs_cors_for_pcomplexes,abs_cors_for_pcomplexes_all)
#pcomplex_tab$no_gene_used %>% unique # 2  3  4  5  6  7  9  10 12 14 28 

##store a table to connect names to the numbers and replace names with number  
pcomplex_name_table=data.frame(no=as.factor(as.numeric(pcomplex_tab$pcomplex)) ,pcomplex=pcomplex_tab$pcomplex)
pcomplex_tab=cbind(no=pcomplex_name_table$no,pcomplex_tab)


tab=pcomplex_tab[,c(1,3,4,5)]

#abs_cors_for_pcomplexes


no_complex=pcomplex_tab[,c("no","pcomplex")] %>% unique
pcomplex_fdr=fdr_based_p_val_list[[2]][match(as.character(no_complex$pcomplex),abs_cors_for_pcomplexes$pcomplex)]

significance_numeric=pcomplex_fdr


p_Pcomplex_list=list()
count=1
for(no_of_gene in list(2,3,c(4,5,6,7,9,10,12,14,28))){
  
  no_no_gene_used=tab[tab$no_gene_used %in% no_of_gene, c("no","no_gene_used")] %>% unique
  
  
  one_star=which(significance_numeric<=0.05 & significance_numeric>0.01)
  one_star=one_star[one_star %in% as.numeric(no_no_gene_used$no)]
  
  two_star=which(significance_numeric<=0.01 & significance_numeric>0.001)
  two_star=two_star[two_star %in% as.numeric(no_no_gene_used$no)]
  
  three_star=which(significance_numeric<=0.001)
  three_star=three_star[three_star %in% as.numeric(no_no_gene_used$no)]
  
  star_list=list(one_star,two_star,three_star)
  label=c("*","**","***")
  
  
  annot_list=list()
  for(i in 1:3){
    
    if(identical(star_list[[i]],integer(0))){ #exception handler (sometimes there won't be any significance found)    
      annot_list[[i]]=NULL
      next
    } 
    
    
    label_df = data.frame(no = factor(star_list[[i]],levels=as.character(1:366)),
                          no_gene_used=no_no_gene_used$no_gene_used[which(as.numeric(no_no_gene_used$no) %in% star_list[[i]])],
                          abs_pcc = 1) 
    
    annot_list[[i]]=geom_text(data=label_df,label = label[i], angle = 90)
    
    
  }
  
  p_Pcomplex_list[[count]]=ggplot(tab[tab$no_gene_used %in% no_of_gene,], aes(x=no, y=abs_pcc)) +  
    theme_classic()+
    theme(plot.title= element_text(size = 20, hjust = 0.5),
          axis.text.x = element_text(angle = 90, hjust = 1,size = text_size,face="bold"),legend.position="none",
          axis.text.y=element_text(size=10),
          axis.title=element_text(size=20))+
    facet_grid(. ~ no_gene_used, scales="free_x", space="free_x")+
    geom_boxplot()+ylab("")+xlab("")+
    annot_list[[1]]+
    annot_list[[2]]+
    annot_list[[3]]+
    random_line+random_line_2+random_line_3+
    ylim(0,1)
  
  count=count+1
}



p_Pcomplex_correction_list=list()
count=1
for(no_of_gene in list(2,3,c(4,5,6,7,9,10,12,14,28))){
  
  no_no_gene_used=tab[tab$no_gene_used %in% no_of_gene, c("no","no_gene_used")] %>% unique
  
  
  one_star=which(significance_numeric<=0.05 & significance_numeric>0.01)
  one_star=one_star[one_star %in% as.numeric(no_no_gene_used$no)]
  
  two_star=which(significance_numeric<=0.01 & significance_numeric>0.001)
  two_star=two_star[two_star %in% as.numeric(no_no_gene_used$no)]
  
  three_star=which(significance_numeric<=0.001)
  three_star=three_star[three_star %in% as.numeric(no_no_gene_used$no)]
  
  star_list=list(one_star,two_star,three_star)
  label=c("*","**","***")
  
  
  annot_list=list()
  for(i in 1:3){
    
    if(identical(star_list[[i]],integer(0))){ #exception handler (sometimes there won't be any significance found)    
      annot_list[[i]]=NULL
      next
    } 
    
    
    label_df = data.frame(no = factor(star_list[[i]],levels=as.character(1:366)),
                          no_gene_used=no_no_gene_used$no_gene_used[which(as.numeric(no_no_gene_used$no) %in% star_list[[i]])],
                          abs_pcc = 1) 
    
    annot_list[[i]]=geom_text(data=label_df,label = label[i], angle = 90)
    
    
  }
  
  p_Pcomplex_correction_list[[count]]=ggplot(tab[tab$no_gene_used %in% no_of_gene,], aes(x=no, y=abs_pcc)) + 
    theme_classic()+
    theme(plot.title= element_text(size = 20, hjust = 0.5),
          axis.text.x = element_text(angle = 90, hjust = 1,size = text_size,face="bold"),legend.position="none",
          axis.text.y=element_text(size=10),
          axis.title=element_text(size=20))+
    facet_grid(. ~ no_gene_used, scales="free_x", space="free_x")+
    geom_point(color="black",alpha=0.5)+ylab("")+xlab("")+
    annot_list[[1]]+
    annot_list[[2]]+
    annot_list[[3]]+
    random_line+random_line_2+random_line_3+
    ylim(0,1)
  
  count=count+1
}


setwd("new_exps/permutation_based_p_val")


#save the figures as pdfs
ggsave("Pwy.pdf", arrangeGrob(grobs = p_Pwy_list,ncol=1),width=15,height=15)
ggsave("Pwy_correction.pdf", arrangeGrob(grobs = p_Pwy_correction_list,ncol=1),width=15,height=15)

ggsave("Pcomplex.pdf", arrangeGrob(grobs = p_Pcomplex_list,ncol=1),width=15,height=15)
ggsave("Pcomplex_correction.pdf", arrangeGrob(grobs = p_Pcomplex_correction_list,ncol=1),width=15,height=15)


#percentage of significant groups
mean(pwy_fdr>0.05)
mean(pcomplex_fdr>0.05)

