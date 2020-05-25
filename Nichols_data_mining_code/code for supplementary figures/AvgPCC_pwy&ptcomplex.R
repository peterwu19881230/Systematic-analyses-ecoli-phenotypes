#Random expectation
random_expectation=mean(1-strain1strain2_allAnnotations_allDistances$pcc) #avg |PCC| of all pairs

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



p_Pwy_list=list()
count=1
for(no_of_gene in list(2,c(3,4),c(5,6),c(7,8,9,10),c(11,12,13,14,15,16,17,18,20,21,22,26,27,28,31,41,48))){
  p_Pwy_list[[count]]=ggplot(tab[tab$no_gene_used %in% no_of_gene,], aes(x=no, y=abs_pcc)) +  
    theme_minimal()+
    #ggtitle("Average |PCC| for gene pairs in EcoCyc Pathways")+
    theme(plot.title= element_text(size = 20, hjust = 0.5),
          axis.text.x = element_text(angle = 90, hjust = 1,size = text_size),legend.position="none",
          axis.text.y=element_text(size=10),
          axis.title=element_text(size=20))+
    facet_grid(. ~ no_gene_used, scales="free_x", space="free_x")+
    geom_boxplot()+ylab("")+xlab("")+
    #ylab("|PCC|")+xlab("EcoCyc Pathway Annotations")+
    geom_hline(yintercept = random_expectation,linetype="longdash",colour="black",size=1,alpha=0.5)+
    ylim(0,1)
  ##Note: add fill=no_gene_used to aes() to make it colorful
  
  count=count+1
}



p_Pwy_correction_list=list()
count=1
for(no_of_gene in list(2,c(3,4),c(5,6),c(7,8,9,10),c(11,12,13,14,15,16,17,18,20,21,22,26,27,28,31,41,48))){
  p_Pwy_correction_list[[count]]=ggplot(tab[tab$no_gene_used %in% no_of_gene,], aes(x=no, y=abs_pcc)) +  
    theme_minimal()+
    #ggtitle("Average |PCC| for gene pairs in EcoCyc Pathways")+
    theme(plot.title= element_text(size = 20, hjust = 0.5),
          axis.text.x = element_text(angle = 90, hjust = 1,size = text_size),legend.position="none",
          axis.text.y=element_text(size=10),
          axis.title=element_text(size=20))+
    facet_grid(. ~ no_gene_used, scales="free_x", space="free_x")+
    geom_point(color="grey")+ylab("")+xlab("")+
    #ylab("|PCC|")+xlab("EcoCyc Pathway Annotations")+
    geom_hline(yintercept = random_expectation,linetype="longdash",colour="black",size=1,alpha=0.5)+
    ylim(0,1)
  ##Note: add colour=no_gene_used to aes and remove color="grey" in geom_point to make it colorful
  
  count=count+1
}




#Pcomplex
#=========================================================================================================



pcomplex_tab=pre_process_for_ggplot(abs_cors_for_pcomplexes,abs_cors_for_pcomplexes_all)
pcomplex_tab$no_gene_used %>% unique # 2  3  4  5  6  7  9  10 12 14 28 

##store a table to connect names to the numbers and replace names with number  
pcomplex_name_table=data.frame(no=as.factor(as.numeric(pcomplex_tab$pcomplex)) ,pcomplex=pcomplex_tab$pcomplex)
pcomplex_tab=cbind(no=pcomplex_name_table$no,pcomplex_tab)


tab=pcomplex_tab[,c(1,3,4,5)]




p_Pcomplex_list=list()
count=1
for(no_of_gene in list(2,3,c(4,5,6,7,9,10,12,14,28))){
  p_Pcomplex_list[[count]]=ggplot(tab[tab$no_gene_used %in% no_of_gene,], aes(x=no, y=abs_pcc)) +  
    theme_minimal()+
    #ggtitle("Average |PCC| for gene pairs in EcoCyc Pathways")+
    theme(plot.title= element_text(size = 20, hjust = 0.5),
          axis.text.x = element_text(angle = 90, hjust = 1,size = text_size),legend.position="none",
          axis.text.y=element_text(size=10),
          axis.title=element_text(size=20))+
    facet_grid(. ~ no_gene_used, scales="free_x", space="free_x")+
    geom_boxplot()+ylab("")+xlab("")+
    #ylab("|PCC|")+xlab("EcoCyc Pathway Annotations")+
    geom_hline(yintercept = random_expectation,linetype="longdash",colour="black",size=1,alpha=0.5)+
    ylim(0,1)
  ##Note: add fill=no_gene_used to aes() to make it colorful
  
  count=count+1
}



p_Pcomplex_correction_list=list()
count=1
for(no_of_gene in list(2,3,c(4,5,6,7,9,10,12,14,28))){
  p_Pcomplex_correction_list[[count]]=ggplot(tab[tab$no_gene_used %in% no_of_gene,], aes(x=no, y=abs_pcc)) +  
    theme_minimal()+
    #ggtitle("Average |PCC| for gene pairs in EcoCyc Pathways")+
    theme(plot.title= element_text(size = 20, hjust = 0.5),
          axis.text.x = element_text(angle = 90, hjust = 1,size = text_size),legend.position="none",
          axis.text.y=element_text(size=10),
          axis.title=element_text(size=20))+
    facet_grid(. ~ no_gene_used, scales="free_x", space="free_x")+
    geom_point(color="grey")+ylab("")+xlab("")+
    #ylab("|PCC|")+xlab("EcoCyc Pathway Annotations")+
    geom_hline(yintercept = random_expectation,linetype="longdash",colour="black",size=1,alpha=0.5)+
    ylim(0,1)
  ##Note: add colour=no_gene_used to aes and remove color="grey" in geom_point to make it colorful
  
  count=count+1
}




#save the figures as pdfs
ggsave("Pwy.pdf", arrangeGrob(grobs = p_Pwy_list,ncol=1),width=15,height=15)
ggsave("Pwy_correction.pdf", arrangeGrob(grobs = p_Pwy_correction_list,ncol=1),width=15,height=15)

ggsave("Pcomplex.pdf", arrangeGrob(grobs = p_Pcomplex_list,ncol=1),width=15,height=15)
ggsave("Pcomplex_correction.pdf", arrangeGrob(grobs = p_Pcomplex_correction_list,ncol=1),width=15,height=15)



#save the tables used to generate the figures
##write.csv(pcomplex_tab,file="pcomplex_tab.csv",quote=F, row.names=F)
##write.csv(pwy_tab,file="pwy_tab.csv",quote=F, row.names=F)


#add common names as additinal column and output the table: no. - pwy/pcomplex id - pwy/pcomplex common name

##save the tables: id_name
pwy_name_table=unique(pwy_name_table)
pcomplex_name_table=unique(pcomplex_name_table)

##write.csv(pwy_name_table,file="pwy_name_table.csv",quote=F, row.names=F)
##write.csv(pcomplex_name_table,file="pcomplex_name_table.csv",quote=F, row.names=F)


source("parse_pathways.dat.R") # get: pwyID_pwyName
source("parse_protcplxs.col.R") # get: pcomplexID_pcomplexName

pwy_name_table$pwy=as.character(pwy_name_table$pwy)
no_pwy_name=left_join(pwy_name_table,pwyID_pwyName,by=c("pwy"="pwyID"))[,c("no","pwy","pwyName")]
names(no_pwy_name)=c("no.","Pathway ID","Common name") #modify the colnames

pcomplex_name_table$pcomplex=as.character(pcomplex_name_table$pcomplex)
no_pcomplex_name=left_join(pcomplex_name_table,pcomplexID_pcomplexName,by=c("pcomplex"="pcomplexID"))[,c("no","pcomplex","pcomplexName")]
names(no_pcomplex_name)=c("no.","Protein complex ID","Common name") #modify the colnames

write.csv(no_pwy_name,file="no_pwy_name.csv",row.names = F)
write.csv(no_pcomplex_name,file="no_pcomplex_name.csv",row.names = F)  #??? Why are there monomers???








