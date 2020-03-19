#Random expectation
random_expectation=mean(1-strain1strain2_allAnnotations_allDistances$pcc)

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

#Pwy

#no. of genes in pathways= c(3:22,27,28,31,41,47,48)
#I am currently using this: no matter how many genes are used, 
# pathways of same No. of genes defined by EcoCyc should be put into 1 box 
#=========================================================================================================

#Note: up to this point there are 282 pathways being used
load("Data/cors_for_Pwy_Pcomplex.RData")
pwy_tab=pre_process_for_ggplot(abs_cors_for_pwys,abs_cors_for_pwys_all)
tab=pwy_tab


##for n=2,3,4,5
p_Pwy1=ggplot(tab[tab$no_gene_used %in% c(2,3,4,5),], aes(x=Pwy, y=abs_pcc)) +  
  theme_minimal()+
  #ggtitle("Average |PCC| for gene pairs in EcoCyc Pathways")+
  theme(plot.title= element_text(size = 20, hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1,size = 5),legend.position="none",
        axis.text.y=element_text(size=20),
        axis.title=element_text(size=20))+
  facet_grid(. ~ no_gene_used, scales="free_x", space="free_x")+
  geom_boxplot()+ylab("")+xlab("")+
  #ylab("|PCC|")+xlab("EcoCyc Pathway Annotations")+
  geom_hline(yintercept = random_expectation,linetype="longdash",colour="black",size=1,alpha=0.5)
##Note: add fill=no_gene_used to aes() to make it colorful


#I need this to correct n=2 into a dot plot
p_Pwy1_forCorrection_1=ggplot(tab[tab$no_gene_used %in% c(2,3,4,5),], aes(x=Pwy, y=abs_pcc))+ #I don't understand why "fill=no_gene_used" will give all black dots
  theme_minimal()+
  #ggtitle("Average |PCC| for gene pairs in EcoCyc Pathways")+
  theme(plot.title= element_text(size = 20, hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1,size = 5),legend.position="none",
        axis.text.y=element_text(size=20),
        axis.title=element_text(size=20))+
  facet_grid(. ~ no_gene_used, scales="free_x", space="free_x")+
  geom_point(color="grey")+ylab("")+xlab("")+
  #ylab("|PCC|")+xlab("EcoCyc Pathway Annotations")+
  geom_hline(yintercept = random_expectation,linetype="longdash",colour="black",size=1,alpha=0.5)
##Note: add colour=no_gene_used to aes and remove color="grey" in geom_point to make it colorful


#I need this to correct n=3 into a dot plot
p_Pwy1_forCorrection_2=ggplot(tab[tab$no_gene_used %in% c(2,3,4,5),], aes(x=Pwy, y=abs_pcc))+ #I don't understand why "fill=no_gene_used" will give all black dots
  theme_minimal()+
  #ggtitle("Average |PCC| for gene pairs in EcoCyc Pathways")+
  theme(plot.title= element_text(size = 20, hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1,size = 5),legend.position="none",
        axis.text.y=element_text(size=20),
        axis.title=element_text(size=20))+
  facet_grid(. ~ no_gene_used, scales="free_x", space="free_x")+
  geom_point(color="grey")+ylab("")+xlab("")+
  stat_summary(fun.y = median,fun.ymin=median,fun.ymax = median,color="black",
               geom = "crossbar", size = 0.3,width=1)+
  #ylab("|PCC|")+xlab("EcoCyc Pathway Annotations")+
  geom_hline(yintercept = random_expectation,linetype="longdash",colour="black",size=1,alpha=0.5)
##Note: add colour=no_gene_used to aes and remove color="grey" in geom_point to make it colorful



##for n=6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21, 22, 26, 27, 28, 31, 41, 48
p_Pwy2=ggplot(tab[tab$no_gene_used %in% c(6,7,8,9,10,11,12,13,14,15,16,17,18,20,21,22,26,27,28,31,41,48),], 
       aes(x=Pwy, y=abs_pcc)) +
  theme_minimal()+
  theme(plot.title= element_text(size = 20, hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1,size = 5),legend.position="none",
        axis.text.y=element_text(size=20),
        axis.title=element_text(size=20))+
  facet_grid(. ~ no_gene_used, scales="free_x", space="free_x")+
  geom_boxplot()+ylab("")+xlab("")+
  #ylab("|PCC|")+xlab("EcoCyc Pathway Annotations")+
  geom_hline(yintercept = random_expectation,linetype="longdash",colour="black",size=1,alpha=0.5)
##Note: add fill=no_gene_used to aes() to make it colorful


#=========================================================================================================



#Pcomplex
#=========================================================================================================

pcomplex_tab=pre_process_for_ggplot(abs_cors_for_pcomplexes,abs_cors_for_pcomplexes_all)
tab=pcomplex_tab

tab$no_gene_used %>% unique # 2  3  4  5  6  7  9  10 11 12 27 


#I need this to correct n=2 into a dot plot
p_Pcomplex_forCorrection_1=ggplot(tab[tab$no_gene_used %in% c(2,3),], aes(x=pcomplex, y=abs_pcc))+ #I don't understand why "fill=no_gene_used" will give all black dots
  theme_minimal()+
  #ggtitle("Average |PCC| for gene pairs in EcoCyc Protein Complex")+
  theme(plot.title= element_text(size = 20, hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1,size = 5),legend.position="none",
        axis.text.y=element_text(size=20),
        axis.title=element_text(size=20))+
  facet_grid(. ~ no_gene_used, scales="free_x", space="free_x")+
  geom_point(color="grey")+ylab("")+xlab("")+
  #ylab("|PCC|")+xlab("EcoCyc Protein Complex Annotations")+
  geom_hline(yintercept = random_expectation,linetype="longdash",colour="black",size=1,alpha=0.5)
##Note: add colour=no_gene_used to aes and remove color="grey" in geom_point to make it colorful


#I need this to correct n=3 into a dot plot
p_Pcomplex_forCorrection_2=ggplot(tab[tab$no_gene_used %in% c(2,3),], aes(x=pcomplex, y=abs_pcc))+ #I don't understand why "fill=no_gene_used" will give all black dots
  theme_minimal()+
  #ggtitle("Average |PCC| for gene pairs in EcoCyc Protein Complex")+
  theme(plot.title= element_text(size = 20, hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1,size = 5),legend.position="none",
        axis.text.y=element_text(size=20),
        axis.title=element_text(size=20))+
  facet_grid(. ~ no_gene_used, scales="free_x", space="free_x")+
  geom_point(color="grey")+ylab("")+xlab("")+
  stat_summary(fun.y = median,fun.ymin=median,fun.ymax = median,color="black",
               geom = "crossbar", size = 0.3,width=1)+
  #ylab("|PCC|")+xlab("EcoCyc Protein Complex Annotations")+
  geom_hline(yintercept = random_expectation,linetype="longdash",colour="black",size=1,alpha=0.5)
##Note: add colour=no_gene_used to aes and remove color="grey" in geom_point to make it colorful



p_Pcomplex=ggplot(tab[tab$no_gene_used %in% c(4, 5, 6,  7,  9, 10, 11, 12, 27),], aes(x=pcomplex, y=abs_pcc)) + 
  theme_minimal()+
  #ggtitle("Average |PCC| for gene pairs in EcoCyc Protein Complex")+
  theme(plot.title= element_text(size = 20, hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1,size = 5),legend.position="none",
        axis.text.y=element_text(size=20),
        axis.title=element_text(size=20))+
  facet_grid(. ~ no_gene_used, scales="free_x", space="free_x")+
  geom_boxplot()+ylab("")+xlab("")+
  #ylab("|PCC|")+xlab("EcoCyc Protein Complex Annotations")+
  geom_hline(yintercept = random_expectation,linetype="longdash",colour="black",size=1,alpha=0.5)
##Note: add fill=no_gene_used to aes() to make it colorful


#=========================================================================================================

#Save the graph in pdf

#Get the directory of current working script (Used when in Rstudio). Ref: https://stackoverflow.com/questions/47044068/get-the-path-of-current-script/47045368
'
dir_of_workingScript=dirname(rstudioapi::getSourceEditorContext()$path) 

pdf(file=paste(dir_of_workingScript,"/pwy_pcom_avgPCC.pdf",sep=""),width=15,height=5) 
p_Pwy1
p_Pwy1_forCorrection_1
p_Pwy1_forCorrection_2
p_Pwy2
p_Pcomplex_forCorrection_1
p_Pcomplex_forCorrection_2
p_Pcomplex
dev.off()
'





