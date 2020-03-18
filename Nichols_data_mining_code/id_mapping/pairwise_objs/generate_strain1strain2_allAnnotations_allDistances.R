#Create the obj: strain1 - strain 2 - (all annotations) - (all distances)

#Execute by: Rscript "path to Nichols_Data_Mining"
args = commandArgs(TRUE)
setwd(args)

source("Nichols_preload.R")

load("Data/strain1_strain2.RData")
load("Data/strain1strain2_allAnnotations.RData")
load("Data/strain1strain2_allDistances.RData")

#Check if strain1 - strain2 are all the same for strain1strain2_allAnnotations, strain1strain2_allDistances 
##identical(strain1strain2_allAnnotations[,1:2],strain1strain2_allDistances[,1:2])
##I don't know why the above doesn't work (may be due to some hidden properties of the 2 objs)
##(In both obj: 1. strain1, strain2 are both integer 2. Rownames are the same)

identical(as.integer(strain1strain2_allDistances$strain1),strain1strain2_allAnnotations$strain1)
identical(as.integer(strain1strain2_allDistances$strain2),strain1strain2_allAnnotations$strain2)


strain1strain2_allAnnotations_allDistances=cbind(strain1strain2_allAnnotations,strain1strain2_allDistances[,-(1:2)],stringsAsFactors=F)

#Create additional T,F column where T stands for: 2 strains both have at least 1 significant phenotype
dim(sigQuantitative_Data_324cutoff_condCollapsed)

TF=apply(sigQuantitative_Data_324cutoff_condCollapsed,1,FUN=function(row){
  sum(row!=0)>0
})

index=rownames(sigQuantitative_Data_324cutoff_condCollapsed)[TF] %>% as.integer

sigOrNot=( (strain1_strain2$strain1 %in% index) & (strain1_strain2$strain2 %in% index) ) #2438736

strain1strain2_allAnnotations_allDistances$sigPhenoInBoth=sigOrNot


#create a TF column: whether strain pairs exist in Price et al or not
load("Data/NichDeut_mi_ternary_allAnnot.RData") 

tab=cbind(NichDeut_mi_ternary_allAnnot[,c("strain1","strain2")],inPriceOrNot=TRUE)
temp=left_join(strain1strain2_allAnnotations_allDistances,tab,by=c("strain1","strain2"))
temp$inPriceOrNot[is.na(temp$inPriceOrNot)]=FALSE
strain1strain2_allAnnotations_allDistances=temp


#create a TF column: whether strain pairs exist in Fuhrer et al or not
gene_names=readxl::read_xls("Data/genomwide_metabolomic/sample_id_zscore.xls",sheet="Sheet1",col_names = F)

#str(gene_names)
#sum(duplicated(gene_names$associated_gene_names)) #no duplications in gene names
names(gene_names)="associated_gene_names"


#assign Nichols' id to Fuhrer's data (I can remap Fuhrer's id to all the annotation sets instead of using Nichols' id to help, but I am not sure it's worth it)
geneName_NicholsID=id_allAttributes[,c("associated_gene_names","ids")] %>% unique
#str(geneName_NicholsID)

df=left_join(gene_names,geneName_NicholsID,by="associated_gene_names")
#str(df) #increased no. of ID (3807 -> 3828) because of duplications in Nichols

df=df[!is.na(df$ids),]#remove NA (remove genes that are in Fuhrer's but not in Nichols)
#str(df) #3453 genes remaining

#Remove duplicates
df_clean=df[! (duplicated(df$associated_gene_names,fromLast = T) | duplicated(df$associated_gene_names,fromLast = F)), ]
#str(df_clean)

#See what those duplicates are
df_removed=df[ (duplicated(df$associated_gene_names,fromLast = T) | duplicated(df$associated_gene_names,fromLast = F)), ]
#str(df_removed) 
##there should only be 12*2 duplicated strains removed, but here it shows 42 
##-> There are genes with different original names but same gene names in Nichols -> For simplification, I think it's ok to remove all of them

strain1_strain2=strain1strain2_allAnnotations_allDistances[,c("strain1","strain2")]

FuhrerPairs=as.data.frame(t(combn(df_clean$ids %>% as.numeric,2)))
FuhrerPairs$inFuhrer=T
names(FuhrerPairs)[1:2]=c("strain1","strain2")

#dim(FuhrerPairs)[1] #5815755

strain1_strain2_inFuhrer=left_join(strain1_strain2,FuhrerPairs,by=c("strain1","strain2"))

swapedFuhrerPairs=FuhrerPairs[,c("strain2","strain1","inFuhrer")]
names(swapedFuhrerPairs)[1:2]=c("strain1","strain2")
strain1_strain2_inFuhrer=left_join(strain1_strain2_inFuhrer,swapedFuhrerPairs,by=c("strain1","strain2"))



inFuhrer=strain1_strain2_inFuhrer$inFuhrer.x | strain1_strain2_inFuhrer$inFuhrer.y

strain1_strain2_inFuhrer=cbind(strain1_strain2_inFuhrer[,c("strain1","strain2")],
                               inFuhrer)


#double check that all the intersected id pairs are matched
#sum(!is.na(strain1_strain2_inFuhrer$inFuhrer))==dim(FuhrerPairs)[1]
strain1_strain2_inFuhrer$inFuhrer[is.na(strain1_strain2_inFuhrer$inFuhrer)]=F


temp=left_join(strain1strain2_allAnnotations_allDistances,strain1_strain2_inFuhrer)
strain1strain2_allAnnotations_allDistances=temp

#double check that all the intersected id pairs are matched
#sum(strain1strain2_allAnnotations_allDistances$inFuhrer)==dim(FuhrerPairs)[1]


#delete so the next time these 2 lazy loading objs can be regenerated
file1="Data/sourced/strain1strain2_allAnnotations_allDistances.rdb"
file2="Data/sourced/strain1strain2_allAnnotations_allDistances.rdx"
if(file.exists(file1) & file.exists(file2)){
  file.remove(file1)
  file.remove(file2)
} 


save(strain1strain2_allAnnotations_allDistances,file="Data/sourced/strain1strain2_allAnnotations_allDistances.RData")









