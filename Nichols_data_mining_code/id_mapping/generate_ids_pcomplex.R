#Goal: map protein complex to ids

currentDir=dirname(rstudioapi::getActiveDocumentContext()$path)
ptcomplex_table=read.table(paste(currentDir,"protcplxs.col",sep="/"),comment.char = "#",fill=T,sep="\t",quote="",header=T,check.names = F)
ptcomplex_table_small=ptcomplex_table[,c(1,grep("GENE-ID",names(ptcomplex_table),fixed=T))]

##use this strategy to transform to the new format: cbind(c("geneID 1","geneID 2"),"complex 1") => bind to yield the complete table
EcoCycID_pcomplex=data.frame()
for(i in 1:dim(ptcomplex_table_small)[1]){
  complex_id=ptcomplex_table_small[i,1]
  gene_id=ptcomplex_table_small[i,-1] %>% as.vector; gene_id=gene_id[gene_id!=""]
  
  
  #These 3 complex don't have any EcoCyc gene id in the table: NQOR-CPLX, TRANSENOYLCOARED-CPLX, CPLX0-7978. Therefore I write an exception to ignore them
  if(identical(gene_id,character(0))) next
  
  
  EcoCycID_pcomplex=rbind(EcoCycID_pcomplex,
                     cbind(gene_id,complex_id))
}

names(EcoCycID_pcomplex)[1]="EcoCycID"


#join id + EcoCyc id with EcoCyc id + protein complex identifier using: 1. ECK_1st_table 2. protcplxs.col
load("Data/ECK_1st_table.RData")


#Some ECK in Nichols doesn't map to EcoCyc IDs because genes.dat doesn't have them () -> If the newest genes.dat from EcoCyc 23.0 is used those EcoCyc IDs are mapped
id_EcoCycID=ECK_1st_table[,c("ids","EcoCycID")] %>% unique 

temp=inner_join(id_EcoCycID,EcoCycID_pcomplex,by="EcoCycID")


ids.pcomplex=temp[,c("ids","complex_id")]
names(ids.pcomplex)[2]="pcomplex"

save(ids.pcomplex,file="Data/ids.pcomplex.RData")

