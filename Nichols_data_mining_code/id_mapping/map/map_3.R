#Create a master object for Nichols' ids


#The following objs are used
##ECK_1st_table #name and strain availability
load("Data/ECK_1st_table.RData")
apply(ECK_1st_table,2,class)
class(ECK_1st_table)


##GO term
load("Data/id.GO") 
used_id.GO=id.GO[,-3]
names(used_id.GO)[2]="GO"
apply(used_id.GO,2,class)
class(used_id.GO)


##UniProt ID 
load("Data/id.UniProt.RData") 
used_id.UniProt=id.UniProt
names(used_id.UniProt)[2]="UniProtID"
apply(used_id.UniProt,2,class)
class(used_id.UniProt)


##Pathway 
load("Data/id_pathwayID.RData")
names(id_pathwayID)[2]="Pwy"
apply(id_pathwayID,2,class)
class(id_pathwayID)


#protein complex (from proteinComplex.R) 
load("Data/ids.pcomplex.RData")
apply(ids.pcomplex,2,class)
class(ids.pcomplex)


#Regulon
load("Data/used_id.regulatorID.RData")

#Operon (from generate_id_operon.R)
load("Data/EcoCycGene_EcoCycOperon.RData")




#kegg modules (from scrapeKEGG.R)
load("Data/KEGGmodulesForNichols.RData")
used_id.KEGGmodules=data.frame()
for(i in 1:length(KEGGmodulesForNichols)){
  id=names(KEGGmodulesForNichols)[i]
  
  for(j in 1:length(KEGGmodulesForNichols[[i]])){
    if(is.na(KEGGmodulesForNichols[[i]][j])) break
    used_id.KEGGmodules=rbind(used_id.KEGGmodules,c(id,KEGGmodulesForNichols[[i]][j]))
  }
}
colnames(used_id.KEGGmodules)=c("ids","kegg_modules")
str(used_id.KEGGmodules)


#Map the ids from this paper: Genomewide landscape of gene–metabolome associations in Escherichia coli

# A column that tells whether the strain has any significant phenotypes (added 3/21/19)
TF=apply(Ternary_Data_324cutff_NAremoved,1,FUN = function(row){
  sum(row!=0)>0  
})
id.TF=data.frame(ids=as.character(1:3979), AnySigFitness=TF)


#Merge them altogether (!)Note: for the protein complex, the unwanted annotations for 
id_allAttributes=list(ECK_1st_table,used_id.GO,used_id.UniProt,id_pathwayID,ids.pcomplex,used_id.KEGGmodules,used_id.regulatorID,id.TF) %>%
  Reduce(function(df1,df2) left_join(df1,df2,by="ids"),.)

##operon uses EcoCycOD instead of the numeric id, so I am binding it separately here
id_allAttributes=left_join(id_allAttributes,EcoCycGene_EcoCycOperon,by="EcoCycID")


#Check if it is correct before saving 
apply(id_allAttributes,2,FUN=function(attr){
  return(c(class(attr),attr %>% unique %>% length)) #(!) class(attr) doesn't necessary return the right type and I don't know why
})

str(id_allAttributes)


#Add common names that are searchable on regulonDB website (ECK12 identifiers don't find stuff when searching on regulonDB)

#Regulators' name
generegulation_tmp=read.table("Data/regulon/generegulation_tmp.txt",sep="\t",quote="",stringsAsFactors = F) #no. of rows match the file: 4986-37+1=4950
ECK12_regulator=generegulation_tmp[,c(1,4)] %>% unique; names(ECK12_regulator)=c("regulator","regulator_name")
##Note: ECK12_regulator is a many-to-many-relationship table

ECK12_regulator_duplicated=ECK12_regulator[duplicated(ECK12_regulator$regulator),]



#join the regulator_name column from RegulonDB
temp=left_join(id_allAttributes,ECK12_regulator,by="regulator")



#rearrange some order for columns
temp=temp[,c("ids","ECK","JW","JW2","associated_gene_names","position","originalECKs","EcoCycID","bNumber","other.synonyms",
             "strain_availibility","GO","UniProtID","Pwy","pcomplex","kegg_modules","operon","regulator_name","AnySigFitness")] 


id_allAttributes=temp

#delete so the next time these 2 lazy loading objs can be regenerated
file1="Data/sourced/id_allAttributes.RData.rdb"
file2="Data/sourced/id_allAttributes.RData.rdx"
if(file.exists(file1) & file.exists(file2)){
  file.remove(file1)
  file.remove(file2)
} 



  
save(id_allAttributes,file="Data/sourced/id_allAttributes.RData")



