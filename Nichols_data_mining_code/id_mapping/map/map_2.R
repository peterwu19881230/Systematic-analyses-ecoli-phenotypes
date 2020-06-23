#Map the GO IDs to Nichols' strains by the .GAF file Suzi made (Note that the GO annotations here are different than if retrieved from UniProt)
#Map the UniProt IDs to Nichols' strains
#Map the pathways based on pathwayscol (a file from EcoCyc) to Nichols' strains


#Use ECK_1st_table created by id_mapping/1st_table.R 
load("Data/ECK_1st_table.RData")



#Map the GO IDs to Nichols' strains by the .GAF file Suzi made

##read.table doesn't work for .txt files that have NA values in some rows. I have to convert the file "2017_05_ECgene_association.ecocyc.txt" into .csv and then use read.csv

##Read the file
ECgene_association.ecocyc=read.csv("Data/2017_05_ECgene_association.ecocyc.csv",stringsAsFactors=F) 
##There are newer versions. Probably don't differ much



##Give the column names according to the teamwiki page: https://hexamer.tamu.edu/team/wiki/index.php/Making_GAF
colnames(ECgene_association.ecocyc)<-c("DB","DB Object ID","DB Object Symbol","Qualifier","GO ID","DB:Reference (IDB:Reference)","Evidence Code","With (or) From","Aspect","DB Object Name","DB Object Synonym (ISynonym)","DB Object Type","Taxon(ITaxon)","Date","Assigned By","Annotation Extension","Gene Product Form ID")


##map based on DB Object ID (EcoCyc symbol)
tab1=ECK_1st_table[,c("ids","EcoCycID")]
tab2=ECgene_association.ecocyc[,c("DB Object ID","GO ID")]
temp=inner_join(tab1,tab2,by=c("EcoCycID"="DB Object ID"))

names(temp)[3]="Data"
temp$'Data Type'="GO ID"
temp=temp[,c("ids","Data","Data Type")]
temp=temp[order(temp$ids %>% as.numeric),]
id.GO=temp


save(id.GO,file="Data/id.GO")





#Map the pathways based on pathwayscol (a file from EcoCyc) to Nichols' strains

pathwayscol=read.csv("Data/pathwayscol.csv",comment.char="#",stringsAsFactors=F)[-1,c(1:2,112:220)]
names(pathwayscol)=pathwayscol[1,]
pathwayscol=pathwayscol[-1,]
##1074 genes are in at least 1 pathway, others are not


##EcoCyc ID - Pathway ID - Pathway Name based on pathwayscol
EcoCycID.Pwys=data.frame()
for(i in 1:dim(pathwayscol)[1]){
  
  j=3 #Gene IDs start from the 3rd column
  while(pathwayscol[i,j]!=""){ ##empty cells were read as "" from read.csv()
    EcoCycID=pathwayscol[i,j]
    
    
    EcoCycID.Pwys=rbind(EcoCycID.Pwys,
                        c(EcoCycID,pathwayscol$`UNIQUE-ID`[i],pathwayscol$NAME[i]),
                        stringsAsFactors=F)
    
    j=j+1
    if(j>dim(pathwayscol)[2]) break
  }
  
}
colnames(EcoCycID.Pwys)=c("EcoCycID","Pathway ID","Pathway Name")

##Combine EcoCycID.Pwys with ids-ECK

###genes to pathways is a many to many relationship:
###length(unique(EcoCycID.Pwys$`EcoCyc ID`))  ###1074 unique genes
###length(EcoCycID.Pwys$`EcoCyc ID`) ###2681
###length(unique(EcoCycID.Pwys$`Pathway ID`)) ###424 unique pathways
###length(EcoCycID.Pwys$`Pathway ID`) ###2681


##map pathwayscol to Nichols id using EcoCyc ID
id_pathwayID=left_join(ECK_1st_table,EcoCycID.Pwys,by="EcoCycID")[,c("ids","Pathway ID")] %>% unique


save(id_pathwayID,file="Data/id_pathwayID.RData")



