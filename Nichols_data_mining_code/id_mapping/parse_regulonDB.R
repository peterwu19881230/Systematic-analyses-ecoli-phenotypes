#Parse files from regulonDB
# Step I: ECK_1st_table + regulon/object_synonym.txt to get long ECKs using b number
# Step II: For the genes not yet mapped with the long ECK, use ECK_1st_table + regulon/object_synonym.txt get long ECKs by ECK number
# Step III: For the genes not yet mapped with the long ECK, use ECK_1st_table + gene.txt to get long ECKs using gene names
# Step IV: manually resove the remaining ones (still 8 unresolved ones)
# Step V: Get the regulon annotations 
# Step VI: Get the operon annotations

##Note: I wasn't able to find clear descriptions for each file, but I suspect that: 
##1. generegulation_tmp.txt contains the details for each gene being regulated
##2. regulon_d_tmp.txt and regulon_tmp contains all regulon ids but regulon_d_tmp.txt has 2 additional columns: RE_ID, REGULON_TF_GROUP
##3. gene.txt is the file that has long ECKs connected to gene names
##4. object_synonym.txt is the file that links gene identifiers (eg. ECK120000045) in generegulation_tmp.txt to what's in Nichols'


##Note: I am not sure which strain(s) of E. coli they used (on the home page it says K-12)


object_synonym=read.table("Data/regulon/object_synonym.txt",sep="\t",quote="")[,1:2] 
##no. of rows match the file: 27757-32+1=27726 




# Step I: ECK_1st_table + regulon/object_synonym.txt to get long ECKs using b number

#Mapping gene ids: Use b number first (object_synonym doesn't have EcoCyc ID)

names(object_synonym)=c("longECK_by_regulonDB","bNumber") 
##3rd column aren't all bNumbers. However, merge() requires 2 columns with identical name
merged_ECK_1st_table=merge(ECK_1st_table,object_synonym,by="bNumber",all.x = T,all.y = F)

##Problem: no. of rows increases from 3986 to 3991. 5 additional rows: some bNumbers got more than 1 long ECK:
##print duplicated rows including the original ones
merged_ECK_1st_table[duplicated(merged_ECK_1st_table[,-13]) | duplicated(merged_ECK_1st_table[,-13],fromLast = T),] 
##Ref: https://stackoverflow.com/questions/7854433/finding-all-duplicate-rows-including-elements-with-smaller-subscripts


##Find out which ones to delete:
##ECK120008976 => an operon (Found from object_synonym.txt. They use the same bNumber: b1826, for both the gene and the operon). Also, no results found from gene.txt!
##ECK120007130 => No results found from gene.txt!
##ECK120015146 => No results found from gene.txt!

corrected_merged_ECK_1st_table=merged_ECK_1st_table[!as.character(merged_ECK_1st_table$longECK_by_regulonDB) %in% 
                                                      c("ECK120008976","ECK120007130","ECK120015146"),]


corrected_merged_ECK_1st_table$longECK_by_regulonDB %>% is.na %>% sum #86 of them don't have a long ECK by_regulonDB
ECK_1st_table$bNumber %>% is.na %>% sum #86 
corrected_merged_ECK_1st_table$bNumber %>% is.na %>% sum #86  
##These guys should also not find a long ECK => at least 86 of Nichols' should not have a long ECK => the 86 found above 
##(86=86, meaning every Nichols' strain that has a b number has a long ECK)






# Step II: For the genes not yet mapped with the long ECK, use ECK_1st_table + regulon/object_synonym.txt to get long ECKs by ECK number

#Use ECK to map long ECK to those 86 genes 
copy_object_synonym=object_synonym
names(copy_object_synonym)=c("longECK_by_regulonDB","ECK")
temp=corrected_merged_ECK_1st_table[is.na(corrected_merged_ECK_1st_table$bNumber),-13]
temp2=merge(temp,copy_object_synonym,by="ECK",all.x=T,all.y=F) #Still lots of them are missing => have to reluctantly use gene names







# Step III: For the genes not yet mapped with the long ECK, use ECK_1st_table + gene.txt to get long ECKs using gene names

#For the above 86 genes, I have to reluctantly use gene names to map because there are no other identifiers found for these genes
#Luckily, these genes don't have gene name synonyms (NA in other.synonyms column that I mapped before)
gene=read.table("Data/regulon/gene.txt",sep="\t",quote="",fill=T)[,1:2] #fill=T => there is at least 1 row that has unequal row length
##no. of rows match the file: 4718-40+1=4679
names(gene)=c("longECK_by_regulonDB","associated_gene_names")
gene$associated_gene_names=as.character(gene$associated_gene_names)
gene$longECK_by_regulonDB=as.character(gene$longECK_by_regulonDB)


##Are all the gene names unique in gene.txt? =>No
gene$associated_gene_names %>% duplicated %>% sum #26

##What are the duplicates?
gene$associated_gene_names[duplicated(gene$associated_gene_names) | duplicated(gene$associated_gene_names,fromLast=T)]
###sokE is duplicated once (and it's not used in Nichols'). Others are blank (no gene names)
###=> These suggests that I can still use gene names in gene.txt to map long ECK to Nichols'


corrected_merged_ECK_1st_table$associated_gene_names=as.character(corrected_merged_ECK_1st_table$associated_gene_names) 
##factor -> character to prevent merge() from not working properly
temp3=merge(temp,gene,by="associated_gene_names",all.x=T,all.y=F)


temp4=merge(temp2,temp3,by=intersect(names(temp2),names(temp3))[-13],all.x=T,all.y=T)
temp4$longECK_by_regulonDB.x==temp4$longECK_by_regulonDB 
##This is to confirm that there is no conflict between "mapping by ECK" and "mapping by gene names"


long=sapply(1:dim(temp4)[1],FUN=function(i){
  if(!is.na(as.character(temp4$longECK_by_regulonDB.x[i]))){
    return(as.character(temp4$longECK_by_regulonDB.x[i]))
  }
  
  if(!is.na(as.character(temp4$longECK_by_regulonDB.y[i]))){
    return(as.character(temp4$longECK_by_regulonDB.y[i]))
  }
  
  return(NA)
  
})

temp5=cbind(temp4[,-c(13,14)],longECK_by_regulonDB=long)






# Step IV: manually resove the remaining ones (still 8 unresolved ones)

is.na(temp5$longECK_by_regulonDB) %>% sum #19. Some of them are duplicates
unmapped19=temp5[is.na(temp5$longECK_by_regulonDB),]
unique(unmapped19$ECK) %>% length #14
##After this, still there are 14 unmapped genes 


#Manually map the above 
unmapped=unmapped19
unmapped$longECK_by_regulonDB=as.character(unmapped$longECK_by_regulonDB)

#peaD ECK0532 => (by searching regulonDB) ybcD . EcoliWiki doesn't have peaD as a synonym => from gene.txt -> ECK120026446
unmapped$longECK_by_regulonDB[unmapped$associated_gene_names=="peaD"]="ECK120026446"

#ycgH => (by searching regulonDB) ycgH_1, ycgH_2. They have the same ECK: ECK1157 => from gene.txt -> ycgH_1: ECK120003204 ycgH_2: ECK120003205
##Good news -> the 2 genes are not regulated (ref: generegulation_tmp.txt)

#gapC => (by searching regulonDB) gapC_1, gapC_2. They have the same ECK: ECK1409 but different b numbers (gapC_1: b1417 gapC_2: b1416) 
## => from gene.txt -> gapC_1: ECK120004456 gapC_2: ECK120002008
##Good news -> the 2 genes are not regulated (ref: generegulation_tmp.txt)

#arpB => (by searching regulonDB) arpB_1, arpB_2. They have the same ECK: ECK1718 but different b numbers (arpB_1: b1720 arpB_2: b1721) 
## => from gene.txt -> arpB_1: ECK120003514 arpB_2: ECK120003515
##Good news -> the 2 genes are not regulated (ref: generegulation_tmp.txt)

#yeeL => (by searching regulonDB) yeeL_1, yeeL_2. They have the same ECK: ECK1975 but different b numbers (yeeL_1: b1980 yeeL_2: b1979) 
## => from gene.txt -> yeeL_1: ECK120003652 yeeL_2: ECK120003651
##Good news -> the 2 genes are not regulated (ref: generegulation_tmp.txt)

#wbbL => (by searching regulonDB) wbbL_1, wbbL_2. They have the same ECK: ECK2025 but different b numbers (wbbL_1: b2031 wbbLL_2:b4540 ) 
## => from gene.txt -> wbbL_1:  ECK120001911 wbbL_2: ECK120026451
##Good news -> the 2 genes are not regulated (ref: generegulation_tmp.txt)

#yghX=> (by searching regulonDB) yghX_1, yghX_2. They have different ECKs: yghX_1: 	ECK2993,  yghX_2: ECK2994. In Nichols' it uses ECK2993
##=> from gene.txt -> yghX_1: ECK120004116
unmapped$longECK_by_regulonDB[unmapped$associated_gene_names=="yghX"]="ECK120004116"

#yghY => (by searching regulonDB) yghY=yghX_2=ECK2994 => from gene.txt ->  ECK120004117
unmapped$longECK_by_regulonDB[unmapped$associated_gene_names=="yghY"]="ECK120004117"

#phnE => (by searching regulonDB) phnE_1, phnE_2. They have different ECKs: phnE_1: 	ECK4097,  phnE_2: ECK4096. In Nichols' it uses ECK4096
##=> from gene.txt -> phnE_2: ECK120004320
unmapped$longECK_by_regulonDB[unmapped$associated_gene_names=="phnE"]="ECK120004320"

#ylbI and ECK4441 returns nothing from regulonDB and ecoliWiki
##=> For now I just assume that it's not regulated

#istR (In Nichols' original names for the data file it is: istR-1) => (by searching regulonDB) istR-1 (a dash, not underscore), istR: ECK4425 istR-1: no ECK (Phantom Gene). In Nichols' their ECK is: ECK5002
##=> For now I just assume that it's not regulated

#rybD and ECK5003 returns nothing from regulonDB and ecoliWiki
##=> For now I just assume that it's not regulated

#ryeF => (by searching regulonDB) micL => from gene.txt -> ECK125165376
unmapped$longECK_by_regulonDB[unmapped$associated_gene_names=="ryeF"]="ECK125165376"

#tpkE (In Nichols' original names for the data file it is: tpke70) => from gene.txt -> ECK120002617
unmapped$longECK_by_regulonDB[unmapped$associated_gene_names=="tpkE"]="ECK120002617"


##Note: If I want to figure out the unresolved 8 long ECK (which might be unnecessary), I have to look into the original Keio paper to see which gene they actually deleted


#Merge all the tables above into one:
tab1=corrected_merged_ECK_1st_table[!is.na(corrected_merged_ECK_1st_table$longECK_by_regulonDB),]
tab2=temp5[!is.na(temp5$longECK_by_regulonDB),]
tab3=unmapped
##dim(tab1)[1]+dim(tab2)[1]+dim(tab3)[1] #3986 #Check to see if I am using the right tables to merge

ECK_1st_table_withLongECK=rbind(tab1,tab2,tab3)


#Step V: Get the regulon annotations
generegulation_tmp=read.table("Data/regulon/generegulation_tmp.txt",sep="\t",quote="",stringsAsFactors = F) #no. of rows match the file: 4986-37+1=4950
class(generegulation_tmp) #data.frame

##For each long ECK of Nichols' strains get the regulator gene ID 
ECK_1st_table_withLongECK$longECK_by_regulonDB=ECK_1st_table_withLongECK$longECK_by_regulonDB %>% as.character

regulonAnnotation=sapply(ECK_1st_table_withLongECK$longECK_by_regulonDB,FUN=function(longECK){
  generegulation_tmp$V1[which(longECK %in% generegulation_tmp$V7)]
})


###Find a way to reduce 3986 back to 3979 (3986-3979=7 -> 7 ECKs have more than 1 JW names (according to the file: inline-supplementary-material-4(Baba et al., 2006).xls))
reduced=ECK_1st_table[,-which(names(ECK_1st_table) %in% c("JW","JW2","position","strain_availability"))]
reduced[duplicated(reduced),]

unique(ECK_1st_table$ids) %>% length
unique(ECK_1st_table$ECK) %>% length

###From the following object I know that evry id in Nichols' only maps to 1 long ECK
tmp=ECK_1st_table_withLongECK[ duplicated(ECK_1st_table_withLongECK$ids) | duplicated(ECK_1st_table_withLongECK$ids,fromLast=T),]

###Get unique id (total no.=3979) - long ECK
id_longECK=ECK_1st_table_withLongECK[,which(names(ECK_1st_table_withLongECK) %in% c("ids","longECK_by_regulonDB"))]
unique_id_longECK=id_longECK[!duplicated(id_longECK$ids),]
unique_id_longECK$longECK_by_regulonDB=as.character(unique_id_longECK$longECK_by_regulonDB)

idInNichols_regulatorID=sapply(unique_id_longECK$longECK_by_regulonDB,FUN=function(longECK){
  generegulation_tmp$V1[generegulation_tmp$V7==longECK]
})
names(idInNichols_regulatorID)=unique_id_longECK$ids

idInNichols_regulatorID=lapply(idInNichols_regulatorID,FUN=function(regulator){ #obj type conversion so it can be used for dist_TF_cumsum()
  if(identical(regulator,character(0))) regulator=as.character(NA)
  
  return(regulator)
})


used_id.regulatorID=data.frame()
for(i in 1:length(idInNichols_regulatorID)){
  id=names(idInNichols_regulatorID)[i]
  
  for(j in 1:length(idInNichols_regulatorID[[i]])){
    if(is.na(idInNichols_regulatorID[[i]][j])) break
    used_id.regulatorID=rbind(used_id.regulatorID,c(id,idInNichols_regulatorID[[i]][j]))
  }
}
colnames(used_id.regulatorID)=c("ids","regulator")
str(used_id.regulatorID)

save(used_id.regulatorID,file="Data/used_id.regulatorID.RData")





#Step VI: Get the transcription unit(operon) annotations
tu_objects_tmp=read.table("Data/regulon/tu_objects_tmp.txt",sep="\t",quote="",stringsAsFactors = F) #no. of rows match the file: 14013-46+1=13968
class(tu_objects_tmp) #data.frame

idInNichols_operonID=sapply(unique_id_longECK$longECK_by_regulonDB,FUN=function(longECK){
  tu_objects_tmp$V1[tu_objects_tmp$V7==longECK]
})
names(idInNichols_operonID)=unique_id_longECK$ids

idInNichols_operonID=lapply(idInNichols_operonID,FUN=function(operon){ #obj type conversion so it can be used for dist_TF_cumsum()
  if(identical(operon,character(0))) operon=as.character(NA)
  
  return(operon)
})




used_id.operonID=data.frame()
for(i in 1:length(idInNichols_operonID)){
  id=names(idInNichols_operonID)[i]
  
  for(j in 1:length(idInNichols_operonID[[i]])){
    if(is.na(idInNichols_operonID[[i]][j])) break
    used_id.operonID=rbind(used_id.operonID,c(id,idInNichols_operonID[[i]][j]))
  }
}
colnames(used_id.operonID)=c("ids","operon")

##I have decided to use operon ids from EcoCyc instead of from regulonDB (ECK12 numbers are not searchable from regulonDB)
#save(used_id.operonID,file="Data/used_id.operonID.RData")
