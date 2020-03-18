#1. Render the original ECKs from Nichols' to a table: "ids" - "originalECKs" -"correctedECKs" -"associated_gene_names"
#2. Map other identifiers: JW, position, EcoCycID, b number, other synonyms, strain availability


#Sort the name column by ECK. Label duplicated strains with -1, -2

###Note: allData.csv is directly from from Nichols' supplemental
ECKs_rows_fixed=read.csv(file="Data/allData.csv",colClasses=c("character", rep("NULL", 324))) %>% t %>% as.vector ##colClasses is used to skip all the other columns
ECKs_rows_fixed=sort(ECKs_rows_fixed)

dup1=ECKs_rows_fixed[duplicated(ECKs_rows_fixed,fromLast=T)] %>% paste0("-1")
dup2=ECKs_rows_fixed[duplicated(ECKs_rows_fixed,fromLast=F)] %>% paste0("-2")

index_dup_1=which(duplicated(ECKs_rows_fixed,fromLast=T))
index_dup_2=which(duplicated(ECKs_rows_fixed,fromLast=F))

ECKs_rows_fixed[index_dup_1]=dup1
ECKs_rows_fixed[index_dup_2]=dup2




#Correct the strain names


#With SPA tag
SPAtag<-ECKs_rows_fixed[grep(pattern="- SPA",x=ECKs_rows_fixed)]
#With DAS tag
DAStag<-ECKs_rows_fixed[grep(pattern="- DAS",x=ECKs_rows_fixed)]
#9 point-mutant alleles of essential genes (plus corresponding ‘‘linked’’ strains for each of the 9 alleles, in which the antibiotic-resistance cassette was linked to the wild-type allele as a control)
mt_linked<-ECKs_rows_fixed[grep(pattern="Linked",x=ECKs_rows_fixed)]
#2 truncation of essential genes
truncations<-ECKs_rows_fixed[grep(pattern="Truncation",x=ECKs_rows_fixed)]


#others
isnot<-which(
  ECKs_rows_fixed %in%
    
    c(
      SPAtag,
      DAStag,
      mt_linked,
      truncations
    ) 
)

other_genes<-ECKs_rows_fixed[-isnot]
rm(isnot)

##In other_genes, I want to cut them into sub-groups

##Most common expression in the whole data
most_common_genes<-other_genes[
  c(
    grep(pattern="^ECK[0-9]{4}-[a-zA-Z]{3}[a-zA-Z0-9]$",x=other_genes)
  )
  ]




##Only ECKXXXX <= This is where we have to add gene names. I don't know why there weren't gene names.
onlyECK<-other_genes[
  grep(pattern="^ECK[0-9][0-9][0-9][0-9]$",x=other_genes)
  ]

##Add gene names for ECK_only manully: from OMP shared/Nichols_phenotypic-landscape/Nichols-Strains-Sorted.xlsx (search for all tabs)
##Some gene names might not be in Nichols-Strains-Sorted.xlsx?
##format: ECKXXXX->ECKXXXX-YYYY

onlyECK_gene_name_added<-c(
  "ECK0012-HTGA", "ECK0017-HOKC", "ECK0266-YKGN", "ECK0320-YAHH", "ECK0359-YAIF", "ECK0367-YKIB", "ECK0369-YAIU", "ECK0503-YBBV",
  "ECK0619-YBEM", "ECK0679-YBFH", "ECK1128-YMFH", "ECK1132-CROE", "ECK1159-YMGG", "ECK1160-YMGH", "ECK1453-YNCM", "ECK1933-YEDM",
  "ECK1990-YOEE", "ECK2132-YOHH", "ECK2331-YCFT", "ECK2636-YPJM", "ECK2637-YPJM", "ECK2647-YPJC", "ECK2650-YGAR", "ECK2651-YGAC",
  "ECK2652-YGAD", "ECK2675-YGAY", "ECK2854-YGEL", "ECK2856-YGEN", "ECK2859-YGEQ", "ECK2994-YGHY", "ECK3474-YHIK", "ECK3672-YSDC",
  "ECK3675-GLVC", "ECK3769-YIFN", "ECK3802-YZCX", "ECK4097-PHNE", "ECK4219-YZFA", "ECK4265-YJGW", "ECK4330-YJIQ", "ECK4334-YJIV",
  "ECK4426-TISA"
)

'
Note:
*No gene name is found for ECK0012 within Nichols-Strains-Sorted.xlsx. However I did find it inside OMP shared/Nichols_phenotypic-landscape/Sorted_ECK_by_Peter/Keio_strains_with_verified_JW. But I fotgot where I got the JW numbers and the matched genes from 
*ECK1132 has gene name synonyms: croE, ymfT (http://ecoliwiki.net/colipedia/index.php/croE:Gene)
*Ecoliwiki indicates that ECK2636 and ECK2637 are the same 
'

##remaining
isnot<-which(
  
  other_genes %in%
    c(
      most_common_genes,
      onlyECK
    ) 
)

remaining_genes<-other_genes[-isnot]

rm(isnot)

##check the total length (should be 3979):
##length(SPAtag)+length(DAStag)+length(mt_linked)+length(truncations)+length(most_common_genes)+length(onlyECK)+length(remaining_genes)


##Start making columns

##Sorted by category
sorted_ECK<-c(SPAtag,DAStag,mt_linked,truncations,most_common_genes,onlyECK,remaining_genes)
sorted_ECK_missing_gene_names_added<-c(SPAtag,DAStag,mt_linked,truncations,most_common_genes,onlyECK_gene_name_added,remaining_genes)




##ECK only
library(stringr)
ECK_only<-str_extract(sorted_ECK_missing_gene_names_added,"^ECK[0-9][0-9][0-9][0-9]")

##gene name only (This is based on the assumption that all the full names start with ECKXXXX-)
ECKXXXXremoved<-str_extract(sorted_ECK_missing_gene_names_added,"-(.*)") 


associated_gene_names<-str_extract(ECKXXXXremoved,"[A-Z0-9]{3,}") 
###These are just gene names, not genotypes (Genes that are associated with the genotypes. Some genotypes are actually more like wild-type )
###This is a little dangerous since whether the names are actual gene names cannot be verified. Will verify after mapping to gene names in the .GAF file

#manualy correct: 
##istR -> istR-1 (this is a phantom gene)
associated_gene_names[3977]="ISTR-1"

##ECK4418/2590-RYFD(WITHCLPB) involves 2 genes
associated_gene_names[3974]="RYFD,CLPB"

rm(ECKXXXXremoved)

##bind them 
ECK_name_columns<-cbind(sorted_ECK,ECK_only,associated_gene_names)


##Replace names for exceptions (The rules above don't apply)
###ECK0005-TP2
ECK_name_columns[ECK_name_columns[,"sorted_ECK"]=="ECK5005-TP2","associated_gene_names"]<-"TP2"

###ECK0086-A-Linked: This is just the wild-type E.coli that has the antibiotic resistance gene linked to MURE
ECK_name_columns[ECK_name_columns[,"sorted_ECK"]=="ECK0086-A - Linked","associated_gene_names"]<-"MURE"

###ECK0086-B-Allele: This is the actual mutant (the mutant is called MurE1)
ECK_name_columns[ECK_name_columns[,"sorted_ECK"]=="ECK0086-B - Allele","associated_gene_names"]<-"MURE1"

###ECK0086-C-SPA: This is also a kind of wt that has SPA+KAN after the orf of interest
ECK_name_columns[ECK_name_columns[,"sorted_ECK"]=="ECK0086-C - SPA","associated_gene_names"]<-"MURE"


###ECK0055-D-IMP4213 - Allele
ECK_name_columns[ECK_name_columns[,"sorted_ECK"]=="ECK0055-D-IMP4213 - Allele","associated_gene_names"]<-"imp"

###ECK4142-4143-ECNAB
ECK_name_columns[ECK_name_columns[,"sorted_ECK"]=="ECK4142-4143-ECNAB","associated_gene_names"]<-"ecnA,ecnB"


#Step2: Corrected strain names (Labeled duplicated strains with -1, -2. Added missing ECK names)


#Original strain names (created with an id column in which ids serve as identifiers)
ids=1:3979 ##identifiers
originalECKs=read.csv(file="Data/allData.csv",colClasses=c("character", rep("NULL", 324))) %>% t %>% as.vector ##colClasses is used to skip all the other columns
id_ECKs<-cbind(ids,originalECKs) #Bind 2 columns




##Get the indices from ECK_name_column and bind -- this sapply() take about 10 sec
indices<-sapply(id_ECKs[,2],FUN=function(x){
  
  ###Escape all metas in regex
  x<-gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", x)
  
  length<-length(grep(x,ECK_name_columns[,"sorted_ECK"]))
  if(length<2&identical(grep(x,ECK_name_columns[,"sorted_ECK"]), integer(0))==FALSE){  
    ###if there are >=2 (duplicated ECK names), I want to deal with it later
    ###integer(0) is what gets returned when grep() finds nothing
    return(grep(x,ECK_name_columns[,"sorted_ECK"]))   
  }else(return(paste("Have to fix->",x))) 
}
)
indices<-unlist(indices) #simplify the result so all the integers are all in 1 list instead of many
indices<-unname(indices) #Get rid of names 

###manaully clean the duplicated strains (Total 12. But ECK2019 - HISA was found twice because there are: ECK2019 - HISA, ECK2019-HISA - SPA)
###
indices[indices=="Have to fix-> ECK1323-YMJC'"]<-c(3739,3740) ###The first index is from -1, the second from -2. 
indices[indices=="Have to fix-> ECK1824-MGRB"]<-c(3789,3790)
indices[indices=="Have to fix-> ECK3357-YHFL"]<-c(3901,3902)
indices[indices=="Have to fix-> ECK3531-DPPA"]<-c(3916,3917)
indices[indices=="Have to fix-> ECK2613-SMPA"]<-c(3855,3856)
indices[indices=="Have to fix-> ECK2019-HISA"]<-1708
indices[indices=="Have to fix-> ECK0295-YKGO"]<-c(3660,3661)
indices[indices=="Have to fix-> ECK1544-GNSB"]<-c(3763,3764)
indices[indices=="Have to fix-> ECK4410-YDGU"]<-c(3968,3969)
indices[indices=="Have to fix-> ECK4415-YPFM"]<-c(3970,3971)
indices[indices=="Have to fix-> ECK1556-HOKD"]<-c(3766,3767)
indices[indices=="Have to fix-> ECK4416-RYFB"]<-c(3972,3973)
indices[indices=="Have to fix-> ECK2593-A-YFIO\\* - Truncation"]<-c(132,133)

###convert indices to numberic
indices<-as.numeric(indices)


###order + cbind() + remove unnecessary columns from ECK_name_columns + correct the gene names (All capital -> the first 3 small)
ordered_ECK_name_columns<-ECK_name_columns[indices,]
id_ECKs_ordered_ECK_name_columns<-cbind(id_ECKs, ordered_ECK_name_columns)
ids_originalECK_geneName<-id_ECKs_ordered_ECK_name_columns[,c("ids","originalECKs","associated_gene_names")]


ids_originalECK_geneName[,"associated_gene_names"]<-sapply(ids_originalECK_geneName[,"associated_gene_names"],function(x) ###(All capital -> the first 3 small)
{
  sub("[A-Za-z]{3}",str_extract(x,"[A-Za-z]{3}") %>% tolower,x)
}
)

ids_originalECK_geneName[,"associated_gene_names"][ids_originalECK_geneName[,"associated_gene_names"]=="TP2"]="tp2"  ###clean the exception
ids_originalECK_geneName=as.data.frame(ids_originalECK_geneName)



###(Updated 9/3/2019) manually correct the fake ECKs made by Nichols (search for gene names in EcoCyc and put back the right ECK)

###====================================================================================================================================================
###The list of 22 fake ECKs (I pulled them out by manually looking at the sorted ECK of the original data file that Nichols' provides) -> correct ECK
###ECK4466-MOKC ECK0018
###ECK4472-YOAI ECK1786
###ECK5000-SROH ECK4505
###ECK5001-SGRT ECK4477
###ECK5002-ISTR-1 G0-10202 This doesn't have ECK number. Use this as the EcoCyc ID
###ECK5003-RYBD ECK4621
###ECK5004-RYEF ECK4574
###ECK5005-TP2 G0-8894 This doesn't have ECK number. Use this as the EcoCyc ID
###ECK5006-TPKE70 G0-8906 This doesn't have ECK number. Use this as the EcoCyc ID
###ECK5007-YKGR ECK4486
###ECK5008-YMIB ECK4487
###ECK5009-YMJD ECK4488
###ECK5010-YNBG ECK4489
###ECK5011-YOAJ ECK4490
###ECK5012-YOAK ECK4491
###ECK5013-YOBI ECK4492
###ECK5014-YOEI ECK4493
###ECK5015-YOHP ECK4494
###ECK5016-YPDK ECK4495
###ECK5017-YQCG ECK4497
###ECK5018-YQEL ECK4498
###ECK5019-YQFG ECK4499


###correct the rows of ECK_1st_table for the above strains

##correct the 2nd, 3rd columns of "ids_originalECK_geneName" generated at the end of clean_names.R (associated_gene_names shouldn't have to change)
fake_ECK_genes=c("ECK4466-MOKC","ECK4472-YOAI","ECK5000-SROH","ECK5001-SGRT","ECK5002-ISTR-1","ECK5003-RYBD","ECK5004-RYEF","ECK5005-TP2","ECK5006-TPKE70",
                 "ECK5007-YKGR","ECK5008-YMIB","ECK5009-YMJD","ECK5010-YNBG","ECK5011-YOAJ","ECK5012-YOAK","ECK5013-YOBI","ECK5014-YOEI","ECK5015-YOHP",
                 "ECK5016-YPDK","ECK5017-YQCG","ECK5018-YQEL","ECK5019-YQFG")

corrected_ECK=c("ECK0018","ECK1786","ECK4505","ECK4477","","ECK4621","ECK4574","","","ECK4486","ECK4487","ECK4488","ECK4489","ECK4490","ECK4491","ECK4492","ECK4493","ECK4494","ECK4495",
                "ECK4497","ECK4498","ECK4499")

temp=ids_originalECK_geneName
temp$correctedECKs=temp$originalECKs
temp=temp[,c("ids","originalECKs","correctedECKs","associated_gene_names")]

for(i in 1:length(fake_ECK_genes)){
  index_=grep(fake_ECK_genes[[i]],ids_originalECK_geneName$originalECKs)
  replacement=str_replace(fake_ECK_genes[[i]],"^ECK[0-9]{4}",corrected_ECK[[i]])
  if(corrected_ECK[[i]]=="") replacement=NA #if no ECK,give NA
  temp[index_,"correctedECKs"]=replacement
}



###give NA for genes that have multiple deleted genes. This prevents any annotations being wrongly attached later
list_=c("ECK4418/2590-RYFD(WITHCLPB)","ECK1205/1207/1209-RDLABC","ECK1205/1207/1209/3525-RDLABCD","ECK2068/2069/2908/3041/4420-SIBABCDE","ECK4142-4143-ECNAB")
temp$correctedECKs[temp$correctedECKs %in% list_]=NA



##Check to see if the table is correct using Supplemental Information: Bacterial Strains from Nichols' paper
## 3737 Keio, 117 SPA, 5 DAS, 9 alleles, 9 linked controls, 2 truncations, and 100 sRNAs
#sum(grepl(" - SPA",temp$sorted_ECK_missing_gene_names_added)) 
#sum(grepl(" - DAS",temp$sorted_ECK_missing_gene_names_added)) 
#sum(grepl(" - Linked",temp$sorted_ECK_missing_gene_names_added))
#sum(grepl(" - Truncation",temp$sorted_ECK_missing_gene_names_added))



ids_originalECK_geneName=temp
###====================================================================================================================================================



#Start from here.......




#Add Unique IDs and all synonyms based on genes.dat (from EcoCyc)
##Extract info from  genes.dat object created by parse.genes.dat.R
##(!)Note that some of the ECKs are pseudogenes. They are not in genes.dat

load("Data/genes.dat.RData")

uniqueID=list()
synonyms=list()
accessions=list()
accessions_bNumber=list()
for(i in 1:length(genes.dat)){
  
  ###retrieve Unique ID
  uniqueID[i]=names(genes.dat)[[i]]
  
  ###retrieve b number
  accessions_bNumber[i]=paste(genes.dat[[i]][["ACCESSION-1 - "]],collapse = ",")
  
  ###retrieve ECK 
  accessions[i]=paste(genes.dat[[i]][["ACCESSION-2 - "]],collapse = ",") 
  ###subsetting a list element that doesn't exist won't return "subscript out of bounds". Instead, it returns NULL
  
  ###retrieve synonyms if any
  synonyms[i]=paste(genes.dat[[i]][["SYNONYMS - "]],collapse = ",") 
  
  
}


indices=list()
for(i in 1:3979){
  ECK=str_extract(ids_originalECK_geneName$correctedECKs[i],"^ECK[0-9]{4}") 
  
  index=which(accessions %in% ECK)
  
  if( index %>% length >= 1 &is.numeric(index)==T){
    indices[[i]]=index
  }else{
    indices[[i]]=""
  } 
  
}
rm(ECK,index) 

###This loop checks if there are more than 2 ECKs found
two.or.more=list()
j=1
for(i in 1:length(indices)){
  if(length(indices[[i]])>=2){
    two.or.more[j]=indices[i]
    names(two.or.more)[j]=paste("index from indices= ",i)
    j=j+1
  }
}
two.or.more ###indices[687]=1859 2015  (ECK1366 was found twice in genes.dat) 

###This block prepares data to cbind(). Mainly because ECK1366 found 2 data in genes.dat
EcoCycID=c()
other.synonyms=c()
bNumber=c()
for(i in 1:3979){
  EcoCycID[i]=paste(uniqueID[indices[[i]]],collapse=",") ###Allowing concatnating a vector into 1 string
  other.synonyms[i]=paste(synonyms[indices[[i]]],collapse=",") 
  if(grepl("^,",other.synonyms[i])==T | other.synonyms[i]=="NULL" | other.synonyms[i]==""){
    other.synonyms[i]=NA ###If there are no synomym, the result will be NA
  }  ###remove "," without any synonym attached
  bNumber[i]=paste(accessions_bNumber[indices[[i]]],collapse=",")
}



##The result obj. Note that ECK1366-LOMR' has 2 EcoCyc identifiers => we should use G6692 because LOMR' means the disruption on the N terminal (if it's 'LOMR it's C terminal)
ECK_1st_table=cbind(ids_originalECK_geneName,
                    EcoCycID=EcoCycID,
                    other.synonyms=other.synonyms,bNumber,stringsAsFactors=F)

##Correct the row: ECK1366-LOMR'
ECK_1st_table[ECK_1st_table$ids=="687",]$EcoCycID="G6692"
ECK_1st_table[ECK_1st_table$ids=="687",]$bNumber="b1369"

##Add Ecocyc ID where strains don't have ECK number: tp2: G0-8894, tpkE70: G0-8906 , istR-1: G0-10202
ECK_1st_table[ECK_1st_table$associated_gene_names=="tp2","EcoCycID"]="G0-8894"
ECK_1st_table[ECK_1st_table$associated_gene_names=="tpkE70","EcoCycID"]="G0-8906"
ECK_1st_table[ECK_1st_table$associated_gene_names=="istR-1","EcoCycID"]="G0-10202"


ECK_1st_table=as.data.frame(ECK_1st_table)




#The following code: 
##1. add a column for JW numbers by an original file provided by Nichols' paper
##2. add columns for the link to CGSC synonym, strain availability, position of deletion
##(problems are explained in each block)


library(stringr)
##1
strainInfo=read_xls("Data/inline-supplementary-material-4(Baba et al., 2006).xls",sheet=1,skip=3,col_names = F)[,1:3] #it says "...8 more problems" but the data look fine
colnames(strainInfo)=c("ECK","gene","JW")
ECK_1st_table$ECK=str_extract(ECK_1st_table$correctedECKs,"^ECK[0-9]{4}")

temp=merge(ECK_1st_table,strainInfo[,c("ECK","JW")],by="ECK",all.x=T,all.y=F)
temp$ids=as.character(temp$ids)
temp=temp[order(as.numeric(temp$ids)),]

### A problem:
sum(duplicated(temp$ids)) #7 ECKs have more than 1 JW names (according to the file: inline-supplementary-material-4(Baba et al., 2006).xls). 

##This table gives the info of all those ECKs that have 2 JWs => The genes deleted are the same but the positions are different
duplicated.ECK=temp[temp$ECK %in% c("ECK0614","ECK1007","ECK1157","ECK1409","ECK1718","ECK1975","ECK2025"),]
##For example ECK0614: 
## http://cgsc2.biology.yale.edu/Strain.php?ID=115890 
## http://cgsc2.biology.yale.edu/Strain.php?ID=115891

##2
##Add 3 columns (CGSC synonym, strain availability, position of deletion) from: JWstrainInfoFromCGSC.txt created by parse_CGSC_Strains.R
JWstrainInfoFromCGSC=read.table("Data/JWstrainInfoFromCGSC.txt")
colnames(JWstrainInfoFromCGSC)=c("JW","JW2","position","strain_availibility")
temp2=merge(temp,JWstrainInfoFromCGSC,by="JW",all.x=T,all.y=F)


temp2=temp2[,c("ids","ECK","JW","JW2","associated_gene_names","position","originalECKs","EcoCycID","bNumber","other.synonyms","strain_availibility")]

### A problem:
##JW1878 has 2 forms from CGSC: JW1878-4, JW1878-4/pWT2594 => I should remove the one with plasmid. It's unlikely Nichols used the strain with a plasmid encoding gfp and other things
JWstrainInfoFromCGSC$JW[which(duplicated(JWstrainInfoFromCGSC$JW))]

temp2=temp2[-which(temp2$JW2=="JW1878-4/pWTZ594"),]
temp2$bNumber[temp2$bNumber=="NULL"]=NA #clean the bNumber column

ECK_1st_table=temp2


## I have to correct things that are not in the 3737 Keio (They do have ECK at the front but they weren't from Keio. However, the code above considered them to be):
hobbs=read_xls("Data/hobbs2009_tables1_smallrnasandproteins.xls",sheet=1,skip=9,col_names=F)[,1:2]
colnames(hobbs)=c("EcoCycID","genotype")
hobbsInNichols=merge(ECK_1st_table,hobbs,by="EcoCycID")

##Can consider adding strain background column here


save(ECK_1st_table,file="Data/ECK_1st_table.RData")








##(!)Some of the ECKs don't match an EcoCyc gene ID. Here they are.

###The table
#NULLgenes=ECK_1st_table[ECK_1st_table$EcoCycID=="NULL",]

###A list of ECKs for those genes
#Nullgenes.list=sapply(NULLgenes$`Original Name`,str_extract,pattern="^ECK[0-9][0-9][0-9][0-9]") %>% as.character

###Output the list to a text file
###write.table(Nullgenes.list,"Nullgenes.list.txt",quote=F,row.names = F,col.names = F)

###Create a smart table and map them to EcoCyc gene ID (Object ID) and pathways. Luckily, none of them are in a pathway
#eightyNulls=read.table("Data/80NULLs.txt")

##Also, I Got all the ECKs and made a smart table on EcoCyc 


###A list of all ECKs: ECKs_rows_fixed (defined in clean_names.R)

#ECKs.list=sapply(ECKs_rows_fixed,str_extract,pattern="^ECK[0-9][0-9][0-9][0-9]") 
###write.table(ECKs.list,"Data/ECKs.list.txt",quote=F,row.names = F,col.names = F)
