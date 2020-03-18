#use the mapped EcoCyc id to get UniProt ID

#Map the UniProt IDs to Nichols' strains


##Get product from genes.dat obj. 
###Some of them doesn't have product (makes sense because not all genes produce proteins)
###Some of them have more than 2 products
load("Data/ECK_1st_table.RData")
load("Data/genes.dat.RData")
EG.product.table=data.frame()
for(i in 1:dim(ECK_1st_table)[1]){ #here i is not Nichols' ID
  
  EcoCycID=ECK_1st_table$EcoCycID[i]
  
  product=paste("ECOCYC:",genes.dat[[EcoCycID]]$`PRODUCT - `,sep="")
  
  if(length(product)>=2) print(product)
  
  
  if(is.null(product)==T){  ##If there is no product
    #print(i); print("=> Doesn't have product") 
    next
  }
  
  
  EG.product.table=rbind(EG.product.table,
                         cbind(ECK_1st_table$ids[i],EcoCycID,product), 
                         stringsAsFactors=F
  )
  
}

#clean the propagated rows caused by propagated rows in ECK_1st_table (ECK_1st_table has 3986 rows, not 3979)
EG.product.table=unique(EG.product.table)
names(EG.product.table)[1]="ids"

##This shows ECK with more than 2 products
###EG.product.table$ids[duplicated(EG.product.table$ids)]
###ECK.id=EG.product.table$ids[duplicated(EG.product.table$ids)] %>% unique %>% as.numeric
###EG.MoreThan2=EG.product.table[EG.product.table$ids %in% ECK.id,]




#Retrieve UniProt ID based on an R package
#=========================================================================================================
##ref: https://bioconductor.org/packages/release/bioc/vignettes/UniProt.ws/inst/doc/UniProt.ws.pdf
##ref: https://bioconductor.org/packages/release/bioc/manuals/UniProt.ws/man/UniProt.ws.pdf

##BiocManager::install("UniProt.ws")
library(UniProt.ws)

up=UniProt.ws(taxId=83333) #83333 is the taxon id for E. coli K-12


#I have to divide EG.product.table first so every entry is <100 rows. Otherwise the result cannot be retrieved
keys=EG.product.table$product
divide_by=50
key_list=list()
for(i in 1:(length(keys)%/%divide_by)){
  key_list[[i]]=keys[1:divide_by+divide_by*(i-1)]
}
key_list=c(key_list,list(tail(keys,length(keys)%%divide_by))) #append the last chunk
  

biocyc_uniprot=data.frame()
for(key_fragment in key_list){
  df=select(up, keys=key_fragment, columns=c("UNIPROTKB"),
                        keytype="BIOCYC")
  biocyc_uniprot=rbind(biocyc_uniprot,df)
}
#=========================================================================================================


id.UniProt=data.frame(EG.product.table$ids,Data=biocyc_uniprot$UNIPROTKB)


save(id.UniProt,file="Data/id.UniProt.RData")
