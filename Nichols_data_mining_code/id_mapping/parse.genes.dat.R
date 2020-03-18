#Goal: Parse genes.dat. Resulting data structure: genes.dat$"gene accession"$"different types of data"

##Warnings:
##"KNOCKOUT-GROWTH-OBSERVATIONS,PRODUCT are not listed in the attribute list". Why? And there might be more.
##Attributes "HIDE-SLOT?" and "INTERRUPTED?" are not in any of the genes. What are they?



#first I have to change the format of the file to Unicode(UTF-8) (from BBEdit). Otherwise there will be warnings

library(stringr) ##load to use the str_extract()


#Define the attribute vector. There are 23.
att=c("TYPES - ","COMMON-NAME - ","ACCESSION-1 - ","ACCESSION-2 - ","CENTISOME-POSITION - ","CITATIONS - ","COMMENT - ","COMPONENT-OF - ",
      "COPY-NUMBER - ","CREDITS - ","DATA-SOURCE - ","DBLINKS - ","DOCUMENTATION - ","LAST-UPDATE - ","LEFT-END-POSITION - ",
      "MEMBER-SORT-FN - ","PATHOLOGIC-NAME-MATCHER-EVIDENCE - ", "PATHOLOGIC-PWY-EVIDENCE - ","REGULATED-BY - ","RIGHT-END-POSITION - ",
      "SYNC-W-ORTHOLOG - ","SYNONYMS - ","UNMAPPED-COMPONENT-OF - ","PRODUCT - ")

#Process the file
con = file("id_mapping/genes.dat", "r")

genes.dat=list()
i=1
temp.att.found=NULL
while ( TRUE ) {
    line = readLines(con, n = 1) ##readLines will put escape characters if the line being read contains any meta characters
    print(line) ##print the line being processed
    
    
    ##stop parsing at the end of the file
    if ( length(line) == 0 ) {
      break
    }
    
    ##store if an ID is found
    ID=grepl("^UNIQUE-ID",line)
    if(ID==T){
      genes.dat[[i]]=list()
                          
      names(genes.dat)[[i]]=sub("^UNIQUE-ID - ","",line)
      next
      }
    
    ##store the data if they are found. 
    att.found=str_extract(line,".+ - ")
    
    NoOfAtt=which(att %in% att.found)
    
    
    ###concatnate if the attributes is the same as the previous one(s)
    if(identical(att.found,temp.att.found)){
      genes.dat[[i]][[temp.att.found]]=c(genes.dat[[i]][[temp.att.found]],sub(att.found,"",line))
      next
      }
    
    ###store various attributes of the data
    if(identical(NoOfAtt,integer(0))==FALSE){
      genes.dat[[i]][[att.found]]=sub(att.found,"",line)  
      temp.att.found=att.found ###temporalily store the att.found in case the next line is a continuation of the same attribute
      next  
    }
    
    
    
    ###append new lines if they are the continuation of the attribute
    if(grepl("^/[^/]",line)==TRUE){
      genes.dat[[i]][[temp.att.found]]=paste(genes.dat[[i]][[temp.att.found]],line,sep="\n") 
      next
    }
    
    
    ##Go to the next gene if // is found
    if(grepl("^//",line)==T){
      i=i+1
      next
      }
  
}
close(con)

save(genes.dat,file="Data/genes.dat.RData")





