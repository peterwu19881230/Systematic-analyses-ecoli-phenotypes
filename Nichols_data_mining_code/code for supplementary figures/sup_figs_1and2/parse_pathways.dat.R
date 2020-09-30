#Goal: Parse pathways.dat and map the pathway names

currentDir=dirname(rstudioapi::getActiveDocumentContext()$path) #only works in RStudio console
setwd(currentDir)

##Warnings:
##"KNOCKOUT-GROWTH-OBSERVATIONS,PRODUCT are not listed in the attribute list". Why? And there might be more.
##Attributes "HIDE-SLOT?" and "INTERRUPTED?" are not in any of the genes. What are they?



#first I have to change the format of the file to Unicode(UTF-8) (from BBEdit). Otherwise there will be warnings

library(stringr) ##load to use the str_extract()


#Define the attribute vector. I only care about some that I need:
att=c("TYPES - ","COMMON-NAME - ","SUB-PATHWAYS - ")



#Process the file
con = file("pathways.dat", "r")

pathways.dat=list()
i=1
temp.att.found=NULL
while ( TRUE ) {
  line = readLines(con, n = 1) ##readLines will put escape characters if the line being read contains any meta characters
  #print(line) ##print the line being processed
  
  
  ##stop parsing at the end of the file
  if ( length(line) == 0 ) {
    break
  }
  
  ##store if an ID is found
  ID=grepl("^UNIQUE-ID",line)
  if(ID==T){
    pathways.dat[[i]]=list()
    
    names(pathways.dat)[[i]]=sub("^UNIQUE-ID - ","",line)
    next
  }
  
  ##store the data if they are found. 
  att.found=str_extract(line,".+ - ")
  
  NoOfAtt=which(att %in% att.found)
  

  
  
  ###concatnate if the attributes is the same as the previous one(s)
  if(identical(att.found,temp.att.found)){
    pathways.dat[[i]][[temp.att.found]]=c(pathways.dat[[i]][[temp.att.found]],sub(att.found,"",line))
    next
  }
  
  ###store various attributes of the data
  if(identical(NoOfAtt,integer(0))==FALSE){
    pathways.dat[[i]][[att.found]]=sub(att.found,"",line)  
    #print(pathways.dat[[i]][[att.found]])
    
    temp.att.found=att.found ###temporalily store the att.found in case the next line is a continuation of the same attribute
    next  
  }
  
  
  
  ###append new lines if they are the continuation of the attribute
  #if(grepl("^/[^/]",line)==TRUE){
  #  pathways.dat[[i]][[temp.att.found]]=paste(pathways.dat[[i]][[temp.att.found]],line,sep="\n") 
  #  next
  #}
  
  
  ##Go to the next gene if // is found
  if(grepl("^//",line)==T){
    i=i+1
    next
  }
  
}
close(con)


#render part of the list into a table
pwyID_pwyName=data.frame()
for(i in seq(pathways.dat)){
  pwyID_pwyName=rbind(pwyID_pwyName,c(names(pathways.dat[i]),pathways.dat[[i]]$`COMMON-NAME - `))
}

names(pwyID_pwyName)=c("pwyID","pwyName")

#clean italized words
for(i in 1:dim(pwyID_pwyName)[1]){
  
  while(grepl("</(.*)>",pwyID_pwyName[i,2])){
    
    for(to_be_removed in c("<i>","</i>","<sup>","</sup>","<I>","</I>","<SUP>","</SUP>","<sub>","</sub>","<small>","</small>","<SUB>","</SUB>")){
      pwyID_pwyName[i,2]=str_replace(pwyID_pwyName[i,2],to_be_removed,"")
    }
  }
  
  #cat(pwyID_pwyName[i,2]); cat("\n") #print each cleaned common name
}




