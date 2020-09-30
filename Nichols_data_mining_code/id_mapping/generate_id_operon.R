#Goal: Parse parse.transunits.dat

#first I have to change the format of the file to Unicode(UTF-8) (from BBEdit). Otherwise there will be warnings

library(stringr) ##load to use the str_extract()

#Define the attribute vector. 
att=c("UNIQUE-ID - ","TYPES - ","COMMON-NAME - ","CITATIONS - ","COMMENT - ","COMPONENT-OF - ","COMPONENTS - ","CREDITS - ","DATA-SOURCE - ","DBLINKS - ","DOCUMENTATION - ",
      "EXTENT-UNKNOWN? - ","HIDE-SLOT? - ","INSTANCE-NAME-TEMPLATE - ","LEFT-END-POSITION - ","MEMBER-SORT-FN - ","PATHOLOGIC-NAME-MATCHER-EVIDENCE - ",
      "PATHOLOGIC-PWY-EVIDENCE - ","REGULATED-BY - ","RIGHT-END-POSITION - ","SYNONYMS - ","UNMAPPED-COMPONENT-OF - ")

#Process the file
con = file("id_mapping/transunits.dat", "r")

transunits.dat=list()
i=1
temp_att_found=NULL
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
    transunits.dat[[i]]=list()
    
    names(transunits.dat)[[i]]=sub("^UNIQUE-ID - ","",line)
    next
  }
  
  ##store the data if they are found. 
  att.found=str_extract(line,".+ - ")
  
  
  NoOfAtt=which(att %in% att.found)
  
  
  ###concatnate if the attributes is the same as the previous one(s)
  if(identical(att.found,temp_att_found)){
    transunits.dat[[i]][[temp_att_found]]=c(transunits.dat[[i]][[temp_att_found]],sub(att.found,"",line))
    next
  }
  
  ###store various attributes of the data
  if(identical(NoOfAtt,integer(0))==FALSE){
    transunits.dat[[i]][[att.found]]=sub(att.found,"",line)  
    temp_att_found=att.found ###temporalily store the att.found in case the next line is a continuation of the same attribute
    next  
  }
  
  ###append new lines if they are the continuation of the attribute
  if(grepl("^/[^/]",line)==TRUE){
    transunits.dat[[i]][[temp_att_found]]=paste(transunits.dat[[i]][[temp_att_found]],line,sep="\n") 
    next
  }
  
  
  ##Go to the next gene if // is found
  if(grepl("^//",line)==T){
    i=i+1
    next
  }
  
}
close(con)

save(transunits.dat,file="Data/transunits.dat.RData")







load("Data/transunits.dat.RData")


EcoCycGene_EcoCycOperon=data.frame()
for(i in seq(transunits.dat)){
  
  operon=names(transunits.dat[i])
  
  for(EcoCycID in transunits.dat[[i]]$`COMPONENTS - `){
    
    EcoCycGene_EcoCycOperon=rbind(EcoCycGene_EcoCycOperon,data.frame(EcoCycID=EcoCycID,operon=operon))
  }
}



save(EcoCycGene_EcoCycOperon,file="Data/EcoCycGene_EcoCycOperon.RData")