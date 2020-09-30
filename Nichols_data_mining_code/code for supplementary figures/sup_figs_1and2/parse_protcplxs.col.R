#Goal: Parse pathways.dat and map the pathway names

currentDir=dirname(rstudioapi::getActiveDocumentContext()$path) #only works in RStudio console
setwd(currentDir)

ptcomplex_table=read.table(paste(currentDir,"protcplxs.col",sep="/"),comment.char = "#",fill=T,sep="\t",quote="",header=T,check.names = F)
pcomplexID_pcomplexName=ptcomplex_table[,1:2]


names(pcomplexID_pcomplexName)=c("pcomplexID","pcomplexName")

#clean italized words
for(i in 1:dim(pcomplexID_pcomplexName)[1]){
  
  while(grepl("</(.*)>",pcomplexID_pcomplexName[i,2])){
    
    for(to_be_removed in c("<i>","</i>","<sup>","</sup>","<I>","</I>","<SUP>","</SUP>","<sub>","</sub>","<small>","</small>")){
      pcomplexID_pcomplexName[i,2]=str_replace(pcomplexID_pcomplexName[i,2],to_be_removed,"")
    }
  }
  
  #cat(pcomplexID_pcomplexName[i,2]); cat("\n") #print each cleaned common name
}


#save(pcomplexID_pcomplexName,file="pcomplexID_pcomplexName.RData")
