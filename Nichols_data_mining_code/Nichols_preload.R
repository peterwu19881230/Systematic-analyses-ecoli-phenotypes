#Goal:Load commonly used stuff for the project of Nichols'

#Set stringsAsFactors = FALSE so things won't be automatically converted to factor
options(stringsAsFactors = FALSE) 




#Load commonly used packages. If not installed, install.
if(!require(pacman)){
  install.packages("pacman")
  library(pacman)
}
#Tidyverse: https://www.tidyverse.org/packages/
#Note: BiocManager has to be installed before its downstream packages
pacman::p_load(tidyverse,pbapply,BiocManager,readxl,factoextra,pheatmap,ComplexHeatmap,XML,reshape2,googlesheets,googledrive,
               circlize,ape,dendextend,sqldf,RMySQL,AnnotationHub,infotheo,stringi,purrr,ggfortify) #Note ggplot2 is loaded when tidyverse or factoextra is loaded 

#Uncomment the following if scripts using the following library need to be rerun

#if(!require(GOSemSim)){
#  BiocManager::install("GOSemSim")
#}

#if(!require(org.EcK12.eg.db)){
#  BiocManager::install("org.EcK12.eg.db")
#}




##Source all the self-defined functions
functions=dir("functions/")
paths=paste("functions/",functions,sep="")
for(func in paths){
  source(func)
}

#Source some of my utility functions from public GitHub
source("https://github.com/peterwu19881230/R_Utility/raw/master/general_r/functions/EscapeMeta.R")
source("https://github.com/peterwu19881230/R_Utility/raw/master/general_r/functions/tableSMY.R")




##Load useful data in "Data" directory
#Lazy loading ref: https://www.r-bloggers.com/lazy-load-with-archivist/
#Note: if the .RData is changed, delete the .rdb and ,rdx files so they can be regenerated (Otherwise, old info will be load from the .rdb and .rdx files)
Data=dir("Data/sourced")
paths=paste("Data/sourced/",Data[!( grepl(".rdx",Data,fixed = T) | grepl(".rdb",Data,fixed = T)) ],sep="")
for(path in paths){
  if(class(try(lazyLoad(path)))=="try-error"){
    lazyLoad = local({load(path); 
      environment()})
    tools:::makeLazyLoadDB(lazyLoad, path)
    lazyLoad(path)
  }else lazyLoad(path)
}










