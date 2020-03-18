#Goal: parse the excel file from copying-pasting from CGSC: http://cgsc2.biology.yale.edu/StrainQuery.php
##(The following link won't lead to the correct page. Only by searching from the query page will)


# I have to do setwd() to change to the current working directory first


library(readxl)
library(stringr)
#I can copy&paste from the search result of CGSC directly, and then make an excel file
currentDir=dirname(rstudioapi::getActiveDocumentContext()$path)
strainInfo=read_excel(paste(currentDir,"keioStrainsFromCGSC.xlsx",sep="/"),sheet=1)
strainInfo$allele=sapply(strainInfo$Genotype,FUN=function(genotype){ #This guy gets like: "xerC757(del)::kan". 757 is the allele No.
  geneWithKan=str_extract(genotype,"Δ[a-zA-Z]{2,5}[A-Z]{0,1}-{0,1}[0-9]{1,4}::kan")  
  #This should only search and replace the 1st pattern found (because I put ::kan at the end)
  geneWithKan=sub("Δ","",geneWithKan) #This removes "Δ"
  geneWithKan=sub("::kan","(del)::kan",geneWithKan) #This adds (del) infront of "::kan"
  })

strainInfo$synonym=strainInfo$Name
strainInfo$Name=str_extract(strainInfo$Name,"JW[0-9]{4}")


strainInfo=strainInfo[,c("Name","synonym","allele")]


#The following gets the hyperlink burried inside the excel file
##Ref: https://stackoverflow.com/questions/24149821/extract-hyperlink-from-excel-file-in-r/24152082
library(XML)

## rename file to .zip (I deleted this .zip file after running this script)
my.excel.file="keioStrainsFromCGSC.xlsx"
my.zip.file <- sub("xlsx", "zip", my.excel.file)
file.copy(from = paste(currentDir,my.excel.file,sep="/"), to = paste(currentDir,my.zip.file,sep="/"))

## unzip the file (this creates a bunch of files that is not needed later, so I deleted it after finishing this script)
unzip(paste(currentDir,my.zip.file,sep="/")) 

## unzipping produces a bunch of files which we can read using the XML package
## assume sheet1 has our data
xml <- xmlParse("xl/worksheets/sheet1.xml")

## finally grab the hyperlinks
hyperlinks <- xpathApply(xml, "//x:hyperlink/@display", namespaces="x")

links=sapply(hyperlinks,FUN=function(link){link[[1]][1]})

strainInfo=cbind(strainInfo,url=links)


#Quality check:

##Check that there is no redundancy for synonyms
sum(strainInfo$synonym==strainInfo$Name)

##How many NA are there for each column
sum(is.na(strainInfo$allele)) #11
NAalleleData=strainInfo[is.na(strainInfo$allele),] #These guys don't have allele

sum(is.na(strainInfo$url)) #0



##The only weird synonym I found (The others are just JWXXXX-X):
##JW1878-4/pWTZ594

write.table(strainInfo,paste(currentDir,"JWstrainInfoFromCGSC.txt",sep="/"),row.names = F,col.names = F,quote=F,sep="\t")




