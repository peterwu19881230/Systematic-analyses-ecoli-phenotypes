#scrape KEGG modules + map to Nichols' ids
#Refs:
#https://www.analyticsvidhya.com/blog/2017/03/beginners-guide-on-web-scraping-in-r-using-rvest-with-hands-on-knowledge/
#http://bradleyboehmke.github.io/2015/12/scraping-html-text.html

library(rvest)
library(stringr)
library(dplyr)

## Get url for all Nichols' strains (genes)
baseURL="https://www.kegg.jp/dbget-bin/www_bget?eco:"
load("Data/ECK_1st_table.RData")
ids_bNumber=ECK_1st_table[,c("ids","bNumber")] %>% distinct #This is for mapping to Nichols' ids
  
## Start scraping for modules of each gene
start.time=Sys.time()



library(foreach)
library(doParallel)
library(parallel)
library(doSNOW)

numCores=detectCores()-1
cl=makeCluster(numCores)
registerDoSNOW(cl)
clusterExport(cl,c("baseURL","ids_bNumber"))

modules=foreach(bNumber=ids_bNumber$bNumber) %dopar%{
  url=paste(baseURL,bNumber,sep="")
  webpage=xml2::read_html(url)
  table_text=rvest::html_text(rvest::html_nodes(webpage,"table")) #Modules can be obtained from tables
  info=stringr::str_extract(table_text,"^eco_M[0-9]{1,6}(.*)") # (.*) captures for anything that follows
  modules=info[!is.na(info)]
  
  modules
}

stopCluster(cl)

end.time=Sys.time()
end.time-start.time #Time difference of 4.449409 mins

#(11/9/2019)why is this shit keep happening? I tried both on my mac and PC. This script used to work.
#Error in { : 
#task 1113 failed - "Timeout was reached: [www.kegg.jp] Operation timed out after 10003 milliseconds with 0 out of 0 bytes received"
#(11/10/2019) I tried one more time this day and it worked but it seems modules from KEGG have drastically changed


names(modules)=ids_bNumber$ids #This is mapping to Nichols' ids
KEGGmodulesForNichols=modules
save(KEGGmodulesForNichols,file="Data/KEGGmodulesForNichols.RData")


#How many modules are pathways? (Some of them are protein complexes)
modules=unlist(KEGGmodulesForNichols) %>% unique #total No. = 196 (new no. from 11/10/2019: 90)



















