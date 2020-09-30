#Ran this on my 15' mac

#install.packages("infotheo")
library(dplyr)


melt_dist=function(distance){
  m1=as.matrix(distance)
  m1[upper.tri(m1)]=NA #I suppose that no distance can be NA, so I can use this to do filtering
  diag(m1)=NA
  library(reshape2)
  m2=melt(m1) #I suppose that no distance can be NA, so I can use this to do filtering
  if(class(m2$Var1)!="character") m2$Var1=as.character(m2$Var1) #This is to prevent numeric names (the class should still be "character") being converted to "numeric" by melt()
  if(class(m2$Var2)!="character") m2$Var2=as.character(m2$Var2) #This is to prevent numeric names (the class should still be "character") being converted to "numeric" by melt()
  
  m2=m2[!is.na(m2[,3]) | is.nan(m2[,3]),] # "| is.nan(m2[,3])" is used because I want to keep the NaN values
  
  m2=m2[,c(2,1,3)] ##reorder the columns
  names(m2)=c("object_1","object_2","value") ##name the columns
  return(m2)
}


load("Ternary_Data_324cutff_NAremoved.RData")
str(Ternary_Data_324cutff_NAremoved)
Ternary_Data_324cutff_NAremoved_5_FDR=Ternary_Data_324cutff_NAremoved

load("ternary_data_various_FDR.RData")
##str(Ternary_Data_324cutff_NAremoved_1_FDR) #I don't necessary have to use this
str(Ternary_Data_324cutff_NAremoved_10_FDR)
str(Ternary_Data_324cutff_NAremoved_15_FDR)
str(Ternary_Data_324cutff_NAremoved_20_FDR)


#for 3 bins:
##just use Ternary_Data_324cutff_NAremoved_5_FDR
ternary_3_bin=Ternary_Data_324cutff_NAremoved_5_FDR

#for 5 bins
ternary_5_bin=Ternary_Data_324cutff_NAremoved_5_FDR+
              Ternary_Data_324cutff_NAremoved_10_FDR

#for 7 bins
ternary_7_bin=Ternary_Data_324cutff_NAremoved_5_FDR+
              Ternary_Data_324cutff_NAremoved_10_FDR+
              Ternary_Data_324cutff_NAremoved_15_FDR


#for 9 bins
ternary_9_bin=Ternary_Data_324cutff_NAremoved_5_FDR+
              Ternary_Data_324cutff_NAremoved_10_FDR+
              Ternary_Data_324cutff_NAremoved_15_FDR+
              Ternary_Data_324cutff_NAremoved_20_FDR



#calculate pairwise MI
##Mutual Information (Sho)
##Ref: https://cran.r-project.org/web/packages/infotheo/infotheo.pdf




#=======The following takes >1 hour=======


library(infotheo)

nbin=c(5,7,9)
i=1
for(ternary_dat in list(ternary_5_bin,ternary_7_bin,ternary_9_bin)){
  start.time = Sys.time()
  mi_ternary=melt_dist(natstobits(mutinformation(as.data.frame((t(ternary_dat))))))
  
  end.time = Sys.time()
  end.time-start.time #Time difference of 1.059961 hours
  
  save(mi_ternary,file=paste0("mi_ternary_",nbin[i],"_bin",".RData"))
  i=i+1
  
  cat(paste0("i=",i," is done"))
}






