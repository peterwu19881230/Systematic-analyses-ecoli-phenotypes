#Compute different types of Nichols' phenotype data


#Primary data preparation
All_Data=read.csv("Data/allData.csv")[,-1]
Ternary_Data_NAnotimputed=TernaryConvert(matrix=All_Data,thresh=3.463)




## Impute all NA with mean of all scores
avg=mean(as.numeric(as.matrix(All_Data)),na.rm=T)
All_Data_NAimputed=All_Data
indexForNA=which(is.na(All_Data_NAimputed),arr.ind=T)
for(i in 1:dim(indexForNA)[1]){
  row=indexForNA[i,1]
  col=indexForNA[i,2]
  All_Data_NAimputed[row,col]=avg
}
which(is.na(All_Data_NAimputed),arr.ind=T)##Verify that there is no NA anymore
save(All_Data_NAimputed,file="Data/sourced/All_Data_NAimputed.RData")


##Convert All_Data_NAimputed to Ternary
 

#Ternary data using 324 cutoffs for each condition
start.time=Sys.time()

FDR_cond=list()

total=3979 
for(i in 1:324){
  
  n=1
  dat=All_Data_NAimputed[,i] #Some conditions have NAs, but they constitute only a small amount. I think ignoring them would be fine.
  fdr=numeric(total)
  for(q in dat){ #Iterate through all fitness scores can find me all possible fractions
    
    if(q>=0) leng=sum(dat>=q) 
    if(q<0)  leng=sum(dat<=q)
    
    fdr[n]=(1-pnorm(abs(q)))/(leng/total)
    n=n+1
  }
  
  FDR_cond[[i]]=cbind(score=dat,FDR=fdr) %>% as.data.frame
}

end.time=Sys.time()
end.time-start.time #Time difference of 22.55891 secs
save(FDR_cond,file="Data/sourced/FDR_cond.RData") 


min_forEachCond=sapply(FDR_cond,FUN=function(q_fdr){
  q_fdr$score[q_fdr$FDR<=0.05] %>% abs %>% min   
})


dat=All_Data
Ternary_Data_324cutff=matrix(numeric(3979*324),nrow=3979)
for(i in 1:324){
  Ternary_Data_324cutff[,i]=TernaryConvert(dat[,i] %>% as.matrix,thresh =min_forEachCond[i] )
}


rownames(Ternary_Data_324cutff)=1:3979



#Impute the NA values with 0. Ref: https://gist.github.com/Jfortin1/d4888d68359a36fbda60
impute_matrix=function(matrix){
  missing=which(is.na(matrix), arr.ind=T)
  if (length(missing)!=0){
    for (j in 1:nrow(missing)){
      matrix[missing[j,1],missing[j,2]]=0
    }
  }
  matrix
}

Ternary_Data_324cutff_NAremoved=impute_matrix(Ternary_Data_324cutff)
rownames(Ternary_Data_324cutff_NAremoved)=1:3979
save(Ternary_Data_324cutff_NAremoved,file="Data/sourced/Ternary_Data_324cutff_NAremoved.RData")


##Add various ternary data based on differnt FDR: 1%, 10%, 15%, 20%
##========================================================================

convertToTernary_byFDR=function(dat=All_Data,FDR_cond=FDR_cond,FDR=0.05){ #FDR_cond is defined in phenotypeCutoff.R
  
  #determine the cutoff for each condition
  min_forEachCond=sapply(FDR_cond,FUN=function(q_fdr){
    q_fdr$score[q_fdr$FDR<=FDR] %>% abs %>% min   
  })
  
  #generate the ternary data by using the cutoffs to filter
  Ternary_Data_324cutff=matrix(numeric(3979*324),nrow=3979)
  for(i in 1:324){
    Ternary_Data_324cutff[,i]=TernaryConvert(dat[,i] %>% as.matrix,thresh =min_forEachCond[i] )
  }
  
  rownames(Ternary_Data_324cutff)=1:3979
  
  #impute NAs with 0. Ref: https://gist.github.com/Jfortin1/d4888d68359a36fbda60
  matrix=Ternary_Data_324cutff
  
  missing=which(is.na(matrix), arr.ind=T)
  if (length(missing)!=0){
    for (j in 1:nrow(missing)){
      matrix[missing[j,1],missing[j,2]]=0
    }
  }
  
  Ternary_Data_324cutff_NAremoved=matrix
  
  return(Ternary_Data_324cutff_NAremoved)
}

Ternary_Data_324cutff_NAremoved_1_FDR=convertToTernary_byFDR(dat=All_Data,FDR_cond=FDR_cond,FDR=0.01)
Ternary_Data_324cutff_NAremoved_10_FDR=convertToTernary_byFDR(dat=All_Data,FDR_cond=FDR_cond,FDR=0.1)
Ternary_Data_324cutff_NAremoved_15_FDR=convertToTernary_byFDR(dat=All_Data,FDR_cond=FDR_cond,FDR=0.15)
Ternary_Data_324cutff_NAremoved_20_FDR=convertToTernary_byFDR(dat=All_Data,FDR_cond=FDR_cond,FDR=0.20)

save(Ternary_Data_324cutff_NAremoved_1_FDR,
     Ternary_Data_324cutff_NAremoved_10_FDR,
     Ternary_Data_324cutff_NAremoved_15_FDR,
     Ternary_Data_324cutff_NAremoved_20_FDR,
     file="Data/ternary_data_various_FDR.RData")

##========================================================================




#Collapse 324 Conditions in Nichols' by the 114 treatments, and create necessary objs for further CorrVSAnnot analysis

#Criteria: if there is either 1 or -1, use that. If both are found or there are only 0s + NAs , switch to 0

dat=Ternary_Data_NAnotimputed

#All_Data obj contains the conditions in the same order as in allData.csv (allData.csv is directly derived from Nichols' supplementary data that contain scores)
condName=names(All_Data)

#Manual collapse of conditions by indices
##Note: the names of those indices are exactly the same as those in this Nichols' supplementary: TableS1_ListOfConditions.xls
uniqueChemIndex=list(
  'SDS+EDTA'=1:3,'Cold shock'=4:6,'Heat shock'=7:10,
  A22=11:14,'Actinomycin D'=15:18,'Acetate (M9)'=19,
  Acriflavine=20:21,Amikacin=22:24,Amoxicillin=25:28,
  Ampicillin=29:32,Anaerobic=33,Azidothymidine=34:36,
  Azithromycin=37:39,Benzalkonium=40:42,Bicyclomycin=43:44,
  'Bile salts'=45:48,Bleomycin=49:52,'CCCP (Carbonyl cyanide 3-chlorophenylhydrazone)'=53:55,
  'CHIR-090'=56:60,Cefoxitin=61:64,'Calcofluor (F3543    
Fluorescent Brightener 28) '=65, #I intended to leave this guy like this (yes, an Enter after F3543) to match the name of the original excel file
  'Cecropin B'=66:67,Cefsulodin=68:71,Cerulenin=72:75,Chlorpromazine=76:79,Cholate=80:83,
  Cisplatin=84:86,Clarythromycin=87:90,'Cobalt stress-CoCl2'=91:92,
  'Copper stress-CuCl2'=93:95,Deoxycholate=96:98,Dibucaine=99:101,
  Doxorubicin=102:103,EDTA=104:106,'Epigallocatechin gallate (EGCG)'=107:109,
  Epinephrine=110:112,Erythromycin=113:116,'Ethidium Bromide'=117:119,
  Fosfomycin=120,'Fusidic acid'=121:124,'N-acetyl Glucosamine'=125,
  'Fosfomycin +Glucose 6P'=126:127,'Glucosamine (M9)'=128,'Glucose (M9)'=129,
  'Glycerol (M9)'=130,Hydroxyurea=131:133,Indolicidin=134,Isoniazid=135:137,
  Levofloxacin=138,MMS=139,'Maltose (M9)'=140, Mecillinam=141:144,
  Methotrexate=145:146,Minocycline=147:149,'Mitomycin C'=150,'NH4Cl (MOPS)'=151,
  'Nickel stress-NiCl2'=152:153,Nigericin=154:156,Nitrofurnatoin=157:161, #Their typo: Nitrofurnatoin -> Nitrofurantoin
  Norepinephrine=162:163,'Propidium iodide'=164:166,'Phenazine methosulfate (PMS)'=167:169,
  'Paraquat dichloride'=170:174,'Hydrogen peroxide'=175:178,Phleomycin=179:181,Procaine=182:185,
  Puromycin=186:188,Pyocyanin=189:191,Radicicol=192:194,
  SDS=195:199,Streptomycin=200,Streptonigrin=201:203,
  'Succinate (M9)'=204,Sulfamonomethoxine=205:206,Taurocholate=207:209,
  Theophylline=210:211,Thiolactomycin=212:214,Tobramycin=215:218,
  'Triclosan/Irgasan'=219,'Triton X-100'=220:222,Tunicamycin=223:225,
  Vancomycin=226:228,Verapamil=229:231,Aztreonam=232:233,
  Bacitracin=234:236,Carbenicillin=237:239,'Cefsulodin + Mecillinam'=240,
  Cefaclor=241:243,Ceftazidime=244:246,Chloramphenicol=247:250,
  Ciprofloxacin=251:253,'Cycloserine D'=254,Doxycycline =255:258,
  EGTA=259:262,EtOH=263:265,Gentamicin=266:267,
  'Iron excess-FeSO4'=268,'Iron starvation-FeSO4'=269,NaCl=270:273,
  'Nalidixic acid'=274:277,Norfloxacin=278:280,Novobiocin=281:286,
  Oxacillin=287:289,'basic pH (TAPS)'=c(290,295:297),'acidic pH (MES-HOMOPIPES)'=291:294,
  'Polymyxin B'=298:301,Rifampicin=302:303,Spectinomycin=304:305,
  Spiramycin=306:308,Sulfamethizole=309:311,Tetracycline=312:315,
  Trimethoprim=316:319,'Trimethoprim + Sulfamethizole'=320,UV=321:324
)

save(uniqueChemIndex,file="Data/sourced/uniqueChemIndex.RData")





##Collapse 324 Conditions in Nichols' by the 114 treatments


##Criteria: if there is either 1 or -1, use that. If both are found or there are only 0s + NAs , switch to 0

dat=Ternary_Data_324cutff

Ternary_Data_324cutff_condCollapsed=matrix(,nrow=3979,ncol=length(uniqueChemIndex))
rownames(Ternary_Data_324cutff_condCollapsed)=1:3979

for(i in 1:length(uniqueChemIndex)){
  
  conds=dat[,uniqueChemIndex[[i]]] %>% as.data.frame #as.data.frame is to prevent treatment with only 1 condition from crashing
  
  treatment=apply(conds,1,FUN = function(row){
    if( (1 %in% row) && !(-1 %in% row) ) return(1) 
    if( (-1 %in% row) && !(1 %in% row) ) return(-1) 
    return(0) # If neither of the above critera were met, use 0 
    
  })
  
  Ternary_Data_324cutff_condCollapsed[,i]=treatment
  
}


save(Ternary_Data_324cutff_condCollapsed,file="Data/sourced/Ternary_Data_324cutff_condCollapsed.RData")
