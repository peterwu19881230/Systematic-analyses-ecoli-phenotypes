#the code is accoring to the reviewer:

# One potential way to implement a formal test would be, for each grouping of interest (pathways for Fig S1, complexes for Fig S2) 
# to perform 1000 random draws of expression profiles for the same number of genes as are actually in that pathway/complex, and calculate
# the fraction of those draws that have a higher |PCC| than the actual value for that pathway -- this would provide a permutation-based p value 
# (which could then be subjected to FDR correction). This would allow flagging of significant groupings (pathways/complexes) just to see how 
# many there are. 


setwd("~/Dropbox/OMP\ shared/peter_wu/1st_paper_nichols_reanalysis/github/Systematic-analyses-ecoli-phenotypes/Nichols_data_mining_code")
source("Nichols_preload.R")


#get all |pcc|
similarity_column="pcc"
pcc=abs(strain1strain2_allAnnotations_allSimilarities[[similarity_column]])

#get annotation data (pwy, complex)
load("Data/cors_for_Pwy_Pcomplex.RData") 

##There are multiple data in it. Only these will be used: abs_cors_for_pwys, abs_cors_for_pcomplexes



#precompose the large no. of sampling pcc from complete data
#===========================================
start.time=Sys.time()

uniq_no_of_members=sort(unique(c(abs_cors_for_pwys$no_gene_used,abs_cors_for_pcomplexes$no_gene_used)))

large_N=10^3 #no. of random draws

pcc_s_for_uniq_no_of_members=list()
for(n in uniq_no_of_members){
  mean_pcc_s=c()
  for(i in 1:1000){
    mean_pcc=sample(pcc,n) %>% mean
    mean_pcc_s=c(mean_pcc_s,
                 mean_pcc)
  }
  
  pcc_s_for_uniq_no_of_members=c(pcc_s_for_uniq_no_of_members,
                                 list(mean_pcc_s))
}


names(pcc_s_for_uniq_no_of_members)=uniq_no_of_members


end.time=Sys.time()
end.time-start.time #Time difference of 2.520909 mins
#===========================================


random_pcc=mean(pcc) #H0


fdr_based_p_val_list=list()
for(annot in list(abs_cors_for_pwys,abs_cors_for_pcomplexes)){
  
  pcc_groups_of_interest=annot
  
  
  p_val=c()
  for(i in 1:dim(pcc_groups_of_interest)[1]){
    
    pcc_group_of_interest=mean(pcc_groups_of_interest$avg_pcc[i])  #Ha
    no_of_members=pcc_groups_of_interest$no_gene_used[i]
    
    
    pcc_s=pcc_s_for_uniq_no_of_members[[which(names(pcc_s_for_uniq_no_of_members)==as.character(no_of_members))]]
    
    diff=abs(pcc_s)-pcc_group_of_interest #mu_ha - mu_h0 (vecorized operation)
    
    
    p_val=c(p_val,
            mean(diff>=0))
    
    paste0(i," out of ", dim(pcc_groups_of_interest)[1]," is complete","\n") %>% cat
  }
  
  
  
  fdr_based_p_val=p.adjust(p_val,method="fdr") 
  fdr_based_p_val_list=c(fdr_based_p_val_list,list(fdr_based_p_val))
}

str(fdr_based_p_val_list)


save(fdr_based_p_val_list,file="Data/fdr_based_p_val_list.RData")


