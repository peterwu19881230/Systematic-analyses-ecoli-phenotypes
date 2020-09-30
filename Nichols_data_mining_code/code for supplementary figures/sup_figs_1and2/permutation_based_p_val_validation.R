#validates that random pcc from generate_permutation_based_p_val.R is equivalent to: taking random n pairs -> calculate pcc -> take mean


setwd("~/Dropbox/OMP\ shared/peter_wu/1st_paper_nichols_reanalysis/github/Systematic-analyses-ecoli-phenotypes/Nichols_data_mining_code")
source("Nichols_preload.R")



#get annotation data (pwy, complex)
load("Data/cors_for_Pwy_Pcomplex.RData") 

##There are multiple data in it. Only these will be used: abs_cors_for_pwys, abs_cors_for_pcomplexes


if(!file.exists("Data/pcc_s_for_uniq_no_of_members.RData")){

#precompose the large no. of sampling pcc from complete data
#===========================================
start.time=Sys.time()

uniq_no_of_members=sort(unique(c(abs_cors_for_pwys$no_gene_used,abs_cors_for_pcomplexes$no_gene_used)))

large_N=10^3 #no. of random draws
pcc_s_for_uniq_no_of_members=list()
for(n in uniq_no_of_members){
  mean_abs_pcc_s=c()
  for(i in 1:1000){
    sub_dat=All_Data_NAimputed[sample(dim(All_Data_NAimputed)[1],n),]
    mean_abs_pcc=(sub_dat %>% t %>% cor %>% abs %>% as.dist %>% melt_dist)[[3]] %>% mean
    
    mean_abs_pcc_s=c(mean_abs_pcc_s,mean_abs_pcc)
    
  }
  
  pcc_s_for_uniq_no_of_members=c(pcc_s_for_uniq_no_of_members,
                                 list(mean_abs_pcc_s))
}

names(pcc_s_for_uniq_no_of_members)=uniq_no_of_members

end.time=Sys.time()
end.time-start.time #Time difference of 12.36659 mins

save(pcc_s_for_uniq_no_of_members,file="Data/pcc_s_for_uniq_no_of_members.RData")
#===========================================

}else{
  load("Data/pcc_s_for_uniq_no_of_members.RData")  
}



str(pcc_s_for_uniq_no_of_members)

#get all pcc
similarity_column="pcc"
pcc=strain1strain2_allAnnotations_allSimilarities[[similarity_column]]
random_abs_pcc=mean(abs(pcc)) #H0
random_abs_pcc



