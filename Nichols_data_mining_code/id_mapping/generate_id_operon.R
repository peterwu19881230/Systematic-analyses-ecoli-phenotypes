load("Data/transunits.dat.RData")


EcoCycGene_EcoCycOperon=data.frame()
for(i in seq(transunits.dat)){
  
  operon=names(transunits.dat[i])
  
  for(EcoCycID in transunits.dat[[i]]$`COMPONENTS - `){
    
    EcoCycGene_EcoCycOperon=rbind(EcoCycGene_EcoCycOperon,data.frame(EcoCycID=EcoCycID,operon=operon))
  }
}



save(EcoCycGene_EcoCycOperon,file="Data/EcoCycGene_EcoCycOperon.RData")