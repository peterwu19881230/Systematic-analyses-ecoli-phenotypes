conn=file("/Users/peterwu/TAMU Google Drive/OMP/Nichols-Data-mining/ECK_name_mapping/pwy-genes.dat")
linn <-readLines(conn) ##Shows incomplete final line but from my inspection the last line is ok
close(conn)

#Concatenate strings into 1 single string (that way I can parse more easily)
pathways.str=c()
for(i in 1:length(linn)){
  pathways.str=paste(pathways.str,linn[i],sep="")
}




