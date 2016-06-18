
setwd("/data/usrhome/LabCCKing/ccking01/Desktop/Ko_DENV/aligned")
subls <-list.files(getwd())

for(i in 1:length(subls)){
  
  t1=Sys.time()  
  ah <- read.csv(subls[i])
  
  account.m <- matrix( 0,6,ncol( ah ) )
  nu.sort  = c( 'A','C','G','N',"TRUE",'T' )
  rownames( account.m ) = nu.sort
  colnames( account.m ) = colnames(ah)
  
  error = c()
  
  for( ll in 2:ncol(ah) ){
    if( length(  table(  ah[,ll] ) ) == 0 ){ error=1 } else {
      account.m[   match( data.frame( table(ah[,ll] ) )$Var1,rownames( account.m ) ),ll   ] = data.frame( table( ah[,ll] ) )$Freq
    } }
  
  final.m <- matrix( 0,5,ncol(ah) )
  rownames(final.m)=c('A','T','C','G','N')
  
  final.m[1,] = account.m[1,]
  final.m[2,] = apply(account.m[5:6,], 2, sum)
  final.m[3,] = account.m[2,]
  final.m[4,] = account.m[3,]
  final.m[5,] = account.m[4,]
  
  colnames( final.m ) = colnames( account.m )
  
  write.csv(final.m, paste0("/data/usrhome/LabCCKing/ccking01/Desktop/Ko_DENV/SUM",subls[i],".csv" ))
  print(Sys.time()-t1)
  print(subls[i])
}