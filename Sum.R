#to creat summary files (SUM)

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



####################### 

#QC

setwd("/data/usrhome/LabCCKing/ccking01/Desktop/Ko_DENV/allDseq/LOG")
subls <-list.files(getwd())

mat=data.frame()

  for(i in 1:length(subls)){
    
    
    ah <- read.csv(subls[i], stringsAsFactors = FALSE)
    
    qfiles<-ah$File[1:4]
    qNin<-sum(!is.na(ah$insert))
    qNdel<-sum(!is.na(ah$deletion))
    qNnull<-sum(!is.na(ah$zero))
    qNout<-sum(!is.na(ah$quality))
    qNTotal<-sum(qNin, qNdel, qNnull, qNout)

    qPin<-qNin/qNTotal*100
    qPdel<-qNdel/qNTotal*100
    qPnull<-qNnull/qNTotal*100
    qPout<-qNout/qNTotal*100

    qMeq<-mean(ah$quality)
    qMxq<-max(ah$quality)
    qMnq<-min(ah$quality)
    
    qDF <- data.frame(qfiles[1], qfiles[2], qfiles[3], qfiles[4], 
                      qNTotal, qPin, qPdel, qPnull, qPout, 
                      qMxq, qMnq, qMeq)
    mat=rbind(mat, qDF)
    
}
    colnames(mat) <- c("File1", "File2", "File3", "File4", "TotalReads", "PIn", "PDel", "PNull", "PReadIn", "MaxQ", "MinQ", "MeanQ")
    write.csv(mat, "QC.csv")
