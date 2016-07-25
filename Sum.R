#to creat summary files (SUM) ----------------------------------------

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

#########################################################################
#### to creat summary files AA (AA SUM) #################################

library(seqinr)
library(dplyr)

setwd("/data/usrhome/LabCCKing/ccking01/Desktop/Ko_DENV/Combined/aligned")
subls <-list.files(getwd())

  for(i in 1:length(subls)){
  
  t1=Sys.time()  
  ah <- read.csv(subls[i])

  account.m =  matrix(0,26, ((dim(ah)[2]-1)/3))
  AA.sort   =  c("G","A","V","L","I","W","Y","D","H"
                 ,"N","E","K","Q","M","R","S","C","P"
                 ,'-',"*","X","F","FALSE","T"," TRUE","TRUE")
  rownames(account.m) = AA.sort
  
  
  for( k in 1:((dim(ah)[2]-1)/3)){
  
  jk=c() # for each aa position k
  h = (k-1)*3 + 1 + 1
  hh = h + 2
  ahk<-ah[, h:hh]
  
  
    ahkk = ahk[complete.cases(ahk),]
    
    if (nrow(ahkk) != 0){
      
      tahkk = as.vector(t(ahkk))
      tahkk = recode(tahkk, "TRUE" = "T", "N" = "*")
      
      jk = translate(tahkk)
    
      account.m[   match( data.frame( table( jk ) )$jk, rownames(account.m) ), k ] =
        
        data.frame(table( jk ))$Freq    
      
    } 

      }

  tf.bind.m         <- matrix(0,23,ncol(account.m))
  tf.bind.m [1:21,] <- account.m[1:21,]
  tf.bind.m [22,]   <- colSums(account.m[22:23,])
  tf.bind.m [23,]   <- colSums(account.m[24:26,])
  
  tfAA.sort =c("G","A","V","L","I","W","Y","D","H"
               ,"N","E","K","Q","M","R","S","C","P",'-',"*","X","F","T")
  rownames(tf.bind.m) = tfAA.sort

  write.csv(tf.bind.m, paste0("/data/usrhome/LabCCKing/ccking01/Desktop/Ko_DENV/Combined/AASUM/AASUM",subls[i],".csv" ))
  
  print(Sys.time()-t1)
  print(subls[i])

}



#Post Summarization QC  ---------------------------------------------

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
    

        
    
#Assay the deepness of the summarized data   --------------------------
    
setwd("/data/usrhome/LabCCKing/ccking01/Desktop/Ko_DENV/allDseq/SUM")
subls <-list.files(getwd())

mat<-matrix(ncol=3000)

for(i in 1:length(subls)){
 
  matt=matrix(ncol=3000)
  ah <- read.csv(subls[i])
  
  for(k in 3:length(ah)){
    
    tck<-sum(ah[,k][1], ah[,k][2], ah[,k][3], ah[,k][4])
    matt[,k-2]=tck
  }
  
  mat=rbind(mat, matt) 
  
}

mat=mat[-1,]
write.csv(mat, "Deep.csv")

#generate the figure ----------------------

ah<-read.csv(file.choose()) #ie, Deep.csv from last section
ahh <- as.matrix(ahh) 
ahh=ah[,-1] #to remove 1st column


#limit region to 880-2525 to reduce compupational space
ahhd <- ahh[,-1:-879]
ahhd <- ahhd[, -1647:-2121] 

library(reshape2)
ahhd.m <- melt(ahhd)

#library(easyGgplot2)
#ggplot2.stripchart(data=ahhd.m, xName='Var2', yName='value', shape=1) + 
#                   theme(axis.text.x=element_text(angle=90)) 
#the result is too complicated


#to generate the mean and minimum in each position
ahh.med<-apply(ahh, 2, median)
ahh.min<-apply(ahh, 2, min)
#ahh.max<-apply(ahh, 2, max)

ahl <- rbind(ahh.med, ahh.min)
ahl <- ahl[,-1:-879]
ahl <- ahl[,-1647:-2121]
ahl.m <- melt(ahl) #melt to fit the ggplot

bk <- attributes(ahhd)$names

xa=seq(0,17) 
xa=xa*100+21
b=seq(9,25)*100 #show labels: 900, 1000, 1100, ... 2500

#to callout the real "names" in original x label, e.g. bk[21] = v900
bkk=c() 
for(i in 1:17){
  bkk[i] <- bk[xa[i]]
} 

#generate figure

Fab<-ggplot(data=ahl.m, aes(Var2, value, group=Var1, color=Var1)) + 
            geom_line(size=2) + scale_x_discrete(breaks=bkk, labels=b) + 
            
            xlab("") +
            ylab("Numbers of total reads per site") + 
  
            theme_bw() +
            theme(axis.text.x=element_text(size=10), 
                  axis.title=element_text(size=16,face="bold")) +
            theme(legend.position="none")


