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
      tahkk = as.character(tahkk)
      
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


########## Post Mapping Quality Control ##################


##### To summarize No. of reads and reads with Indel #####
#


library(pasillaBamSubset)
library(stringr)
library(Rsamtools)


setwd("/data/usrhome/LabCCKing/ccking01/Ko/0810CLCmapping/0810bam/")

ls0 = list.files(getwd())

ls129 = ls0[1:129]
ls129c = strsplit(ls129, split="_con", fixed = T)
ls129clean = sapply(ls129c, function(x) head(x,1))
ls92 = ls0[130:length(ls0)]
ls92c = strsplit(ls92, split= " ", fixed = T)
ls92clean = sapply(ls92c, function(x) head(x,1))

bamls = c(ls129clean, ls92clean)                # modified name

rmcigar = c("N", "P", "=", "X")           # discard 4
rmcigarID = c("I", "D")                   # ID

mat=data.frame()

for (i in 1: length( list.files(getwd())  ) ){
  
  
  r.bam <- scanBam(BamFile( ls0[i] ))
  
  read.no =  length( r.bam[[1]]$seq )
  
  lcigar = strsplit( r.bam[[1]]$cigar, split ="" )
  
  rml1 = which( unlist(lapply(lcigar, function(x){any(rmcigar %in% x)})) == TRUE ) 
  rmlNPX = length(rml1)
  
  rml2 = which( unlist(lapply(lcigar, function(x){any(rmcigarID %in% x)})) == TRUE ) 
  rmlID = length(rml2)
  
  rml3 = which( r.bam[[1]]$pos > 3000 ) 
  rmlPos = length(rml3)
  
  rmlT = length( unique(c(rml1, rml2, rml3)) )
  
  mappedP = ( rmlT/read.no )*100
  
  
  mat.i = data.frame(bamls[i], ls0[i], read.no, rmlNPX, rmlID, rmlPos, rmlT, mappedP)  
  mat = rbind(mat, mat.i)
  
  print(bamls[i])
  
}


write.csv(mat, "/data/usrhome/LabCCKing/ccking01/Desktop/Ko_DENV/Reads.csv")


    
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


