setwd("/data/usrhome/LabCCKing/ccking01/Ko/Combined_E")

# WARNING HERE - each sp has 23 sample = 23*2*2 files
sp = 1

assign(paste0("fd", sp), list.files(getwd())[((sp-1)*92+1):(sp*92)])

fda = paste0("fd", sp) 
fd = get(paste0("fd", sp))

######################

library("ShortRead")
library("Biostrings")
library("seqinr")

######################

Ref <- readDNAStringSet("/data/usrhome/LabCCKing/ccking01/Desktop/Ko_DENV/01/2001Con.fas") 
RefSeq <- Ref[[1]][1:3000]

######################

np = 5000 #number of splits

######################

NoFd <- length(fd)  #number of files in the folder
Rfd = seq(0, (NoFd/4-1))*4+1

######################

for(Rfdi in 1: length(Rfd)) { #for each sample; each Rfdi indicates one sample
  
  
  t1=Sys.time()
  
  ## 4 QC vectors + primers
  
  outlierseqI <- c() #insersion
  outlierseqD <- c() #deletion
  zero <- c() #low quality
  zzero <- c() #valid reads
  
  Pf<-c(884, 901, 1310, 1327, 1616, 1635, 2026, 2044)
  Pr<-c(1357, 1376, 1752, 1770, 2035, 2053, 2504, 2524)
  
  ############  F1 
  
  ## call FsqF 1 
  
  FsqF <- readFastq(fd[Rfd[Rfdi]]) #r1
  
  ## create dummy matrix
  Dum = rep(1:np, length=length(FsqF))
  for(i in 1:np){
    assign(paste("mxF",i,sep=""), matrix( NA, 1, length(RefSeq)))} 
  
  ## FsqF 1 mapping and create mxFF1
  
  for (i in 1:length(FsqF)){
    
    sseq <- sread ((FsqF)[ i ])
    
    BestMatch <- pairwiseAlignment(RefSeq, sseq, type = ("local-global"),
                                   subjectQuality<- PhredQuality(FsqF[i]),
                                   substitutionMatrix = "BLOSUM62", gapOpening = 10, gapExtension = 5)
    
    if ( width( toString( BestMatch@pattern[1])) > width(BestMatch@pattern[1])   ) {   outlierseqI[ length( outlierseqI )+1 ] = i  } else 
      if ( width( toString( BestMatch@subject[1])) > width(BestMatch@subject[1])   ) {   outlierseqD[ length( outlierseqD )+1 ] = i  } else 
        if ( BestMatch@score < 0 ){  zero[ length( zero )+1 ] = as.numeric ( i )  }  else
          
        {
          zzero [ length(zzero)+1 ] <- BestMatch@score
          
          rtt <- c(  rep(  NA,  BestMatch@pattern[1]@range@start-1 ) ,
                     unlist( strsplit( toString( BestMatch@subject[1] ), split='') ) ) 
          
          assign(paste("mxF",Dum[i],sep=""), 
                 rbind(get(paste("mxF",Dum[i],sep="")), 
                       t( matrix( c(rtt, rep( NA, ncol(mxF1) - length(rtt) )) ) )))
          
        }
    
    
    
  }
  
  mxFlist<-mget(paste0("mxF",1:np), env=globalenv())
  mxFF1<-do.call("rbind", mxFlist)
  
  
  ############  F2
  
  ## call FsqF 2 
  
  FsqF <- readFastq(fd[Rfd[Rfdi] + 1 ]) #R1
  
  ## create dummy matrix
  Dum = rep(1:np, length=length(FsqF))
  for(i in 1:np){
    assign(paste("mxF",i,sep=""), matrix( NA, 1, length(RefSeq)))}
  
  ## FsqF 2 mapping and create mxFF2
  
  for (i in 1:length(FsqF)){
    
    sseq <- sread ((FsqF)[ i ])
    
    BestMatch <- pairwiseAlignment(RefSeq, sseq, type = ("local-global"),
                                   subjectQuality<- PhredQuality(FsqF[i]),
                                   substitutionMatrix = "BLOSUM62", gapOpening = 10, gapExtension = 5)
    
    if ( width( toString( BestMatch@pattern[1])) > width(BestMatch@pattern[1])   ) {   outlierseqI[ length( outlierseqI )+1 ] = i  } else 
      if ( width( toString( BestMatch@subject[1])) > width(BestMatch@subject[1])   ) {   outlierseqD[ length( outlierseqD )+1 ] = i  } else 
        if ( BestMatch@score < 0 ){  zero[ length( zero )+1 ] = as.numeric ( i )  }  else
          
        {
          zzero [ length(zzero)+1 ] <- BestMatch@score
          
          rtt <- c(  rep(  NA,  BestMatch@pattern[1]@range@start-1 ) ,
                     unlist( strsplit( toString( BestMatch@subject[1] ), split='') ) ) 
          
          assign(paste("mxF",Dum[i],sep=""), 
                 rbind(get(paste("mxF",Dum[i],sep="")), 
                       t( matrix( c(rtt, rep( NA, ncol(mxF1) - length(rtt) )) ) )))
          
        }
    
    
    
  }
  
  mxFlist<-mget(paste0("mxF",1:np), env=globalenv())
  mxFF2<-do.call("rbind", mxFlist)
  
  mxFFF <- rbind(mxFF1, mxFF2)
  
  
  for(k in 1:4){
    
    kk = (k-1)*2+1
    mxFFF[ ,Pf[kk]:Pf[kk+1]]=NA
    
  }
  
  print(Sys.time()-t1)
  t2=Sys.time()
  
  ######################
  
  ## call FsqR 1 
  
  FsqR <- readFastq(fd[Rfd[Rfdi] + 2 ]) #r2
  
  ## create dummy matrix
  
  Dum = rep(1:np, length=length(FsqR))
  for(i in 1:np){
    assign(paste("mxR",i,sep=""), matrix( NA, 1, length(RefSeq)))}
  
  ## FsqR 1 mapping and create mxRR1
  
  for (i in 1:length(FsqR)){
    
    
    sseq <- sread ((FsqR)[ i ])
    
    BestMatch <- pairwiseAlignment(RefSeq, reverseComplement(sseq), type = ("local-global"),
                                   subjectQuality<- PhredQuality(FsqR[i]),
                                   substitutionMatrix = "BLOSUM62", gapOpening = 10, gapExtension = 5)
    
    if ( width( toString( BestMatch@pattern[1])) > width(BestMatch@pattern[1])   ) {   outlierseqI[ length( outlierseqI )+1 ] = i  } else 
      if ( width( toString( BestMatch@subject[1])) > width(BestMatch@subject[1])   ) {   outlierseqD[ length( outlierseqD )+1 ] = i  } else 
        if ( BestMatch@score < 0 ){  zero[ length( zero )+1 ] = as.numeric ( i )  }  else
          
        {
          zzero [ length(zzero)+1 ] <- BestMatch@score
          
          rtt <- c(  rep(  NA,  BestMatch@pattern[1]@range@start-1 ) ,
                     unlist( strsplit( toString( BestMatch@subject[1] ), split='') ) ) 
          
          assign(paste("mxR",Dum[i],sep=""), 
                 rbind(get(paste("mxR",Dum[i],sep="")), 
                       t( matrix( c(rtt, rep( NA, ncol(mxR1) - length(rtt) )) ) )))
          
        }
    
    
  }
  
  mxRlist<-mget(paste0("mxR",1:np), env=globalenv())
  mxRR1<-do.call("rbind", mxRlist)
  
  ############  R2
  
  ## call FsqR 2 
  
  FsqR <- readFastq(fd[Rfd[Rfdi] + 3 ]) #R2
  
  ## create dummy matrix
  
  Dum = rep(1:np, length=length(FsqR))
  for(i in 1:np){
    assign(paste("mxR",i,sep=""), matrix( NA, 1, length(RefSeq)))}
  
  ## FsqR 2 mapping and create mxRR2
  
  for (i in 1:length(FsqR)){
    
    
    sseq <- sread ((FsqR)[ i ])
    
    BestMatch <- pairwiseAlignment(RefSeq, reverseComplement(sseq), type = ("local-global"),
                                   subjectQuality<- PhredQuality(FsqR[i]),
                                   substitutionMatrix = "BLOSUM62", gapOpening = 10, gapExtension = 5)
    
    if ( width( toString( BestMatch@pattern[1])) > width(BestMatch@pattern[1])   ) {   outlierseqI[ length( outlierseqI )+1 ] = i  } else 
      if ( width( toString( BestMatch@subject[1])) > width(BestMatch@subject[1])   ) {   outlierseqD[ length( outlierseqD )+1 ] = i  } else 
        if ( BestMatch@score < 0 ){  zero[ length( zero )+1 ] = as.numeric ( i )  }  else
          
        {
          zzero [ length(zzero)+1 ] <- BestMatch@score
          
          rtt <- c(  rep(  NA,  BestMatch@pattern[1]@range@start-1 ) ,
                     unlist( strsplit( toString( BestMatch@subject[1] ), split='') ) ) 
          
          assign(paste("mxR",Dum[i],sep=""), rbind(get(paste("mxR",Dum[i],sep="")), 
                                                   t( matrix( c(rtt, rep( NA, ncol(mxR1) - length(rtt) )) ) )))
          
        }
    
    
  }
  
  mxRlist<-mget(paste0("mxR",1:np), env=globalenv())
  mxRR2<-do.call("rbind", mxRlist)
  
  mxRRR <- rbind(mxRR1, mxRR2)
  
  for(k in 1:4){
    
    kk = (k-1)*2+1
    mxRRR[ ,Pr[kk]:Pr[kk+1]]=NA
    
  }
  
  
  mxx=rbind(mxFFF, mxRRR)
  
  write.csv(mxx, paste0("/data/usrhome/LabCCKing/ccking01/Desktop/Ko_DENV/aligned/",fd[Rfd[Rfdi]],".csv"))
  
  
  print(Sys.time()-t2)
  ######################
  
  id<-c(fd[Rfd[Rfdi]], fd[Rfd[Rfdi] + 1 ], fd[Rfd[Rfdi] + 2 ], fd[Rfd[Rfdi] + 3 ])
  td <- Sys.time()-t1
  xx <- list(id, td, outlierseqI, outlierseqD, zero, zzero)
  nxx <- sapply(xx, length)
  seq.xx <- seq_len(max(nxx))
  mat <-sapply(xx, "[", i=seq.xx)
  colnames(mat) <-c("Files","Time","insert","deletion","zero","quality")
  
  write.csv(mat, paste("/data/usrhome/LabCCKing/ccking01/Desktop/Ko_DENV/LOG/","LOG",fd[Rfd[Rfdi]],".csv",sep=""))
  
  
}

