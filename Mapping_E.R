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


##################################################################
###
### .bam mapping
###
##################################################################

setwd("/data/usrhome/LabCCKing/ccking01/Ko/0810CLCmapping/0810bam/")
ls0 = list.files(getwd())

ls129 = ls0[1:129]
ls129c = strsplit(ls129, split="_con", fixed = T)
ls129clean = sapply(ls129c, function(x) head(x,1))

ls92 = ls0[130:length(ls0)]
ls92c = strsplit(ls92, split= " ", fixed = T)
ls92clean = sapply(ls92c, function(x) head(x,1))

bamls = c(ls129clean, ls92clean)

bamseqprep <- function(bam){
  
  library(pasillaBamSubset)
  library(stringr)
  library(Rsamtools)
  
  r.bam <- scanBam(BamFile(bam))
  
  rmcigar = c("N", "P", "=", "X")           # discard 4
  # keep: D, I, S, H, M          
  lcigar = strsplit(r.bam[[1]]$cigar, split ="")
  rml = which(unlist(lapply(lcigar, function(x){any(rmcigar %in% x)})) == TRUE)
  
  if (length(rml) != 0){
    
    lcigar = as.list( r.bam[[1]]$cigar[-rml] )                # lcigar
    lseq = as.list( as.character(r.bam[[1]]$seq) )[-rml]      # lseq
    lpos = as.list( r.bam[[1]]$pos[-rml] )                    # lpos
    lss = as.list(r.bam[[1]]$strand[-rml])                    # sense
    
    ll = Map(c, lpos, lcigar, lss, lseq)
    
  }else{   
    
    lcigar = as.list( r.bam[[1]]$cigar )                # lcigar
    lseq = as.list( as.character(r.bam[[1]]$seq) )      # lseq
    lpos = as.list( r.bam[[1]]$pos )                    # lpos
    lss = as.list(r.bam[[1]]$strand)                    # sense
    
    
    ll = Map(c, lpos, lcigar, lss, lseq)
    
  }
  
  
  mm = do.call("rbind", lapply(ll, function(y){
    
    ## function(y)
    
    cigar.i = y[2]
    Scigar.i = strsplit(cigar.i, split="")[[1]]
    ss.i = y[3]
    
    seq.i = unlist(strsplit(y[4], split = ""))
    
    p.Chr = which(regexpr("[[:upper:]]+", Scigar.i) + 1 > 0)  # Upper Alphabet
    p.Num = which(regexpr("[[:digit:]]+", Scigar.i) + 1 > 0)  # Digit
    
    
    d.unit=c(p.Chr[1]-1)                      # length of each digit
    if(length(p.Chr) > 1){
      for(i in 2:(length(p.Chr))){
        d.unit[length(d.unit) + 1 ] = p.Chr[i] - p.Chr[i-1] -1 }}
    
    
    l.unit = c()                              # exract digit part
    for(i in 1: length(p.Chr)){ 
      l.unit[length(l.unit) + 1] = substr(cigar.i,p.Chr[i]-d.unit[i], p.Chr[i]-1)}
    l.unit=as.numeric(l.unit)
    
    Chri=Scigar.i[p.Chr]
    
    if(Chri[1] == "S"){
      pos.i = as.numeric(y[1]) - l.unit[1]
      Mt.i = matrix(c(rep("NA", pos.i-1), seq.i, rep("NA", 3000-length(seq.i)-pos.i+1)), 1)
      
    }else{
      pos.i = as.numeric(y[1])
      Mt.i = matrix(c(rep("NA", pos.i-1), seq.i, rep("NA", 3000-length(seq.i)-pos.i+1)), 1)
    }
    
    
    ID = c("I", "D")
    
    if (any(ID %in% Chri)){
      
      if (length(which(Chri == "I")) >= 1){
        
        irm = c()                    # deal with I
        
        for(h in 1: length(which(Chri == "I"))){
          
          l = which(Chri == "I")[h] - 1
          
          i.start = pos.i + sum(l.unit[0:l])
          i.stop = pos.i + sum(l.unit[0:l]) + l.unit[which(Chri == "I")[h]] - 1
          
          irm = c(irm, seq(i.start,i.stop)) }
        
        Mt.i = matrix(c(Mt.i[,-(irm)], 
                        rep("NA", sum( l.unit[which(Chri == "I")] ) )), 1) }
      
      
      if (length(which(Chri == "D")) >= 1){
        
        irm = c()                    # deal with D
        
        for(h in 1: length(which(Chri == "D"))){
          
          l = which(Chri == "D")[h] - 1
          
          i.start = pos.i + sum(l.unit[0:l])
          
          Mt.i = matrix(c(Mt.i[1:i.start-1], rep("NA", l.unit[which(Chri == "D")[h]]),
                          Mt.i[i.start:(3000- (l.unit[which(Chri == "D")[h]]))]) ,1) }
        
      }
      
      
      dMI = which(Chri != "M" & Chri != "I" & Chri != "D")
      
      if(length(dMI) >= 1){
        
        l.unit.dI = l.unit[-which(ID %in% Chri)]
        
        for(k in 1: length(dMI)){
          
          l = dMI[k] - 1 
          
          r.start = pos.i + sum(l.unit.dI[0:l])
          r.stop = pos.i + sum(l.unit.dI[0:l]) + l.unit[dMI[k]] - 1
          
          for(p in 1:length(r.start)){
            Mt.i[ , r.start[p]:r.stop[p]]= "NA" }}   }  
      
      
      if(ss.i == "1"){
        
        Pf<-c(884, 901, 1310, 1327, 1616, 1635, 2026, 2044)
        for(k in 1:4){
          
          kk = (k-1)*2+1
          Mt.i[ ,Pf[kk]:Pf[kk+1]]="NA"}
        
        return(Mt.i)}
      
      if(ss.i == "2"){
        
        Pr<-c(1357, 1376, 1752, 1770, 2035, 2053, 2504, 2524)
        for(k in 1:4){
          
          kk = (k-1)*2+1
          Mt.i[ ,Pr[kk]:Pr[kk+1]]="NA"}
        
        return(Mt.i)}
      
    } 
    
    else {
      
      
      dM = which(Chri != "M")
      
      if(length(dM) >= 1){
        
        for(k in 1: length(dM)){
          
          l = dM[k] - 1 
          
          r.start = pos.i + sum(l.unit[0:l])
          r.stop = pos.i + sum(l.unit[0:l]) + l.unit[dM[k]] - 1
          
          for(p in 1:length(r.start)){
            Mt.i[ , r.start[p]:r.stop[p]]= "NA"}}  
      }    
      
      if(ss.i == "1"){
        
        Pf<-c(884, 901, 1310, 1327, 1616, 1635, 2026, 2044)
        for(k in 1:4){
          
          kk = (k-1)*2+1
          Mt.i[ ,Pf[kk]:Pf[kk+1]]="NA"}
        
        return(Mt.i)}
      
      if(ss.i == "2"){
        
        Pr<-c(1357, 1376, 1752, 1770, 2035, 2053, 2504, 2524)
        for(k in 1:4){
          
          kk = (k-1)*2+1
          Mt.i[ ,Pr[kk]:Pr[kk+1]]="NA"}
        
        return(Mt.i)}
      
    } 
    
    
  }))
  
}

for (i in 1: length(list.files(getwd()))){
  
  t1 = Sys.time()
  
  Mt = bamseqprep(ls0[i])
  write.csv(Mt, 
            paste0("/data/usrhome/LabCCKing/ccking01/Desktop/Ko_DENV/bam/mapped/", bamls[i], ".csv"))
  
  print(bamls[i])
  
  t2 = Sys.time() - t1
  print(t2)
  
}











