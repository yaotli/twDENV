##### mapping from .bam files ######
# use fnuction: bamseqprep or bamseqprepID in Function.R
#

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
  rml2 = which(r.bam[[1]]$pos > 3000)
  rml = unique(c(rml, rml2))
  
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
      
      if (length(which(Chri == "D")) >= 1){
        
        irm = c()                    # deal with D
        
        for(h in 1: length(which(Chri == "D"))){
          
          l = which(Chri == "D")[h] - 1
          
          i.start = pos.i + sum(l.unit[0:l])
          
          Mt.i = matrix(c(Mt.i[1:i.start-1], rep("NA", l.unit[which(Chri == "D")[h]]),
                          Mt.i[i.start:(3000- (l.unit[which(Chri == "D")[h]]))]) ,1) }
        
      }
      
      if (length(which(Chri == "I")) >= 1){
        
        irm = c()                    # deal with I
        
        for(h in 1: length(which(Chri == "I"))){
          
          l = which(Chri == "I")[h] - 1
          
          i.start = pos.i + sum(l.unit[0:l])
          i.stop = pos.i + sum(l.unit[0:l]) + l.unit[which(Chri == "I")[h]] - 1
          
          irm = c(irm, seq(i.start,i.stop)) }
        
        Mt.i = matrix(c(Mt.i[,-(irm)], 
                        rep("NA", sum( l.unit[which(Chri == "I")] ) )), 1) }
      
      Chri.dI = Chri[-(which(Chri == "I"))]
      dMI = which( Chri.dI != "M" & Chri.dI != "D" )
      
      if( length(dMI) >= 1 ){
        
        l.unit.dI = l.unit[-which(Chri %in% "I" == TRUE)]
        
        for(k in 1: length(dMI)){
          
          l = dMI[k] - 1 
          
          r.start = pos.i + sum(l.unit.dI[0:l])
          r.stop = pos.i + sum(l.unit.dI[0:l]) + l.unit[dMI[k]] - 1
          
          for(p in 1:length(r.start)){
            Mt.i[ , r.start[p]:r.stop[p]]= "NA" } }   }  
      
      
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
            Mt.i[ , r.start[p]:r.stop[p]]= "NA"} }  }    
      
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

