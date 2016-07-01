#Fuctions to detect minor variant

#1. Proportion of Nucleotide (pon)
pon<-function(csvfile, nucleotide, position){
  
  #ddf = read.csv(files)
  #nucleotide: A=1, T=2, C=3, G=4, N=5
  #position = real nucleotide position of the polyproteion of 2001 DENV
  
  ddf <- read.csv(csvfile)
  xpo = (position)+2
  totalcount<-sum(ddf[,xpo][1], ddf[,xpo][2], ddf[,xpo][3], ddf[,xpo][4],
                  ddf[,xpo][5])
  
  Npo = ddf[,xpo][nucleotide]/totalcount
  return(Npo)
  
}

#2. ATCG specific pon
ponA<-function(csvfile, position){
  
  
  #ddf = read.csv(files)
  #nucleotide = A
  #position = real nucleotide position of the polyproteion of 2001 DENV
  
  ddf <- read.csv(csvfile)
  xpo = (position)+2
  totalcount<-sum(ddf[,xpo][1], ddf[,xpo][2], ddf[,xpo][3], ddf[,xpo][4],
                  ddf[,xpo][5])
  
  Npo = ddf[,xpo][1]/totalcount
  return(Npo)
  
}
ponT<-function(csvfile, position){
  
  
  #ddf = read.csv(files)
  #nucleotide = T
  #position = real nucleotide position of the polyproteion of 2001 DENV
  
  ddf <- read.csv(csvfile)
  xpo = (position)+2
  totalcount<-sum(ddf[,xpo][1], ddf[,xpo][2], ddf[,xpo][3], ddf[,xpo][4],
                  ddf[,xpo][5])
  
  Npo = ddf[,xpo][2]/totalcount
  return(Npo)
  
}
ponC<-function(csvfile, position){
  
  
  #ddf = read.csv(files)
  #nucleotide = A
  #position = real nucleotide position of the polyproteion of 2001 DENV
  
  ddf <- read.csv(csvfile)
  xpo = (position)+2
  totalcount<-sum(ddf[,xpo][1], ddf[,xpo][2], ddf[,xpo][3], ddf[,xpo][4],
                  ddf[,xpo][5])
  
  Npo = ddf[,xpo][3]/totalcount
  return(Npo)
  
}
ponG<-function(csvfile, position){
  
  
  #ddf = read.csv(files)
  #nucleotide = G
  #position = real nucleotide position of the polyproteion of 2001 DENV
  
  ddf <- read.csv(csvfile)
  xpo = (position)+2
  totalcount<-sum(ddf[,xpo][1], ddf[,xpo][2], ddf[,xpo][3], ddf[,xpo][4],
                  ddf[,xpo][5])
  
  Npo = ddf[,xpo][4]/totalcount
  return(Npo)
  
}

#3. Vector-Directed pon (with multiple samples' percentage)
#require pon (*4)
vdpon<-function(site, nameoflist){
  
  pA=c()
  pT=c()
  pC=c()
  pG=c()
  
  for (i in 1: length(nameoflist)){
    
    pA[i] = ponA(nameoflist[i], site)
    pT[i] = ponT(nameoflist[i], site)
    pC[i] = ponC(nameoflist[i], site)
    pG[i] = ponG(nameoflist[i], site)
    
  }
  rrm <- cbind(pA, pT, pC, pG)
  return(rrm)
  
  
}



#Detection of variant (dvariant)

dvariant<-function(ddf, pc){
  
  VarN<-c()
  refN<-c()
  pVarN<-c()
  sVarN<-c()
  
  dddf<-read.csv(ddf)
  lth = length(dddf)
  
  for(i in 3: lth){
    
    if ((dddf[1,i] | dddf[2,i] | dddf[3,i] | dddf[4,i] | dddf[5,i]) != 0 ) {
      
      tcc = sum(dddf[1,i] + dddf[2,i] + dddf[3,i] + dddf[4,i])  #total counts
      tccpc = tcc*pc  #least counts we are interested in
      tccpcp = tcc*(1-pc) #maximum count of a putative consensus when minor population exits
      
      mm = max(dddf[1,i], dddf[2,i], dddf[3,i], dddf[4,i]) #maximum counts (count of consensus seq)
      
      if (mm < tccpcp){
        
        aamatrix=c("A","T","C","G")
        ii = i-2  
        
        if ((dddf[1,i] >= tccpc) & (dddf[1,i] != mm)){
          (VarN[length(VarN)+1] = "A")
          (pVarN[length(pVarN)+1] = (dddf[1,i]/tcc))
          (sVarN[length(sVarN)+1] = (ii))
          (refN[length(refN)+1] = aamatrix[which.max(mm)])}
        
        if ((dddf[2,i] >= tccpc) & (dddf[2,i] != mm)){
          (VarN[length(VarN)+1] = "T")
          (pVarN[length(pVarN)+1] = (dddf[2,i]/tcc))
          (sVarN[length(sVarN)+1] = (ii))
          (refN[length(refN)+1] = aamatrix[which.max(mm)])}
        
        if ((dddf[3,i] >= tccpc) & (dddf[3,i] != mm)){
          (VarN[length(VarN)+1] = "C")
          (pVarN[length(pVarN)+1] = (dddf[3,i]/tcc))
          (sVarN[length(sVarN)+1] = (ii))
          (refN[length(refN)+1] = aamatrix[which.max(mm)])}
        
        if ((dddf[4,i] >= tccpc) & (dddf[4,i] != mm)){
          (VarN[length(VarN)+1] = "G")
          (pVarN[length(pVarN)+1] = (dddf[4,i]/tcc))
          (sVarN[length(sVarN)+1] = (ii))
          (refN[length(refN)+1] = aamatrix[which.max(mm)])}
      }
      
      print(i/lth*100)
    }
  }
  
  rb<-data.frame(sVarN, VarN, pVarN)
  return(rb)
}
