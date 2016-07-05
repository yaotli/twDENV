#Fuctions to detect minor variant

##1. Proportion of Nucleotide (pon)
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


#1.1 ATCG specific pon
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


#1.2 Vector-Directed pon (with multiple samples' percentage)
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


##2 Detection of variant (dvariant) 
#ddf = a csvfile name
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
        nmax = which.max(c(dddf[1,i], dddf[2,i], dddf[3,i], dddf[4,i]))
        
        ii = i-2  
        
        if ((dddf[1,i] >= tccpc) & (dddf[1,i] != mm)){
          (VarN[length(VarN)+1] = "A")
          (pVarN[length(pVarN)+1] = (dddf[1,i]/tcc))
          (sVarN[length(sVarN)+1] = (ii))
          (refN[length(refN)+1] = aamatrix[nmax])}
        
        if ((dddf[2,i] >= tccpc) & (dddf[2,i] != mm)){
          (VarN[length(VarN)+1] = "T")
          (pVarN[length(pVarN)+1] = (dddf[2,i]/tcc))
          (sVarN[length(sVarN)+1] = (ii))
          (refN[length(refN)+1] = aamatrix[nmax])}
        
        if ((dddf[3,i] >= tccpc) & (dddf[3,i] != mm)){
          (VarN[length(VarN)+1] = "C")
          (pVarN[length(pVarN)+1] = (dddf[3,i]/tcc))
          (sVarN[length(sVarN)+1] = (ii))
          (refN[length(refN)+1] = aamatrix[nmax])}
        
        if ((dddf[4,i] >= tccpc) & (dddf[4,i] != mm)){
          (VarN[length(VarN)+1] = "G")
          (pVarN[length(pVarN)+1] = (dddf[4,i]/tcc))
          (sVarN[length(sVarN)+1] = (ii))
          (refN[length(refN)+1] = aamatrix[nmax])}
      }
    }
    
  }
  
  rb<-data.frame(sVarN, VarN, refN, pVarN)
  return(rb)
}

#2.2 Detection of variant (dvariant) with reliable reads number
dvariant1000r<-function(ddf, pc){
  
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
      
      if ((mm < tccpcp) & (tcc > 1000)){
        
        aamatrix=c("A","T","C","G")
        nmax = which.max(c(dddf[1,i], dddf[2,i], dddf[3,i], dddf[4,i]))
        
        ii = i-2  
        
        if ((dddf[1,i] >= tccpc) & (dddf[1,i] != mm)){
          (VarN[length(VarN)+1] = "A")
          (pVarN[length(pVarN)+1] = (dddf[1,i]/tcc))
          (sVarN[length(sVarN)+1] = (ii))
          (refN[length(refN)+1] = aamatrix[nmax])}
        
        if ((dddf[2,i] >= tccpc) & (dddf[2,i] != mm)){
          (VarN[length(VarN)+1] = "T")
          (pVarN[length(pVarN)+1] = (dddf[2,i]/tcc))
          (sVarN[length(sVarN)+1] = (ii))
          (refN[length(refN)+1] = aamatrix[nmax])}
        
        if ((dddf[3,i] >= tccpc) & (dddf[3,i] != mm)){
          (VarN[length(VarN)+1] = "C")
          (pVarN[length(pVarN)+1] = (dddf[3,i]/tcc))
          (sVarN[length(sVarN)+1] = (ii))
          (refN[length(refN)+1] = aamatrix[nmax])}
        
        if ((dddf[4,i] >= tccpc) & (dddf[4,i] != mm)){
          (VarN[length(VarN)+1] = "G")
          (pVarN[length(pVarN)+1] = (dddf[4,i]/tcc))
          (sVarN[length(sVarN)+1] = (ii))
          (refN[length(refN)+1] = aamatrix[nmax])}
      }
    }
    
  }
  
  rb<-data.frame(sVarN, VarN, refN, pVarN)
  return(rb)
}

dvariant2000r<-function(ddf, pc){
  
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
      
      if ((mm < tccpcp) & (tcc > 2000)){
        
        aamatrix=c("A","T","C","G")
        nmax = which.max(c(dddf[1,i], dddf[2,i], dddf[3,i], dddf[4,i]))
        
        ii = i-2  
        
        if ((dddf[1,i] >= tccpc) & (dddf[1,i] != mm)){
          (VarN[length(VarN)+1] = "A")
          (pVarN[length(pVarN)+1] = (dddf[1,i]/tcc))
          (sVarN[length(sVarN)+1] = (ii))
          (refN[length(refN)+1] = aamatrix[nmax])}
        
        if ((dddf[2,i] >= tccpc) & (dddf[2,i] != mm)){
          (VarN[length(VarN)+1] = "T")
          (pVarN[length(pVarN)+1] = (dddf[2,i]/tcc))
          (sVarN[length(sVarN)+1] = (ii))
          (refN[length(refN)+1] = aamatrix[nmax])}
        
        if ((dddf[3,i] >= tccpc) & (dddf[3,i] != mm)){
          (VarN[length(VarN)+1] = "C")
          (pVarN[length(pVarN)+1] = (dddf[3,i]/tcc))
          (sVarN[length(sVarN)+1] = (ii))
          (refN[length(refN)+1] = aamatrix[nmax])}
        
        if ((dddf[4,i] >= tccpc) & (dddf[4,i] != mm)){
          (VarN[length(VarN)+1] = "G")
          (pVarN[length(pVarN)+1] = (dddf[4,i]/tcc))
          (sVarN[length(sVarN)+1] = (ii))
          (refN[length(refN)+1] = aamatrix[nmax])}
      }
    }
    
  }
  
  rb<-data.frame(sVarN, VarN, refN, pVarN)
  return(rb)
}

#2.3 Group detection of dvariant
#ls = list of csvfile
#pc = percentage of curring poing
#nr = number of reads as threshold: 0 = origin, 1 = 1000, 2 = 2000

lsvariant<-function(ls, pc, nr){

#dvariant
  if( nr  == 0 ){    
    for (i in 1:length(ls)){
    
    assign(paste0("sdf",".",i), dvariant(ls[i], pc)) #dataframe of sample
    print((i/length(ls)*100))
  }
    }

#dvariant1000  
  if( nr == 1  ){    
    for (i in 1:length(ls)){
      
      assign(paste0("sdf",".",i), dvariant1000r(ls[i], pc)) #dataframe of sample
      print((i/length(ls)*100))
    }
  }
  
#dvariant2000  
  if( nr == 2   ){    
    for (i in 1:length(ls)){
      
      assign(paste0("sdf",".",i), dvariant2000r(ls[i], pc)) #dataframe of sample
      print((i/length(ls)*100))
    }
  }
  
      
    plist<-mget(paste0("sdf.",1:length(ls)))
    lldvar<-do.call("rbind", plist)
    
  return(lldvar)
}

#3 nonsynonymous mutation 
#dvdf = dataframe generated by #4 dvariant 
#ls = csvfile list
nsvariant<-function(dvdf, ls){
  
  library(seqinr)
  nsv<-c()
  
  for(i in 1:dim(dvdf)[1]){
    
    sVarN = dvdf[[1]][i]
    
    if(sVarN %%3 == 2 ){
      
      nsv[length(nsv)+1] = "TRUE"
      
    }
    
    if(sVarN %%3 == 1 ){
      
      n = strsplit(attributes(dvdf)$row.names, split = ".", fixed = T)[[i]][2] #number in list
      fls <- read.csv(ls[as.numeric(n)])
      
      aamatrix=c("A","T","C","G")
      
      p = sVarN + 2
      
      pN = which.max(fls[[p]])
      p1N = which.max(fls[[p+1]])
      p2N = which.max(fls[[p+2]])
      
      refaa<-c(aamatrix[pN], aamatrix[p1N], aamatrix[p2N])
      varaa<-c(as.character(dvdf[[2]][i]), aamatrix[p1N], aamatrix[p2N])
      
      if(translate(refaa) == translate(varaa)){ nsv[length(nsv)+1] = "FALSE" } else 
      { nsv[length(nsv)+1] = "TRUE" }
      
    }
      
    if(sVarN %%3 == 0 ){
      
      n = strsplit(attributes(dvdf)$row.names, split = ".", fixed = T)[[i]][2] #number in list
      fls <- read.csv(ls[as.numeric(n)])
      
      aamatrix=c("A","T","C","G")
      
      p = sVarN + 2
      
      pm2N = which.max(fls[[p-2]])
      pm1N = which.max(fls[[p-1]])
      pN = which.max(fls[[p]])
      
      refaa<-c(aamatrix[pm2N], aamatrix[pm1N], aamatrix[pN])
      varaa<-c(aamatrix[pm2N], aamatrix[pm1N], aamatrix[dvdf[[2]][i]])
       
      if(translate(refaa) == translate(varaa)){ nsv[length(nsv)+1] = "FALSE" } else 
                                              { nsv[length(nsv)+1] = "TRUE" }

      
    }
    
print((i/dim(dvdf)[1])*100)
    
  }
  
    
  
  return(nsv)
  
}

