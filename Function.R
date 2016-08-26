#Fuctions to detect minor variant

##1. Proportion of Nucleotide (pon) --------------------------------
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


##2 Detection of variant (dvariant, lsvariant) --------------------------------
#ddf = a csvfile name
dvariant<-function(ddf, pc, nr){
  
  VarN<-c()
  refN<-c()
  pVarN<-c()
  sVarN<-c()
  nrv=c(0, 500, 1000, 2000, 2500) #1 = 0; 2 = 500; 3 = 1000; 4 = 2000; 5 = 2500
  
  dddf<-read.csv(ddf)
  lth = length(dddf)
  
  for(i in 3: lth){
    
    if ((dddf[1,i] | dddf[2,i] | dddf[3,i] | dddf[4,i] | dddf[5,i]) != 0 ) {
      
      tcc = sum(dddf[1,i] + dddf[2,i] + dddf[3,i] + dddf[4,i])  #total counts
      tccpc = tcc*pc  #least counts we are interested in
      tccpcp = tcc*(1-pc) #maximum count of a putative consensus when minor population exits
      
      mm = max(dddf[1,i], dddf[2,i], dddf[3,i], dddf[4,i]) #maximum counts (count of consensus seq)
      
      if ((mm <= tccpcp) & (tcc >= nrv[nr])){
        
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


#2.2 Group detection of dvariant (lsvariant)
#ls = list of csvfile
#pc = percentage of curring poing
#nr = number of reads as threshold: 1 = 0; 2 = 500; 3 = 1000; 4 = 2000; 5 = 2500

lsvariant<-function(ls, pc, nr){

    for (i in 1:length(ls)){
    
    assign(paste0("sdf",".",i), dvariant(ls[i], pc, nr)) #dataframe of sample
    print((i/length(ls)*100))
  }
  
    plist<-mget(paste0("sdf.",1:length(ls)))
    lldvar<-do.call("rbind", plist)
    
  return(lldvar)
}


#2.3 Poisson test (PT) version of dvariant

variantPT <- function(df, pc, nr, adjER){
  
  library(dplyr)
  
  # need adjER from ErrorRate.R
  
  nrv = c(0, 500, 1000, 2000, 2500)
  mAA = c("A", "T", "C", "G")
  
  sVarN <- c()
  VarN  <- c()
  refN  <- c()
  pVarN <- c()
  
  for (i in 3: length(df)){
    
    ii = i - 2 
    
    vc = as.vector(df[,i])
    vc = data.frame(vc, c(1:5))
    names(vc) = c("count", "AA")
    
    tcc = sum(vc[,1][1:4])                         # nr
    
    if ( tcc >= nrv[nr] ){
      
      mm = max(vc[,1])
      tcc0 = sum(vc[,1])                         
      pp = qpois(pc, tcc0*adjER[ii], lower.tail = FALSE)
      
      
      vci =   
        vc %>%
        filter(count != mm) %>%
        filter(count > pp )  
      
      if (nrow(vci) > 0){ 
        
        VarNi = mAA[vci[,2]]
        refNi = rep(mAA[vc[,2][which.max(vc[,1])]], nrow(vci))
        pVarNi = c(vci[,1]/tcc)
        
        sVarN = c(sVarN, rep(ii, nrow(vci)))
        VarN  = c(VarN, VarNi)
        refN  = c(refN, refNi)
        pVarN = c(pVarN, pVarNi)
        
      }        
      
    } 
    
  }
  
  rb<-data.frame(sVarN, VarN, refN, pVarN)
  return(rb)
  
  
}

#2.4 ls of variantPT

lsvariantPT <- function(df, pc, nr, adjER, listname){
                                                # data source: cSUM
  
  LSi = split(df, rep(1:length(listname), each = 5))
  
  VR = do.call("rbind" , lapply(LSi, variantPT, pc=pc, nr=nr, adjER=adjER))
  
  n0 = attributes(VR)$row.names
  n1 = strsplit(n0, split=".", fixed=T)
  n1 = as.numeric(sapply(n1, function(x) head(x,1)))
  tt = as.vector(table(cut(n1, breaks=c(1:(length(listname)+1)), right=FALSE)))
  
  dd = data.frame(tt, listname)
  nn = do.call("c", apply(dd, 1, function(x){ 
    
    time = as.numeric(x[1])
    name = as.character(x[2])
    rr = rep(name, time)
    return(rr)      }) )
  
  
  VRn = split(VR, nn)
  
  return(VRn)
  print( listname[which(tt == 0)] )
}


##3 nonsynonymous mutation (nsvariant, nsvariantsite, cdvariant)  --------------------------
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
      varaa<-c(aamatrix[pm2N], aamatrix[pm1N], as.character(dvdf[[2]][i]))
       
      if(translate(refaa) == translate(varaa)){ nsv[length(nsv)+1] = "FALSE" } else 
                                              { nsv[length(nsv)+1] = "TRUE" }

      
    }
    
print((i/dim(dvdf)[1])*100)
    
  }
  
    
  
  return(nsv)
  
}

#3.2
nsvariantsite<-function(df, site, wt, mu){ #wt, mu = 1(A), 2(T), 3(C), 4(G)
  library(seqinr)
  
  ssite = as.numeric(site)
  PoC = ssite %%3 
  
  if(PoC == 2){return("TRUE")} else {
    
    ddf <- read.csv(df)
    aamatrix=c("A","T","C","G")
    p = ssite + 2
    
    if(PoC == 1){
      refaa<-c(aamatrix[as.numeric(wt)], aamatrix[which.max(ddf[[p+1]])], 
               aamatrix[which.max(ddf[[p+2]])])
      varaa<-c(aamatrix[as.numeric(mu)], aamatrix[which.max(ddf[[p+1]])], 
               aamatrix[which.max(ddf[[p+2]])])
      
      if(translate(refaa) == translate(varaa)){ return("FALSE") } else 
      { return("TRUE") }
    }
    
    if(PoC == 0){
      refaa<-c(aamatrix[which.max(ddf[[p-2]])], aamatrix[which.max(ddf[[p-1]])], 
               aamatrix[as.numeric(wt)])
      varaa<-c(aamatrix[which.max(ddf[[p-2]])], aamatrix[which.max(ddf[[p-1]])], 
               aamatrix[as.numeric(mu)])
      
      if(translate(refaa) == translate(varaa)){ return("FALSE") } else 
      { return("TRUE") }
    }
    
    
  }
  
  
}


#3.3 combined 2.3 + 3

lsvariantNS<-function(ls, pc, nr){
  
    for (i in 1:length(ls)){
      
      assign(paste0("sdf",".",i), dvariant(ls[i], pc, nr)) #dataframe of sample
      print((i/length(ls)*100))
    }
  
  plist<-mget(paste0("sdf.",1:length(ls)))
  lldvar<-do.call("rbind", plist)
  
  lsns<-nsvariant(lldvar, ls)
  ffls<-cbind(lldvar, lsns)
  
  return(ffls)
}


#3.4 count only (for individual file)

cdvariant<-function(ddf, pc, nr){
  
  NpV = 0 #number of "position with variant"
  NV = 0  #"number of variant" in all positions
  
  NpVa = 0 #number of nonsynonymous "site" in NpV
  NVa = 0 #number of nonsynonymous "variant" in NV
  
  nrv=c(0, 500, 1000, 2000, 2500) #1 = 0; 2 = 500; 3 = 1000; 4 = 2000; 5 = 2500
  
  dddf<-read.csv(ddf)
  lth = length(dddf)
  
  for(i in 3: lth){ #each position
    
    sitecount = 0
    sitecountaa = 0
    
    if ((dddf[1,i] | dddf[2,i] | dddf[3,i] | dddf[4,i] | dddf[5,i]) != 0 ) {
      
      tcc = sum(dddf[1,i] + dddf[2,i] + dddf[3,i] + dddf[4,i])  #total counts
      tccpc = tcc*pc  #least counts we are interested in
      tccpcp = tcc*(1-pc) #maximum count of a putative consensus when minor population exits
      
      mm = max(dddf[1,i], dddf[2,i], dddf[3,i], dddf[4,i]) #maximum counts (count of consensus seq)
      
      if ((mm <= tccpcp) & (tcc >= nrv[nr])){
        
        nmax = which.max(c(dddf[1,i], dddf[2,i], dddf[3,i], dddf[4,i]))
        ii = i-2  
        
        if ((dddf[1,i] >= tccpc) & (dddf[1,i] != mm)){
          
          sitecount = sitecount + 1
          if ((nsvariantsite(ddf,ii ,nmax ,1)) == "TRUE"){ sitecountaa = sitecountaa + 1}
          }
        
        if ((dddf[2,i] >= tccpc) & (dddf[2,i] != mm)){
          
          sitecount = sitecount + 1
          if ((nsvariantsite(ddf,ii ,nmax ,2)) == "TRUE"){ sitecountaa = sitecountaa + 1}
        }
        
        if ((dddf[3,i] >= tccpc) & (dddf[3,i] != mm)){
          
          sitecount = sitecount + 1
          if ((nsvariantsite(ddf,ii ,nmax ,3)) == "TRUE"){ sitecountaa = sitecountaa + 1}
        }
        
        if ((dddf[4,i] >= tccpc) & (dddf[4,i] != mm)){
          
          sitecount = sitecount + 1
          if ((nsvariantsite(ddf,ii ,nmax ,4)) == "TRUE"){ sitecountaa = sitecountaa + 1}
        }
    }
    
      
  }
  
    
    if (sitecount != 0){ NpV = NpV + 1}
    if (sitecountaa != 0){ NpVa = NpVa + 1}
    
    NV = sum(sitecount) + NV
    NVa = sum(sitecountaa) + NVa
  }
  
  re = c(NpV, NV, NpVa, NVa)
  return(re)
}


#3.5 position recognition (for Heatmap)

pvariant<-function(ls, b, pc, nr, kd){ #kd = kind: 1 = position, 2 = all variant, 
                                       #           3 = position of NS, 4 = variant of NS
                                       # b = breaks, ex: breaks = seq(930, 2430, by=30)
  mx = c()
  nrv=c(0, 500, 1000, 2000, 2500) #1 = 0; 2 = 500; 3 = 1000; 4 = 2000; 5 = 2500
  
  for (k in 1:length(ls)){  #each file
    
    dddf<-read.csv(ls[k])
    lth = length(dddf)
    sitecount = c()
    sitecountaa = c()
    
    
    for(i in 3: lth){ #each position
      
      if ((dddf[1,i] | dddf[2,i] | dddf[3,i] | dddf[4,i] | dddf[5,i]) != 0 ) { 
        
        tcc = sum(dddf[1,i] + dddf[2,i] + dddf[3,i] + dddf[4,i])  #total counts
        tccpc = tcc*pc  #least counts we are interested in
        tccpcp = tcc*(1-pc) #maximum count of a putative consensus when minor population exits
        
        mm = max(dddf[1,i], dddf[2,i], dddf[3,i], dddf[4,i]) #maximum counts (count of consensus seq)
        
        if ((mm <= tccpcp) & (tcc >= nrv[nr])){
        
          nmax = which.max(c(dddf[1,i], dddf[2,i], dddf[3,i], dddf[4,i]))
          ii = i-2
           
            if ((dddf[1,i] >= tccpc) & (dddf[1,i] != mm)){
                
              sitecount[length(sitecount) +1 ] = ii
              
                  if (((nsvariantsite(ls[k],ii,nmax,1)) == "TRUE") & (kd > 2)){ 
                
                    sitecountaa[length(sitecountaa) +1 ] = ii
              } } #A
            
            if ((dddf[2,i] >= tccpc) & (dddf[2,i] != mm)){
                
                sitecount[length(sitecount) +1 ] = ii
                
                if (((nsvariantsite(ls[k],ii,nmax,2)) == "TRUE") & (kd > 2)){ 
                  
                  sitecountaa[length(sitecountaa) +1 ] = ii
                } } #T
              
            if ((dddf[3,i] >= tccpc) & (dddf[3,i] != mm)){
                  
                  sitecount[length(sitecount) +1 ] = ii
                  
                if (((nsvariantsite(ls[k],ii,nmax,3)) == "TRUE") & (kd > 2)){ 
                    
                    sitecountaa[length(sitecountaa) +1 ] = ii
                  } } #C
             
            if ((dddf[4,i] >= tccpc) & (dddf[4,i] != mm)){
                    
                    sitecount[length(sitecount) +1 ] = ii
                    
                 if (((nsvariantsite(ls[k],ii,nmax,4)) == "TRUE") & (kd > 2)){ 
                      
                      sitecountaa[length(sitecountaa) +1 ] = ii
                    }   } #G       
        
           
          
        }

        
      }
      
      
 
      
      }

    if(kd == 1){
      sitecount=unique(sitecount)
      region<-cut(sitecount, b, right = FALSE)
      region<-data.frame(k, table(region))
      
      mx=rbind(mx,region) }
    
    if(kd == 2){
      region<-cut(sitecount, b, right = FALSE)
      region<-data.frame(k, table(region))
      
      mx=rbind(mx,region) }
    
    if(kd == 3){
      sitecountaa=unique(sitecountaa)
      region<-cut(sitecountaa, b, right = FALSE)
      region<-data.frame(k, table(region))
      
      mx=rbind(mx,region)}
    
    if(kd == 4){      
      region<-cut(sitecountaa, b, right = FALSE)
      region<-data.frame(k, table(region))
    
      mx=rbind(mx,region) }  

    print((k/length(ls))*100)
    
  }


return(mx)    

}


#3.6 modified position recognition 


pvariantM<-function(ls, pc, nr, kd){ #kd = kind: 1 = position, 2 = all variant, 
  #           3 = position of NS, 4 = variant of NS
  # b = breaks, ex: breaks = seq(930, 2430, by=30)
  mx = c()
  mxk = c()
  nrv=c(0, 500, 1000, 2000, 2500) #1 = 0; 2 = 500; 3 = 1000; 4 = 2000; 5 = 2500
  
  for (k in 1:length(ls)){  #each file
    
    dddf<-read.csv(ls[k])
    lth = length(dddf)
    sitecount = c()
    sitecountaa = c()
    
    
    for(i in 3: lth){ #each position
      
      if ((dddf[1,i] | dddf[2,i] | dddf[3,i] | dddf[4,i] | dddf[5,i]) != 0 ) { 
        
        tcc = sum(dddf[1,i] + dddf[2,i] + dddf[3,i] + dddf[4,i])  #total counts
        tccpc = tcc*pc  #least counts we are interested in
        tccpcp = tcc*(1-pc) #maximum count of a putative consensus when minor population exits
        
        mm = max(dddf[1,i], dddf[2,i], dddf[3,i], dddf[4,i]) #maximum counts (count of consensus seq)
        
        if ((mm <= tccpcp) & (tcc >= nrv[nr])){
          
          nmax = which.max(c(dddf[1,i], dddf[2,i], dddf[3,i], dddf[4,i]))
          ii = i-2
          
          if ((dddf[1,i] >= tccpc) & (dddf[1,i] != mm)){
            
            sitecount[length(sitecount) +1 ] = ii
            
            if (((nsvariantsite(ls[k],ii,nmax,1)) == "TRUE") & (kd > 2)){ 
              
              sitecountaa[length(sitecountaa) +1 ] = ii
            } } #A
          
          if ((dddf[2,i] >= tccpc) & (dddf[2,i] != mm)){
            
            sitecount[length(sitecount) +1 ] = ii
            
            if (((nsvariantsite(ls[k],ii,nmax,2)) == "TRUE") & (kd > 2)){ 
              
              sitecountaa[length(sitecountaa) +1 ] = ii
            } } #T
          
          if ((dddf[3,i] >= tccpc) & (dddf[3,i] != mm)){
            
            sitecount[length(sitecount) +1 ] = ii
            
            if (((nsvariantsite(ls[k],ii,nmax,3)) == "TRUE") & (kd > 2)){ 
              
              sitecountaa[length(sitecountaa) +1 ] = ii
            } } #C
          
          if ((dddf[4,i] >= tccpc) & (dddf[4,i] != mm)){
            
            sitecount[length(sitecount) +1 ] = ii
            
            if (((nsvariantsite(ls[k],ii,nmax,4)) == "TRUE") & (kd > 2)){ 
              
              sitecountaa[length(sitecountaa) +1 ] = ii
            }   } #G       
          
          
          
        }
        
        
      }
      
      
      
      
    }
    
    if(kd == 1){
      sitecount=unique(sitecount)
      
      kk = c(rep(k, length(sitecount)))
      mx = c(mx, sitecount)
      mxk = c(mxk, kk)
      
      }
    
    if(kd == 2){
      
      kk = c(rep(k, length(sitecount)))
      mx = c(mx, sitecount)
      mxk = c(mxk, kk)
      
       }
    
    if(kd == 3){
      sitecountaa=unique(sitecountaa)

      kk = c(rep(k, length(sitecountaa)))
      mx = c(mx, sitecountaa)
      mxk = c(mxk, kk)
      
      }
    
    if(kd == 4){    
      
      kk = c(rep(k, length(sitecountaa)))
      mx = c(mx, sitecountaa)
      mxk = c(mxk, kk)
       }  
    
    print((k/length(ls))*100)
    
  }
  
  mxx=data.frame(mx,mxk)
  
  return(mxx)    
  
}


#3.7 position recognition fo AA (& for Heatmap)

pvariantAA.E<-function(ls, pc, nr, kd){ #kd = kind: 1 = position, 2 = all variant, 
  
  mx = c()
  nrv=c(0, 500, 1000, 2000, 2500) #1 = 0; 2 = 500; 3 = 1000; 4 = 2000; 5 = 2500
  
  library(dplyr)
  
  for (k in 1:length(ls)){  #each file
    
    dddf<-read.csv(ls[k])
    lth = length(dddf)
    
    # E protein AA region: 313 - 807
    
    ia = 313 + 1
    ib = 807 + 1
    saa = c(1:23)[-(19:21)]           #19- 21 : "*", "-", "X"
    
    for(i in ia: ib){                 #each position
      
      vc = as.vector(dddf[,i])
      vc = data.frame(vc)
      
      tcc = sum(vc[,1][saa]) # nr
      position = i - 313
      
      if ( tcc >= nrv[nr] ){
        
        mm = max(vc[,1][saa])
        
        NOr =      
          vc %>%
          filter(vc != mm) %>%
          filter(vc >= mm*pc) %>%
          nrow()        
        
        if (kd == 1){
          
          if (NOr > 0){ variantaa = 1 } else { variantaa = 0 }
          
          region = data.frame(k, position, variantaa, tcc)
          mx = rbind(mx, region)
          
          
        }else{
          
          variantaa = NOr
          
          region = data.frame(k, position, variantaa, tcc)
          mx = rbind(mx, region)
          
        }
        
        
      }else { variantaa = 0 
      
      region = data.frame(k, position, variantaa, tcc)
      mx = rbind(mx, region)}
      
      
    }     
    print((k/length(ls)*100)) 
    
  }      
  
  return(mx)
}  


#3.8 Poisson test for pvariant

pvariantPT<-function(ls, b, pc, nr, kd){ #kd = kind: 1 = position, 2 = all variant, 
  
  # b = breaks, ex: breaks = seq(930, 2430, by=30)
  # pc = p value ex: 0.05, 0.01, 0.001
  
  mx  = c()
  nrv = c(0, 500, 1000, 2000, 2500)      #1 = 0; 2 = 500; 3 = 1000; 4 = 2000; 5 = 2500
  
  library(dplyr)
  
  for (k in 1:length(ls)){  #each file
    
    variantc = c()
    dddf<-read.csv(ls[k])
    lth = length(dddf)
    
    for(i in 3: lth){                           #each position
      
      ii = i - 2 
      
      vc = as.vector(dddf[,i])
      vc = data.frame(vc)
      
      tcc = sum(vc[,1][1:4])                         # nr
      
      if ( tcc >= nrv[nr] ){
        
        mm = max(vc[,1])
        tcc0 = sum(vc[,1])                         
        pp = qpois(pc, tcc0*0.0113, lower.tail = FALSE)
        
        
        NOr =      
          vc %>%
          filter(vc != mm) %>%
          filter(vc > pp ) %>%
          nrow()        
        
        if (NOr > 0){ variantc = c(variantc, rep(ii, NOr)) }        
      } 
      
    }
    
    if (kd > 0){
      
      variantc = unique(variantc)
      
      region<-cut(variantc, b, right = FALSE)
      region<-data.frame(k, table(region))
      mx = rbind(mx, region)    } else {
        
        region<-cut(variantc, b, right = FALSE)
        region<-data.frame(k, table(region))
        mx = rbind(mx, region)
        
      }
    
    print((k/length(ls))*100) 
  }
  
  return(mx)}



##4 main sequences differences between paired sample 
#(Main-Difference-between-Paired, mdbp & lsmdbp) -----------------------------

mdbp<-function(ddf, ddg, nr){
  
  sVarN<-c()
  N1<-c()
  N2<-c()
  
  nrv=c(0, 500, 1000, 2000, 2500) #1 = 0; 2 = 500; 3 = 1000; 4 = 2000; 5 = 2500
  
  dddf<-read.csv(ddf)
  dddg<-read.csv(ddg)
  
  lth = length(dddf)
  
  for(i in 3: lth){
    
    if ( ( (dddf[1,i] | dddf[2,i] | dddf[3,i] | dddf[4,i] | dddf[5,i]) != 0 ) 
      & ( (dddg[1,i] | dddg[2,i] | dddg[3,i] | dddg[4,i] | dddg[5,i]) != 0 ) ) {
    
      tcc = sum(dddf[1,i] + dddf[2,i] + dddf[3,i] + dddf[4,i])  #total counts
      tccg = sum(dddg[1,i] + dddg[2,i] + dddg[3,i] + dddg[4,i])  #total counts
      
      mm = max(dddf[1,i], dddf[2,i], dddf[3,i], dddf[4,i]) #maximum counts (count of consensus seq)
      mmg = max(dddg[1,i], dddg[2,i], dddg[3,i], dddg[4,i]) #maximum counts (count of consensus seq)
      
      nmax = which.max(c(dddf[1,i], dddf[2,i], dddf[3,i], dddf[4,i]))
      nmaxg = which.max(c(dddg[1,i], dddg[2,i], dddg[3,i], dddg[4,i]))
      
      
      if ((nmax != nmaxg) & ( min(tcc, tccg) >= nrv[nr])){
        
        aamatrix=c("A","T","C","G")
        
        ii = i-2  

        sVarN[length(sVarN)+1] = (ii)
        N1[length(N1)+1] = aamatrix[nmax]
        N2[length(N2)+1] = aamatrix[nmaxg]
        
        
        }
    }
    
  }
  
  rb<-data.frame(sVarN, N1, N2)
  return(rb)
} #ddf = fist file; ddg = second file


#4.2 listed

lsmdbp<-function(ls, lss, nr){  #fist list and second list
  
  for (i in 1:length(ls)){
    
    assign(paste0("pdf",".",i), mdbp(ls[i], lss[i], nr)) #dataframe of sample
    print((i/length(ls)*100))
  }
  
  plist<-mget(paste0("pdf.",1:length(ls)))
  lldvar<-do.call("rbind", plist)
  
  
  return(lldvar)
  
}


##5 extract consens sequence ----------------------------------------

conSUM<-function(ls, nr){
  library(seqinr)
  mx=c()
  nrv=c(0, 500, 1000, 2000, 2500) #1 = 0; 2 = 500; 3 = 1000; 4 = 2000; 5 = 2500
  
  for(i in 1: length(ls)){
   
    df <- read.csv(ls[i])
    lth = length(df)
    
    cons<-c()
    
    for(k in 3: lth){
      
      tcc = sum(df[1,k] + df[2,k] + df[3,k] + df[4,k])  #total counts
      
      if (((df[1,k] | df[2,k] | df[3,k] | df[4,k] | df[5,k]) != 0 ) & (tcc >= nrv[nr])){
          
          aamatrix=c("A","T","C","G")
          nmax = which.max(c(df[1,k], df[2,k], df[3,k], df[4,k]))
          
          cons[length(cons) +1 ] = aamatrix[nmax]  }else{
      
          cons[length(cons) +1 ] = "-"}
        
    }
    
    assign(paste0("cdf.",i), cons) 
    mx=cbind(mx, get(paste0("cdf.",i)))
}
   
   mx=split(mx, rep(1:ncol(mx), each = nrow(mx)))
   
   write.fasta(mx, names=ls, file.out="mx.fasta",
               nbchar = max(sapply(a, function(x) length(x))))
  
   return(mx)
}

##6 extract 2nd last variant (%) -------------------------------------- 


sndvariantp<-function(ls, nr){
  
  pVarN<-c()

  nrv=c(0, 500, 1000, 2000, 2500) #1 = 0; 2 = 500; 3 = 1000; 4 = 2000; 5 = 2500
  
  for(k in 1: length(ls)){ #file
    
    dddf<-read.csv(ls[k])
    lth = length(dddf)
    
    for(i in 3: lth){
      
      tcc = sum(dddf[1,i] + dddf[2,i] + dddf[3,i] + dddf[4,i])  #total counts
      
      if (((dddf[1,i] | dddf[2,i] | dddf[3,i] | dddf[4,i] | dddf[5,i]) != 0 ) & 
          (tcc >= nrv[nr])) {
        
        a=c(dddf[1,i], dddf[2,i], dddf[3,i],dddf[4,i])
        
          if((sort(a, decreasing = TRUE)[2]/tcc) != 0){
            
            pVarN[length(pVarN) +1 ] = (sort(a, decreasing = TRUE)[2]/tcc)
            
          } 
        
      
        }
      }
     
     print((k/length(ls))*100)
    }
   return(pVarN) 
    
  }

##7 for error rate -------------------------------------- 

#1
ff1c = function(v){
  
  mx = which.max(v)
  v = v[-mx]  
  return(v)  
}

ff1t = function(v){
  
  fft1 = sum(v)
  
  mx = which.max(v)
  t = rep(fft1, length(v[-mx]))
  return(t)
}

#0
ff0c = function(v){
  
  mx = which.max(v)
  m0 = which(v == 0)
  
  m = c(mx, m0)
  v = v[-m]  
  return(v)  
}

ff0t = function(v){
  
  fft0 = sum(v)
  
  mx = which.max(v)
  m0 = which(v == 0)
  
  m = c(mx, m0)
  t = rep(fft0, length(v[-m]))
  return(t)
}

##8 for .bam file -------------------------------------- 

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


#8.1
bamseqprepID <- function(bam){
  
  library(pasillaBamSubset)
  library(stringr)
  library(Rsamtools)
  
  r.bam <- scanBam(BamFile(bam))
  
  rmcigar = c("N", "P", "=", "X")           # discard 4
  rmcigarID = c("I", "D")                   # discard 2
  
  # keep: S, H, M          
  lcigar = strsplit(r.bam[[1]]$cigar, split ="")
  rml1 = which(unlist(lapply(lcigar, function(x){any(rmcigar %in% x)})) == TRUE)
  rml2 = which(unlist(lapply(lcigar, function(x){any(rmcigarID %in% x)})) == TRUE)
  rml3 = which(r.bam[[1]]$pos > 3000)
  rml = unique(c(rml1, rml2, rml3))
  
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




