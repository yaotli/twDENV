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
   
   write.fasta(mx, names=c(1: length(ls)), file.out="mx.fasta",
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


