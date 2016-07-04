#setwd("/Users/yaosmacbook/twDENV/allSUM")

lsfull <- c(list.files(getwd())[213:221], list.files(getwd())[1:212])
ls0102 <- c(list.files(getwd())[213:221], list.files(getwd())[1:129])
ls05 <- c(list.files(getwd())[130:212])

#1 C0 sample vs C1 ----------------------------------------------------

# "ALL" C0 vs C1 + varies cutting value (0.05, 0.1, 0.2, 0.25)

lsfull.c <- strsplit(lsfull, split="-", fixed = T) #cut by "-"
cstatus <- sapply(lsfull.c, function(x) length(x) == 1 )
cstatus <- gsub("TRUE", replacement = "C0", cstatus)
cstatus <- gsub("FALSE", replacement = "C1", cstatus)         #all C status

lsfull.C0 <- lsfull[which(cstatus == "C0")]
lsfull.C1 <- lsfull[which(cstatus == "C1")]


rmSUM <-gsub("SUM", replacement = "", lsfull)
rmSUM.c <- strsplit(rmSUM, split="_", fixed = T) #cut by "-"

lsfull.NO <- sapply(rmSUM.c, function(x) head(x,1))           #all sample numbering


for (i in 1:4){
  
  a=c(0.05, 0.1, 0.2, 0.25)
  
  assign(paste0("lsv.C0.r",a[i]), lsvariant(lsfull.C0, a[i]))
  assign(paste0("lsv.C1.r",a[i]), lsvariant(lsfull.C1, a[i]))
}


n1 = c("lsv.C0.r0.05", "lsv.C0.r0.1", "lsv.C0.r0.2", "lsv.C0.r0.25")
n2 = c("lsv.C1.r0.05", "lsv.C1.r0.1", "lsv.C1.r0.2", "lsv.C1.r0.25")


nn1=c() for (i in 1:4){
  nn1[length(nn1)+1] = dim(get(n1[i]))[1]}


nn1n=nn1/151
nn2n=nn2/70

a<-cbind(nn1, nn2, nn1n, nn2n)
rownames(a)=c(0.05, 0.1, 0.2, 0.25)
colnames(a)=c("C0", "C1", "%", "%") #synomynous

#!consider the read number



# paired C0 vs C1 + varies cutting value (0.05, 0.1, 0.2, 0.25)

pairedS.No <- c()                                             #all number of paired C0
for(i in 2: length(lsfull.NO)){
  
  k=i -1
  if (strsplit(lsfull.NO[i], split="-", fixed = T)[[1]][1] == 
      strsplit(lsfull.NO[k], split="-", fixed = T)[[1]][1]){
    
    pairedS.No[length(pairedS.No)+1] = i-1
    
  }
}







#2 Variant sites by region


#2.1 Variants by sites
