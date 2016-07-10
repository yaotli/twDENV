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

lsfull.C0 <- lsfull[which(cstatus == "C0")] #n = 151
lsfull.C1 <- lsfull[which(cstatus == "C1")] #n = 70


rmSUM <-gsub("SUM", replacement = "", lsfull)
rmSUM.c <- strsplit(rmSUM, split="_", fixed = T) #cut by "-"

lsfull.NO <- sapply(rmSUM.c, function(x) head(x,1))           #all sample numbering


### ----- pooled
#syn.r1000.table0605

for (i in 1:4){
  
  a=c(0.05, 0.1, 0.2, 0.25)
  
  assign(paste0("outlsv.C0.1000r",a[i]), lsvariantNS(lsfull.C0, a[i], 1))
  assign(paste0("outlsv.C1.1000r",a[i]), lsvariantNS(lsfull.C1, a[i], 1))
}


n1 = c("outlsv.C0.1000r0.05", "outlsv.C0.1000r0.1", "outlsv.C0.1000r0.2", "outlsv.C0.1000r0.25")
n2 = c("outlsv.C1.1000r0.05", "outlsv.C1.1000r0.1", "outlsv.C1.1000r0.2", "outlsv.C1.1000r0.25")

nn1=c()
nn2=c() 
aa1=c()
aa2=c()
for (i in 1:4){
  nn1[length(nn1)+1] = dim(get(n1[i]))[1]
  nn2[length(nn2)+1] = dim(get(n2[i]))[1]
  aa1[length(aa1)+1] = length(which(get(n1[i]) == "TRUE"))
  aa2[length(aa2)+1] = length(which(get(n2[i]) == "TRUE"))}


nn1n=nn1/151
nn2n=nn2/70
aa1n=aa1/151
aa2n=aa2/70


a<-cbind(nn1, nn2, nn1n, nn2n, aa1, aa2, aa1n, aa2n)  #r1000.table0606am1 = a
rownames(a)=c(0.05, 0.1, 0.2, 0.25)
colnames(a)=c("C0", "C1", "%", "%", "aC0", "aC1", "a%", "a%") 



# paired C0 vs C1 + varies cutting value (0.05, 0.1, 0.2, 0.25) ------------------------

pairedS.No <- c()                                             #all number of paired C0 #n = 65
for(i in 2: length(lsfull.NO)){
  
  k=i -1
  if (strsplit(lsfull.NO[i], split="-", fixed = T)[[1]][1] == 
      strsplit(lsfull.NO[k], split="-", fixed = T)[[1]][1]){
    
    pairedS.No[length(pairedS.No)+1] = i-1
  }
}


### test for consensus sequence
#pairedS.Nom <- pairedS.No[1:64]

b=pairedS.No + 1
a<-lsmdbp(lsfull[pairedS.Nom], lsfull[pairedS.Nom + 1], 2) # table.mdbp.r500.0708am02 = a
                                                           # table.mdbp.r500.0708pm03 = a

### ---- deal with table of consensus change between pairs # table.mdbp.r500.0711am12
  
  a<-c(17, 18, 32, 38, 40, 40, 57, 57)
  aa<-c(1,4,3,4,4,1,3,1)
  aaa<-c(4,1,2,1,2,4,2,4)
  
  b=c()
  bb=c()
  bbb=c()
for(i in 1: length(a)){
  
   ii = i+1
   b[length(b)+1] = lsfull.NO[ pairedS.No [a[i]]]
   bb[length(bb)+1] = pon(lsfull[pairedS.No[a[i]]], aa[i], table.mdbp.r500.0708pm03[,1][i])
   bbb[length(bbb)+1] = pon(lsfull[pairedS.No[a[i]] + 1], aaa[i], table.mdbp.r500.0708pm03[,1][i])}
  
    t<-data.frame(b,table.mdbp.r500.0708pm03[,1], as.character(table.mdbp.r500.0708pm03[,2]), as.numeric(bb), 
                  as.character(table.mdbp.r500.0708pm03[,3]), as.numeric(bbb), c("N", rep("S", 6), "N"))
    colnames(t)=c("ID", "Position", "C0", "%", "C1", "%", "Type") # table.mdbp.r500.0711am12 = t

#### use cdvariant ------------------

for (k in 1:4){
  
  v=c(0.05, 0.1, 0.2, 0.25)

  mx = c(0,0,0,0,0,0,0,0)
  for (i in 1:length(pairedS.No)){
    
    ii = pairedS.No[i]
    iii = ii + 1
    
    a = cdvariant(lsfull[ii], v[k], 2) #2 = 500
    b = cdvariant(lsfull[iii], v[k], 2)
    c = c (a[1], b[1], a[2], b[2], a[3], b[3], a[4], b[4])
    
    mx = rbind(mx, c)
    
    print(i/length(pairedS.No)*100)
  }
  
  mx<-mx[-1,]
  rownames(mx) <- c(1:length(pairedS.No))
  assign(paste0("outlsv.paired.500r.", v[k]), mx)

}

#outlsv.paired.1000r.0.05
#outlsv.paired.1000r.0.1
#outlsv.paired.1000r.0.2
#outlsv.paired.1000r.0.25



a = melt(cbind(outlsv.paired.500r.0.05[,1], outlsv.paired.500r.0.05[,2], 
               outlsv.paired.500r.0.1[,1], outlsv.paired.500r.0.1[,2],
               outlsv.paired.500r.0.2[,1], outlsv.paired.500r.0.2[,2],
               outlsv.paired.500r.0.25[,1], outlsv.paired.500r.0.25[,2]))

#cp 
a[,4] = c(rep(0.05, 130), rep(0.1, 130), rep(0.2, 130), rep(0.25, 130))

#type
a[,5] = rep(c(rep("C0", 65), rep("C1", "65")), 4)
names(a) = c("Var1", "Var2", "count", "Cutpoint", "Type")
aa<-a[-(which(a$Var1== 65)),] #get rid of #65 pair


library(ggplot2)

ggplot(data=a, aes(x=Var2, y=value, group=Var1, color=Var1)) + geom_line()


# Pirate plot ----------

install.packages("devtools")
library("devtools")
install_github("ndphillips/yarrr")
library(yarrr)

#500 r
p=pirateplot(formula = count ~ Type + Cutpoint, data = a, theme.o = 3, 
           main="Comparison of paired samples", inf = "ci",
           ylab="No. of position with variants", gl.col = gray(.8),
           pal="basel")

box(which = "p")

## conventional one ------------------


for (i in 1:4){
  
  a=c(0.05, 0.1, 0.2, 0.25)
  
  assign(paste0("outlsv.pC0.1000r",a[i]), lsvariantNS(lsfull[pairedS.No], a[i], 3))
  assign(paste0("outlsv.pC1.1000r",a[i]), lsvariantNS(lsfull[pairedS.No + 1], a[i], 3))
}


n1 = c("outlsv.pC0.1000r0.05", "outlsv.pC0.1000r0.1", "outlsv.pC0.1000r0.2", "outlsv.pC0.1000r0.25")
n2 = c("outlsv.pC1.1000r0.05", "outlsv.pC1.1000r0.1", "outlsv.pC1.1000r0.2", "outlsv.pC1.1000r0.25")

nn1=c()
nn2=c() 
aa1=c()
aa2=c()
for (i in 1:4){
  nn1[length(nn1)+1] = dim(get(n1[i]))[1]
  nn2[length(nn2)+1] = dim(get(n2[i]))[1]
  aa1[length(aa1)+1] = length(which(get(n1[i]) == "TRUE"))
  aa2[length(aa2)+1] = length(which(get(n2[i]) == "TRUE"))}


nn1n=nn1/65
nn2n=nn2/65
aa1n=aa1/65
aa2n=aa2/65


a<-cbind(nn1, nn2, nn1n, nn2n, aa1, aa2, aa1n, aa2n)  #r1000.paired.table0607am12 = a
rownames(a)=c(0.05, 0.1, 0.2, 0.25)
colnames(a)=c("C0", "C1", "%", "%", "aC0", "aC1", "a%", "a%") 


#2 Repeated sample





#4 DF/ DHF + unset date

