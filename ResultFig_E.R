getwd()
setwd("/Users/yaosmacbook/twDENV/SUM")

#138 list
ls0102 = c(list.files(getwd())[213:221], list.files(getwd())[1:129])

#Site composition matrix --------
SiteMx1073=vdpon(1073,ls0102)
SiteMx1227=vdpon(1227,ls0102)
SiteMx2064=vdpon(2064,ls0102)

#assign groups of virus ---------
mgroup = c(rep("Ia", 35), rep("II", 103))
a=c(3,4,5,6,23,24,25,26,40,41,53,54,57,58,73,74,80,85,86)
for (i in 1:length(a)){
  
  mgroup[a[i]] = "Ib"
  
  
}

#strain name by outer csv
names138 <- read.csv(file.choose(), header = FALSE)
names138 <- as.character(names138[,1])

#to build the Data source matrix
mdf1227=data.frame(names138, mgroup, rmgroup, SiteMx[,1], SiteMx[,2], SiteMx[,3], SiteMx[,4])
attr(mdf1227, "names") = c("Strain", "Group", "rGroup", "pA", "pT", "pC", "pG")

#c1
ctt<-c(2,4,6,7,9,11,13,15,18,20,22,24,26,28,30,32,37,39,41,43,45,47,49,50,52,54,
       56, 58,60,62,64,66,68,70,72,74,76,79,82,84,86,88,90,93,95,97,99,103,105,106,112,113,116,
       119,121,123,126,128,133,136,138)

tc<-c(1,5,7,8,16,17,2,3,22,28,30,38,45,23,24,27,34,44,46,48,50,53,54,57,60,61,66,70,71,74,75,77)

class(mdf1227dT)

mdf1227dT<-mdf1227d[tc,]
mdf1227dT[,8]=seq(1:32)
colnames(mdf1227dT)[8] = "NO"


library(ggplot2)
v = ggplot(data=mdf1227dT, aes(NO, pA)) + geom_line(size=1.5) + geom_point(size=3) + theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, color=coll[mdf1227dT$Group], size=15, face="bold"), 
        axis.title=element_text(size=16,face="bold"), axis.text=element_text(size=15)) + 
  scale_x_continuous(breaks=seq(1,32,1), labels=mdf1227dT[,1]) + xlab("") + ylab("Proportion")
v1 = v+geom_line(data=mdf1227dT, aes(NO, pT), color="#bdbdbd") +
  geom_point(data=mdf1227dT, aes(NO, pT), color="#bdbdbd",size=2)

vv = v1 + annotate("text", x = 30, y = 0.5, label = "1227 A", size=8, fontface="bold") + 
  geom_vline(xintercept=8.5, linetype="dashed", color="#333333")



coll=c("#56B4E9","#009E73","#CC79A7")

