getwd()
setwd("/Users/yaosmacbook/twDENV/SUM")

#138 list (01-03 viruses)
ls0102 = c(list.files(getwd())[213:221], list.files(getwd())[1:129])

#Site composition matrix : vdpon function --------
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
Mx3residue=data.frame(names138, mgroup, 
                      SiteMx1073[,3], SiteMx1073[,2], #C1073T
                      SiteMx1227[,2], SiteMx1227[,1], #T1227A
                      SiteMx2064[,2], SiteMx2064[,3]) #T2064C

attr(Mx3residue, "names") = c("Strain","Group", 
                           "pC1073", "p1073T",
                           "pT1227", "p1227A", 
                           "pT2064", "p2064C")

#to eliminate C1 sample or direct extract favored samples -------
nC1<-c(2,4,6,7,9,11,13,15,18,20,22,24,26,28,30,32,37,39,41,43,45,47,49,50,52,54,
       56, 58,60,62,64,66,68,70,72,74,76,79,82,84,86,88,90,93,95,97,99,103,105,
       106,112,113,116,119,121,123,126,128,133,136,138)

ns32<-c(1,5,7,8,16,17,2,3,22,28,30,38,45,23,24,27,34,44,46,48,50,53,54,57,60,
        61,66,70,71,74,75,77)

Mx3residuedC1<-Mx3residue[-nC1,] 
Mx3residue32s<-Mx3residuedC1[ns32,] #dim(Mx3residue32s)
Mx3residue32s[,9]=seq(1:32)
colnames(Mx3residue32s)[9] = "NO"


#generate the figures ------------------------------------------------

#color adjustment
coll=c("#56B4E9","#009E73","#CC79A7")

library(ggplot2)

#C1073T

f1073T = ggplot(data=Mx3residue32s, aes(NO, p1073T)) + 
                geom_line(size=1.5) + geom_point(size=3) + 
                
                theme_bw() +
                theme(
                      axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, 
                                          color=coll[Mx3residue32s$Group], size=15, face="bold"),
                      axis.title=element_text(size=16,face="bold"), 
                      axis.text=element_text(size=15)) + 
  
                scale_x_continuous(breaks=seq(1,32,1), 
                                   labels=Mx3residue32s[,1]) + xlab("") + ylab("Proportion")

ff1073T = f1073T +  geom_line(data=Mx3residue32s, aes(NO, pC1073), color="#bdbdbd") +
                    geom_point(data=Mx3residue32s, aes(NO, pC1073), color="#bdbdbd",size=2)

fff1073T = ff1073T + 
                    annotate("text", x = 30, y = 0.5, label = "1073 T", size=8, fontface="bold") + 
                    geom_vline(xintercept=8.5, linetype="dashed", color="#333333") #2001 -- 2002


#T1227A

f1227A = ggplot(data=Mx3residue32s, aes(NO, p1227A)) + 
                geom_line(size=1.5) + geom_point(size=3) + 
  
                theme_bw() +
                theme(
                      axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, 
                                  color=coll[Mx3residue32s$Group], size=15, face="bold"),
                      axis.title=element_text(size=16,face="bold"), 
                      axis.text=element_text(size=15)) + 
  
                scale_x_continuous(breaks=seq(1,32,1), labels=Mx3residue32s[,1]) + 
                xlab("") + ylab("Proportion")

ff1227A = f1227A +  geom_line(data=Mx3residue32s, aes(NO, pT1227), color="#bdbdbd") +
                    geom_point(data=Mx3residue32s, aes(NO, pT1227), color="#bdbdbd",size=2)

fff1227A = ff1227A + 
                  annotate("text", x = 30, y = 0.5, label = "1227 A", size=8, fontface="bold") + 
                  geom_vline(xintercept=8.5, linetype="dashed", color="#333333") #2001 -- 2002

#T2064C

f2064C = ggplot(data=Mx3residue32s, aes(NO, p2064C)) + 
                geom_line(size=1.5) + geom_point(size=3) + 
  
                theme_bw() +
                theme(
                      axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, 
                                  color=coll[Mx3residue32s$Group], size=15, face="bold"),
                      axis.title=element_text(size=16,face="bold"), 
                      axis.text=element_text(size=15)) +  
  
                      scale_x_continuous(breaks=seq(1,32,1), labels=Mx3residue32s[,1]) + 
                      xlab("") + ylab("Proportion")

ff2064C = f2064C +  geom_line(data=Mx3residue32s, aes(NO, pT2064), color="#bdbdbd") +
                     geom_point(data=Mx3residue32s, aes(NO, pT2064), color="#bdbdbd",size=2)

fff2064C = ff2064C + 
                   annotate("text", x = 30, y = 0.5, label = "2064 C", size=8, fontface="bold") + 
                   geom_vline(xintercept=8.5, linetype="dashed", color="#333333") #2001 -- 2002



#a few modification ------------------------------------------------

#mulitplot function from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
multiplot(fff1073T,fff1227A,fff2064C)




