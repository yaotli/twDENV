###########################################################
##### This file is for global comparison ##################


# read-in
library(seqinr)
fasaa.g<-read.fasta(file.choose(), seqtype = "AA")
fasmx.g<-matrix("NA", length(fasaa.g), length(getSequence(fasaa.g[[1]])))

# extract seq from fasta
for( i in 1:length(fasaa.g)){
  
  fasmx.g[i,] =  getSequence(fasaa.g[[i]])
  
}

# find E region
fasmx.g.E<-fasmx.g[,281:775]

# polymorphism in each position

library(dplyr)

g.E.poly<-c() 
po<-c()

for(i in 1: dim(fasmx.g.E)[2]){
  
  aa=
  data.frame(table(fasmx.g.E[,i])) %>%
  filter(Var1 != "-") %>%
  filter(Var1 != "X") %>%
  filter(Freq != 1) %>%
  select(Var1, Freq)
  
  n = nrow(aa)
  
  g.E.poly[length(g.E.poly) +1 ] = n
  po[length(po) +1 ] = i
  
}


library(ggplot2)

gp <- data.frame(g.E.poly,as.numeric(po))
ggplot(gp,aes(x=po, y=g.E.poly)) + geom_line()

