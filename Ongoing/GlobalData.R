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


##################################################
##### NS1 

# tip ID cleaning 

library(seqinr) 

file = read.fasta(file.choose())

  seq.name0 = attributes(file)$names
  seq0 = getSequence(file)
  
  duplicated.id = which(duplicated(seq.name0) == "TRUE")
  
  seq.name0 = gsub("\\(", "-", seq.name0   )
  seq.name0 = gsub("\\)", "-", seq.name0   )
  seq.name0 = gsub("\\[", "/", seq.name0   )
  seq.name0 = gsub("\\]", "/", seq.name0   )
  seq.name0 = gsub(" ", "_", seq.name0   )
  seq.name0 = gsub("\\'", "", seq.name0   )
  seq.name0 = gsub(">", "", seq.name0   )
  
# replace ":number..number"
  
  a = "\\:([0-9]+)\\.\\.([0-9]+)\\_"
  seq.name = gsub(a, "_", seq.name0 )
  seq.name = gsub("\\/", "-", seq.name ) 
  
# deal with year
  
  for(i in 1:26){
    seq.name[i] = paste0(seq.name[i],"0")} # write.csv(seq.name, "a.csv")

  write.fasta(seq0, file.out = "revised.fasta", names = seq.name)
  
# extract the time 
  
  install.packages("stringr")
  library(stringr)

  a = "\\_([0-9]+)"                   #   grepl(a, seq.name)
  dvtime <- as.numeric(str_match(seq.name, a)[,2])
  
  a = "\\_([^0-9]+)\\_"                   #   grepl(a, seq.name)
  # seq.name[which(grepl(a, seq.name) == FALSE)]
  dvnation <- str_match(seq.name, a)[,2]
  dvnation[which(grepl(a, seq.name) == FALSE)] = "NA"

###### Tree
  
  library(ggplot2)
  library(ggtree)
  
  library(ape)
  
  eTree <-read.tree(file.choose())
  ggtree(eTree)

  eTreedata <- fortify(eTree)   # tablized tree data
  
  eTree$tip.label[1:10]
  seq.name[1:10]

  # dvtipinfo = data.frame(seq.name, dvnation, dvtime)
  
# target labeling   
  
  dvnation.taiwan = c()
  for(i in 1:length(dvnation)){
    if(dvnation[i] == "Taiwan"){ dvnation.taiwan[i] = 1 } else { dvnation.taiwan[i] = 0  }}
  
# color 
  
  Tcolor = rep("Black", length(seq.name))
  Tcolor[which(dvnation.taiwan == 1)] = "Red"
  Tcolor[1125:1127] = "Green" # 2016
  Tcolor[c(19, 900)] = "Blue" # Others
  
# order code: E
  
  eOrder = c()
  for(i in 1:length(seq.name)){
    
    eOrder[i] = which( seq.name == eTree$tip.label[i] )  }
  
  a=ggtree(eTree) + geom_tippoint(color = Tcolor[eOrder])



# for ns1
  
  ns1Tree <-read.tree(file.choose())
  ggtree(ns1Tree)

  ns1Order = c()
  for(i in 1:length(seq.name)){
    
    ns1Order[i] = which( seq.name == ns1Tree$tip.label[i] )  }
  
  b=ggtree(ns1Tree) + geom_tippoint(color = Tcolor[ns1Order])
  
  multiplot(a,b)
  
  
  


