# readin tree

library(ape)
library(ggtree)
library(stringr)

#   E: 129 tips
# ORF: 44

  E_tree <- read.tree("~/twDENV/Ongoing/E.nwk")
ORF_tree <- read.tree("~/twDENV/Ongoing/ORF.nwk")

  E_treedata <- fortify(E_tree) 
ORF_treedata <- fortify(ORF_tree)

# extract distacne to root 

dis_E <- dist.nodes(E_tree)[133, 1:129]
dis_ORF <- dist.nodes(ORF_tree)[45, 44]

# check the node 
# ggtree(E_tree) + geom_text(aes(label=node), size = 2 ) + geom_tiplab(size = 1)
# adjust the size when export as pdf

# manual enter clade 
  E_clade <- c(rep("2", 69), rep("Ib", 9), rep("Ia", 2), "F",
             rep("Ia", 3), "Ib", rep("Ia", 8), rep("Ib",2),
             rep("Ia", 3), rep("Ib",4), rep("Ia", 7), rep("F", 20))

ORF_clade <- c(rep("2", 28), rep("Ib", 8), rep("Ia", 7), "F")




E_disdata <- data.frame(name = 
                          str_match(
                            gsub(pattern = "'", replacement = "" , 
                                 E_treedata$label[1:129]), "([A-Za-z0-9-/]+)_")[, 2], 
                        
                        )





