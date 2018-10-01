# readin tree

library(ape)
library(ggtree)
library(stringr)
library(ggplot2)
library(dplyr)

#   E: 129 tips
# ORF: 44

  E_tree <- read.tree("~/twDENV/Sources//E.nwk")
ORF_tree <- read.tree("~/twDENV/Sources//ORF.nwk")

  E_treedata <- fortify(E_tree) 
ORF_treedata <- fortify(ORF_tree)

# extract distacne to root 

dis_E <- dist.nodes(E_tree)[133, 1:129]
dis_ORF <- dist.nodes(ORF_tree)[45, 1:44]

# check the node 
# ggtree(E_tree) + geom_text(aes(label=node), size = 2 ) + geom_tiplab(size = 1)
# adjust the size when export as pdf

# manual enter clade 
  E_clade <- c(rep("II", 69), rep("Ib", 9), rep("Ia", 2), "F",
             rep("Ia", 3), "Ib", rep("Ia", 8), rep("Ib",2),
             rep("Ia", 3), rep("Ib",4), rep("Ia", 7), rep("F", 20))

ORF_clade <- c(rep("II", 28), rep("Ib", 8), rep("Ia", 7), "F")


# make datafreame

E_disdata <- data.frame(name = 
                            str_match(
                              gsub(pattern = "'", replacement = "" , 
                                   E_treedata$label[1:129]), "([A-Za-z0-9-/]+)_")[, 2], 
                        clade = E_clade, 
                        distance = dis_E, 
                        
                        time = 
                            as.numeric(str_match(
                              gsub(pattern = "'", replacement = "" , 
                                   E_treedata$label[1:129]), "_([0-9.]+)")[,2]),
                        
                        adjYear = (as.numeric(str_match(
                          gsub(pattern = "'", replacement = "" , 
                               E_treedata$label[1:129]), "_([0-9.]+)")[,2]) - 2000)
                        
                        )


ORF_disdata <- data.frame(name = 
                          str_match(
                            gsub(pattern = "'", replacement = "" , 
                                 ORF_treedata$label[1:44]), "([A-Za-z0-9-/]+)_")[, 2],
                          clade = ORF_clade, 
                          distance = dis_ORF, 
                          
                          time = 
                          as.numeric(str_match(
                            gsub(pattern = "'", replacement = "" , 
                                 ORF_treedata$label[1:44]), "_([0-9.]+)")[,2]), 
                          
                          adjYear = (as.numeric(str_match(
                            gsub(pattern = "'", replacement = "" , 
                                 ORF_treedata$label[1:44]), "_([0-9.]+)")[,2]) - 2000)
                          
                        )


# figure 6 ------------

  E_disdata_sub <- E_disdata[which(E_disdata$clade != "F"), ]
ORF_disdata_sub <- ORF_disdata[which(ORF_disdata$clade != "F"), ]

# E 
ggplot(data = E_disdata_sub, aes(x = adjYear, y = distance, color = clade)) +
  theme_bw() +  
  
  scale_y_continuous(limits = c(0.004,0.008), labels = seq(4,8)) + 
  scale_x_continuous(breaks=seq(1,3.5,0.25), labels = seq(2001, 2003.5, by = 0.25) ) + 
  geom_smooth(method=lm, se=FALSE, size = 2) +
  
  geom_point(size=4, alpha = 0.9) +
  geom_point(color = "black", size=4.01, alpha = 1, shape = 1) +
  scale_color_manual(values = c("#619CFF", "#00BA38", "#F8766D")) +
  

  theme(axis.text = element_text(size=15), 
        axis.title = element_text(size=15,face="bold"), 
        legend.text = element_text(face="bold", size = 13), 
        legend.title = element_text(face="bold", size = 14),
        legend.position = c(1,0), 
        legend.justificatio = c(1,0) ) + 
  
  xlab("Year") + 
  ylab("Genetic distance (sub. / site x 10^3)")
  # ylab( expression(paste("Genetic distance(sub./ sitex", 10^{3}, ")")) )
    
    

# ORF

ggplot(data = ORF_disdata_sub, aes(x = adjYear, y = distance, color = clade)) +
  theme_bw() +  
  
  scale_y_continuous(limits = c(0.0025,0.0045), labels = seq(2.5, 4.5, by = 0.5)) + 
  
  # make the x-axis include 1.75
  expand_limits(x=c(1.75:3.25)) + 
  scale_x_continuous(breaks=seq(1,3.5,0.25), labels = seq(2001, 2003.5, by = 0.25)) + 
  geom_smooth(method = lm, se = FALSE, size = 2) +
  
  geom_point(size=4, alpha = 0.9) +
  geom_point(color = "black", size=4.01, alpha = 1, shape = 1) +
  scale_color_manual(values = c("#619CFF", "#00BA38", "#F8766D")) +
  
  theme(axis.text = element_text(size=15), 
        axis.title = element_text(size=15,face="bold"), 
        legend.text = element_text(face="bold", size = 13), 
        legend.title = element_text(face="bold", size = 14),
        legend.position = c(1,0), 
        legend.justificatio = c(1,0) ) + 
  
  xlab("Year") + 
  ylab("Genetic distance (sub. / site x 10^3)")


# slope ####

v_clade <- c()
v_gene <- c()
v_slope <- c()
v_sd <- c()
v_p <- c()
v_range1 <- c()
v_range2 <- c()

for (i in 1:10){
  
  if (i %in% c(1,2,3,4,5,6)){
    
     cladeg <- c("Ia", "Ia", "II", "II", "Ib", "Ib")
    disdate <- list(E_disdata_sub, ORF_disdata_sub, E_disdata_sub, ORF_disdata_sub, 
                    E_disdata_sub, ORF_disdata_sub)
    
    y = disdate[[i]]$distance[which(disdate[[i]]$clade == cladeg[i] )]
    x = disdate[[i]]$adjYear[which(disdate[[i]]$clade == cladeg[i] )]
    
    range = range( x )
      lm0 = lm( y ~ x )
    
     v_clade <- c(v_clade, cladeg[i])
     v_slope <- c(v_slope, lm0$coefficients[2])
        v_sd <- c(v_sd, summary(lm0)$coefficients[2,2])
         v_p <- c(v_p, summary(lm0)$coefficients[2,4])
    v_range1 <- c(v_range1, range[1])
    v_range2 <- c(v_range2, range[2])
    
  }
  
  if (i %in% c(7,8,9,10)){
    
     cladeg <- "Ib"
    disdate <- list(E_disdata_sub, E_disdata_sub,
                    ORF_disdata_sub, ORF_disdata_sub)
    
    k = i - 6
    
    y = disdate[[k]]$distance[which(disdate[[k]]$clade == cladeg )]
    x = disdate[[k]]$adjYear[which(disdate[[k]]$clade == cladeg )]
    
    if (i %in% c(7,9)){
      
      y = y[which(x < 2.25)]
      x = x[which(x < 2.25)]
      
      range = range( x )
        lm0 = lm( y ~ x )

    }else{
      
      y = y[which(x >= 2.25)]
      x = x[which(x >= 2.25)]
      
      range = range( x )
      lm0 = lm( y ~ x )
      
    }
    
    v_clade <- c(v_clade, cladeg)
    v_slope <- c(v_slope, lm0$coefficients[2])
    v_sd <- c(v_sd, summary(lm0)$coefficients[2,2])
    v_p <- c(v_p, summary(lm0)$coefficients[2,4])
    v_range1 <- c(v_range1, range[1])
    v_range2 <- c(v_range2, range[2])
    
  }
  
  }

slopedata <- data.frame(v_clade, v_slope, v_sd, v_p, v_range1, v_range2)  
slopedata <- data.frame(slopedata, gene = c("E", "ORF", "E", "ORF", "E", "ORF", 
                                            "E", "E", "ORF", "ORF"))

