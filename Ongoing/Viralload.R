qpcr <- read.csv("~/Desktop/qPCR.csv")

         vl_qpcr <- (qpcr[,2] + qpcr[,3]) / 2
      onset_qpcr <- qpcr[,10]
      clade_qpcr <- as.character(qpcr[,11])
   repeated_qpcr <- qpcr[,9]
   
   
   all_qpcr <- data.frame(as.character(qpcr[,1]), 
                          vl_qpcr, 
                          onset_qpcr,
                          clade_qpcr, 
                          repeated_qpcr)
   
   colnames(all_qpcr) = c("ID", "VL", "onset", "clade", "repeated") 
   
   phase_qpcr = c()
   for (i in 1: length(all_qpcr$onset)){
     
     
     if (all_qpcr$onset[i] > 10){
       
       phase_qpcr[i] = "NA"
       
     }
     
     if (all_qpcr$onset[i] <= 3){
       
       phase_qpcr[i] = "0-3"
       
      }
     
     if ( 3 < all_qpcr$onset[i] & all_qpcr$onset[i] <= 10  ){
       
       phase_qpcr[i] = "4-10"
       
     }
     
   }
   
   
   all_qpcr = cbind(all_qpcr, phase_qpcr)
   colnames(all_qpcr)[6] = "phase"
   
   
library("dplyr")
library(ggplot2)
   
  db1<- all_qpcr %>%
     filter(phase != "NA" )  %>%
     select(VL, phase, clade)
   
  # VL vs Clade   
  ggplot(db1, aes( x = clade, y = log10(VL), fill = phase)) + 
    geom_boxplot(alpha = 0.6) + 
    scale_x_discrete(name = "Clade") + 
    theme_bw() + 
    theme(
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          axis.title = element_text(face="bold"),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          axis.text.x = element_text(size = 15), 
          axis.text.y = element_text(size = 15), 
          legend.text = element_text(size = 16)) + 
    
    scale_fill_brewer(palette = "Accent") +
    guides(fill = guide_legend(title=NULL))

    
  # Variation vs VL
  
   
   
   
   
   
   
   
   
   
   