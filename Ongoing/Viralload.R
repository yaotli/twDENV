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
   
  # VL vs Clade   ####
 c1 <- ggplot(db1, aes( x = phase, y = log10(VL), fill = clade)) + 
    geom_boxplot(alpha = 0.6, size = 1.5) + 
    scale_x_discrete(name = "Onset") + 
    theme_bw() + 
    theme(
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          axis.title = element_text(face="bold"),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 25),
          axis.text.x = element_text(size = 15), 
          axis.text.y = element_text(size = 15), 
          legend.text = element_text(size = 15)) + 
    
    scale_fill_manual(values = c("blue", "green", "purple")) +
    
    guides(fill = guide_legend(title=NULL))

    
  # Variation vs VL ####
  # c0.dhfonset.Dcon.Oct.0102
  
  db2<- all_qpcr %>%
    filter(phase != "NA" )  %>%
    select(ID, VL, clade)
   
  
  matchid <- c()
   
  m = gsub("DEN", "", as.character(c0.dhfonset.Dcon.Oct.0102[,1])) 
  for(k in 1: length(db2[,1])){
    
    if ( length(grep(as.character(db2[,1])[k], m)) == 0 ){
      
      matchid[k] = 10000
      
    }else{
      
      matchid[k] <- grep(as.character(db2[,1])[k], m)
      
    }
    
  }
    
  db2[,4] <- c0.dhfonset.Dcon.Oct.0102$out.n.c0[matchid]
  db2_a = db2[-which(is.na(db2[,4])),]
  
  # dot plot 
  ggplot(db2_a, 
         aes(x = log10(VL), y = V4, color = clade)) + 
    
    theme_bw() +
    
    # with ring-like effect
    geom_point(color = "black", size = 3.2, alpha = 0.8, shape = 1) + 
    geom_point(size = 3, alpha = 0.8) + 
    
    ylab("Variation") +  
    
    theme(
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      axis.title = element_text(face="bold"),
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      axis.text.x = element_text(size = 15), 
      axis.text.y = element_text(size = 15), 
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 20))
  
  # Variation vs VL stratified by clade ####
  
  p <- ggplot(db2_a, aes(x = V4, 
                    y = log10(VL), 
                    color = clade, group = clade)) + 
    
    
    geom_point(color = "black", size = 3.2, alpha = 0.8, shape = 1) + 
    geom_point(size = 3, alpha = 0.8) + 
  
    geom_smooth(method=lm, se = FALSE) +
    facet_wrap(~ clade, ncol = 1) +
    theme_bw() +
    xlab("Variation") +
    scale_x_continuous(breaks = seq(0, 60, by = 10), 
                       limits = c(0,60)) + 
    
    theme(
          axis.title = element_text(face="bold"),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          strip.background = element_rect(fill = "white", color = "white"),
          strip.text = element_text(size = 16, face = "bold"), 
          legend.position="none"
          )
  
   # lm
  
  lm1_db2_a <- lm(log10(db2_a$VL[which(db2_a$clade == "Ia")]) ~ 
                   db2_a$V4[which(db2_a$clade == "Ia")])

  lm2_db2_a <- lm(log10(db2_a$VL[which(db2_a$clade == "Ib")]) ~ 
                   db2_a$V4[which(db2_a$clade == "Ib")])
  
  lm3_db2_a = lm(log10(db2_a$VL[which(db2_a$clade == "II")]) ~ 
                   db2_a$V4[which(db2_a$clade == "II")])
  
  
  summary(lm1_db2_a)
  
  # add R-square
  
  Rlendend <- data.frame( x = rep(50, 3), 
                          y = rep(7, 3),
                          clade = c("Ia","Ib","II"),
                          labs = c('R^2 == 0.05', 'R^2 == 0.68', 'R^2 == 0.55' ))
  
  p + geom_text( data = Rlendend,
                 aes(x, y, 
                    label = labs, 
                    group = NULL),
                parse = TRUE,
                color = "Black", size = 5)
  
    
  # Stritification of viral load ####
  
  dbx<- all_qpcr %>%
    filter(phase != "NA" )  %>%
    select(ID, VL, clade, onset)
  
  dbx[,5] <- c0.dhfonset.Dcon.Oct.0102$out.n.c0[matchid]
  dbx_a <- dbx[-which(is.na(dbx[,5])),]
  
  # dot plot (onset vs VL)
  ggplot(dbx_a, 
         aes(x = onset, y = log10(VL), color = clade)) + 
    
    theme_bw() +
    
    # with ring-like effect
    geom_point(color = "black", size = 3.2, alpha = 0.8, shape = 1) + 
    geom_point(size = 3, alpha = 0.8) + 
    
    theme(
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      axis.title = element_text(face="bold"),
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      axis.text.x = element_text(size = 15), 
      axis.text.y = element_text(size = 15), 
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 20)) 
  
  # regressions
  
  lm1 = lm(log10(dbx_a$VL) ~ dbx_a$onset)
  lm2 = lm(log10(dbx_a$VL) ~ dbx_a$V5)
  lm3 = lm(dbx_a$V5 ~ dbx_a$onset)
  lm4 = lm(log10(dbx_a$VL)[ which(dbx_a$onset < 6) ] ~ dbx_a$V5[ which(dbx_a$onset < 6) ])
  
  plot(log10(dbx_a$VL)[ which(dbx_a$onset < 6) ] ~ dbx_a$V5[ which(dbx_a$onset < 6) ])
  
  summary(lm4)
  
  # ploty
  # https://plot.ly/~yaoli/4/ia-ib-ii/
  
   
  # Temporal pattern of VL ####
  
  db3 <- all_qpcr %>%
    filter(repeated != "NA" )  %>%
    filter(onset <= 10) %>%
    select(repeated, VL, onset, clade)
  
  # elimiate repeated 1715 since one time point  
  db3_b <- db3[-which(db3$repeated == 1715),]
  
  ggplot(db3_b, aes(x=onset, 
                  y=log10(VL), 
                  color=clade, group=repeated)) + 
    
    geom_line(size = 1.1) + 
    geom_point() +
    facet_wrap(~ repeated) + 
    scale_x_continuous(breaks = seq(0,10, by = 2)) + 
    scale_y_continuous(breaks = seq(2,8, by = 2)) + 
    theme_bw() + 
    xlab("Onset day") +
    theme(axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          strip.background = element_rect(fill = "white", color = "white"),
          strip.text = element_text(size = 11, face = "bold")) +
    
    scale_color_discrete(name = "Clade")
    
           
   
   
   
   
   
   
   