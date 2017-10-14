# Epi and variants 

# result from lofreq: var_pi_vcf
# result from Poisson distribution: var_pi_PT
#
#

   lsorigin <- list.files("~/twDENV/allSUM/")
  lsarrange <- c(lsorigin[213:221], lsorigin[1:212])
combinedcsv <- do.call("rbind", lapply(  paste0("~/twDENV/allSUM/", lsarrange ) , read.csv))      

# lsvariantPT in Function.R

  # set proposed error for each position
  # length = 3000

  adjER = rep(0.001, 3000)
  lsarrange_clean <- gsub(pattern = "SUM", "", gsub(pattern = ".csv", "", lsarrange))
  
  
  # min coverage = 1000 (3)
  
  PT_e3_var_211 <- lsvariantPT(combinedcsv, 0.001, 3, adjER, lsarrange_clean)
  
  # add pi
  
  PT_e3_var_211_pi <- lapply(PT_e3_var_211, function(x){
    
    y1 <- x[,4]*(1 - x[,4])*2
    x[, 5] = y1
    colnames(x)[5] = "pi"
    return(x)
    
  })
    
    
# extract number of variants
  
  var_pi_PT <- t(sapply(PT_e3_var_211_pi, function(x){
    
    y.1 <- length( which( x[,1] %in% seq(937, 2421) == TRUE ) )
    
    y.2 <- sum(x$pi[ which( x[,1] %in% seq(937, 2421) == TRUE ) ]) / (2421 - 937 + 1)
      
    return(y = c(y.1, y.2))
    
  }))
  colnames(var_pi_PT) <- c("n", "pi")
  
  
# read-in epi date  
  
library(stringr)
    
  epi_demo0 <- read.csv("./Sources/epi_data_csv.csv",stringsAsFactors = FALSE)
  colnames(epi_demo0) = c("ID", "C1", "DHF", "Sex", "Age", 
                          "DoSamp", "IsoD", "TpIll", "Snd", "R", 
                          "X", "Clade", "D", "Info")
  
  # remove previous testing result
  
     epi_demo0 <- epi_demo0[,-c(11, 13)]
  epi_demo0$ID <- gsub(pattern = "C1", replacement = "_C1", epi_demo0$ID)
  
  epi_demo0[,12][ which(epi_demo0[,12] == "no data") ] = "Missing"
  epi_demo0[,12][ which(epi_demo0[,12] == "conta") ] = "Dirty"

  ordertem <- c()
  for(i in 1: dim(var_pi_PT)[1]){
    
    ordertem[ length(ordertem) + 1] = match(x = rownames(var_pi_PT)[i], epi_demo0$ID )
  }
  
  # order as variant data, and keep only 2001-2003 data
  epi_demo_ordered <- epi_demo0[ordertem, ][c(1:138), ]
  
  
  # combine previous data
  
  epi_demo_ordered[, c(13, 14)] = var_pi_PT[, c(1,2)][c(1:138), ]
  epi_demo_ordered[, c(15, 16)] = var_pi_vcf[, c(1,2)]
  
  colnames(epi_demo_ordered)[13:16] = c("PT_n", "PT_pi", "vcf_n", "vcf_pi")
  
  
# Figure: Epi_var
# 
# c1 = 0
# clade 
# IsoD  
  
  
library(ggplot2)  
library(dplyr)
  
epidemo_var_epi_clade_77 <- 
  
  epi_demo_ordered %>%
   filter(C1 == 0) %>%
  select(ID, Clade, IsoD, PT_n, vcf_n)
  
  
  # plotting 

  ggplot(data = epidemo_var_epi_clade_77, 
         
         aes(x=IsoD, y=PT_n, color = Clade)) + 
    
    geom_point(size=4, alpha = 0.9) +
    geom_point(color = "black", size=4.01, alpha = 1, shape = 1) +
    
    expand_limits(x=c(0.75:2.5)) +
    scale_x_continuous(breaks=seq(0.75,2.5,0.25), labels = seq(2001.75, 2003.5, by = 0.25)) + 
    
    ylab("No. Variant") + xlab("Year") + 
    # ggtitle("2001-2003") + 
    theme_bw()+
    theme(axis.title = element_text(size = 20), 
          axis.text.x = element_text(size = 13), 
          axis.text.y = element_text(size = 16), 
          plot.title = element_text(size = 25),
          legend.position = c(0.1, 0.845) ) +
    scale_color_manual(values = c("#619CFF", "#00BA38", "#F8766D")) 
  
  
# genetic diversity across clade and clinical phase
#
  
  
library(ggplot2)  
library(dplyr)

epidemo_var_TpIll_clade_77 =     
  
  epi_demo_ordered %>%
    filter(C1 == 0) %>%
    select(ID, Clade, TpIll, PT_n, PT_pi, vcf_n, vcf_pi)
  
  # define acute phase (0-3)

  epidemo_var_TpIll_clade_77[,8] = 0
  epidemo_var_TpIll_clade_77[, 8][which( epidemo_var_TpIll_clade_77$TpIll < 4)] = "0-3"
  epidemo_var_TpIll_clade_77[, 8][which( epidemo_var_TpIll_clade_77$TpIll >= 4 )] = "4+"

  colnames(epidemo_var_TpIll_clade_77)[8] = "Acute"


  # ggplot 
p=  
  ggplot(epidemo_var_TpIll_clade_77, aes(x = Clade, y = vcf_pi, fill = Acute) ) + 
    
    theme_bw()+
    scale_y_continuous(breaks = c(0, 0.0005, 0.001, 0.0015), 
                       labels = c(0, 0.05, 0.1, 0.15)) + 
    
    ylab("Genetic diversity (pi x 100)") + xlab("") + 
    theme(axis.title = element_text(size = 20), 
          axis.text.x = element_text(size = 15), 
          axis.text.y = element_text(size = 10), 
          plot.title = element_text(size = 25),
          legend.position = c(0.1, 0.845),
          legend.title = element_blank())  +
    geom_dotplot(binaxis = 'y', stackdir = 'center', 
                 position=position_dodge(1), binwidth = 0.00003) 
  
  # add interhost pi
 p + geom_point(x = 1, y = 0.00031, shape = 5, size = 5) + 
   geom_point(x = 2, y = 0.00044, shape = 5, size = 5) + 
   geom_point(x = 3, y = 0.00094, shape = 5, size = 5) 
  
  
  
  

### population and vatiants ----------------
    
csv.dist_var <- read.csv( "~/twDENV/Sources/dist_var_analysis_ed.csv", stringsAsFactors = FALSE)    
colnames(csv.dist_var)[4] <- "Group"
csv.dist_var$Group <- as.character( csv.dist_var$Group )

Period <- rep( "3", dim(csv.dist_var)[1] )
Period[ which(csv.dist_var$IsoD < 1 ) ] = "1"
Period[ which( 1 <= csv.dist_var$IsoD & csv.dist_var$IsoD < 1.7 ) ] = "2"

csv.dist_var <- data.frame( csv.dist_var, Period )

# 1 
a=
ggplot( csv.dist_var, aes( x = Group, y = n_new) ) + theme_bw() + 
  geom_violin( ) + geom_jitter( position=position_jitter(0.05), alpha = 0.5)

# 2
b=
ggplot( csv.dist_var, aes( x = Group, y = Population_density) ) + theme_bw() + 
  geom_violin( ) + geom_jitter( position=position_jitter(0.05), alpha = 0.5)

# 3 
c=
ggplot( csv.dist_var, aes( x = Period, y = n_new) ) + theme_bw() + 
  geom_violin( ) + geom_jitter( position=position_jitter(0.05), alpha = 0.5)




multiplot(a, b, c, cols = 1)  