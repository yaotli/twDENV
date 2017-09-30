require(ggtree)
require(seqinr)
source("./Function.R")

# files
trefile_E    <- "./Sources/E_129.tre"
trefile_subE <- "./Sources/sub_E.tre"
trefile_ORF  <- "./Sources/ORF_44.tre"
fas_aa       <- "./Sources/ORF_aa.fas"


# readin
E_mcc   <- read.beast( trefile_subE )
orf_mcc <- read.beast( trefile_ORF )
orf_aa  <- read.fasta( fas_aa )


## E tree ---------------

# define groups
tre.d <- fortify( E_mcc )
Ia.E  <- c( getDes(143), 143)
Ib.E  <- c( getDes(166), getDes(174), 166, 174)
II.E  <- c( getDes(182), 182)


f.orf <- 
  ggtree( E_mcc, right = TRUE, mrsd = "2009-01-03", size = 1) + theme_tree2() +
  scale_x_continuous( breaks = seq(1990.5, 2009.5, by = 2), 
                      labels = seq(1990, 2008, by = 2) ) + 
  theme( axis.text.x = element_text( size = 16 ), 
         axis.ticks  = element_blank())

rectdf <- data.frame( xstart = seq( 1990, 2008, 2), 
                      xend   = seq( 1991, 2009, 2))

# add timebar
f1.orf <- 
  f.orf + 
  geom_rect(data = rectdf, 
            aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf), 
            fill = "gray", alpha = 0.2, inherit.aes=FALSE) +
  scale_y_continuous( expand = c(0,2) ) 

# color group
tredata.E                              <- tre.d
tredata.E[, ncol(tredata.E) + 1 ]      <- "black"
colnames(tredata.E)[ ncol(tredata.E) ] <- "colorr"

tredata.E$colorr[ Ia.E ] <- "#1f77b4"
tredata.E$colorr[ Ib.E ] <- "#2ca02c"
tredata.E$colorr[ II.E ] <- "#9467bd"

f1.orf %<+% tredata.E + aes( color = I(colorr) )









## ORF tree ---------------
f.orf  <- ggtree( orf_mcc, right = TRUE, mrsd = "2002-12-30") + theme_tree2() 
rectdf <- data.frame( xstart = seq(1999, 2001, 2), xend = seq(2000, 2002, 2))

f.orf.rec <-  
  f.orf +  
  geom_rect(data = rectdf, 
          aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf), 
          fill = "gray", alpha = 0.2, inherit.aes=FALSE) +
  cale_y_continuous( expand = c(0,0.3) ) 

# prep heatmap

orf_max <- as.data.frame( do.call( rbind, getSequence(orf_aa) ), stringsAsFactors = FALSE)
colnames( orf_max ) = c( "46", "271", "357" )
rownames( orf_max ) = attributes( orf_aa )$names

f.orf.rec  %>%
gheatmap( orf_max, width = 0.1, colnames = FALSE) %>%
scale_x_ggtree( breaks = seq(1999.5, 2003.5, by = 1), 
                labels = seq(1999, 2002, by = 1) )



