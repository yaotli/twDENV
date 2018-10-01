require(ggtree)
require(seqinr)
source("./Function.R")

# files
trefile_E    <- "./Sources/E_129.tre"
trefile_subE <- "./Sources/sub_E.tre"
trefile_ORF  <- "./Sources/ORF_44.tre"
fas_aa       <- "./Sources/ORF_aa.fas"
fas_nt_E     <- "./Sources/E_nt.fas"


# readin
E_mcc   <- read.beast( trefile_E )
orf_mcc <- read.beast( trefile_ORF )
orf_aa  <- read.fasta( fas_aa )
e_nt    <- read.fasta( fas_nt_E )


## E tree with nt (1227, 1073, 2064) ---------------

# define groups
tre.d <- fortify( E_mcc )
Ia.E  <- c( getDes(142), 142)
Ib.E  <- c( getDes(242), getDes(166), 242, 166)
II.E  <- c( getDes(174), 174)

e_max <- as.data.frame( do.call( rbind, getSequence(e_nt) ), stringsAsFactors = FALSE)
e_max <- e_max[ c( 137, 291, 1128 ) ]
colnames( e_max ) = c( "1073", "1227", "2064" )
rownames( e_max ) = attributes( e_nt )$names
e_max$`1227`[ which( e_max$`1227` == "w" ) ] = "t"


p1 = 
ggtree( E_mcc, right = TRUE, mrsd = "2009-01-01", size = 1) + theme_tree2() 
p1 %>% gheatmap(e_max, width = 0.1, color = "black") + 
  scale_fill_manual( breaks = c("a", "c", "g", "t"), 
                     values =  c("#D55E00", "#56B4E9", "#E69F00", "#009E73") )
  

f.e <- 
  ggtree( E_mcc, right = TRUE, mrsd = "2009-01-01", size = 1) + theme_tree2() +
  scale_x_continuous( breaks = seq(1998.5, 2009.5, by = 2), 
                      labels = seq(1998, 2009, by = 2), limit = c(1998, 2005)) +
  scale_y_continuous( limit = c(0, 117), expand = c(0,2) ) +
  theme( axis.text.x = element_text( size = 20 ), 
         axis.ticks  = element_blank())

rectdf <- data.frame( xstart = seq( 1990, 2008, 2), 
                      xend   = seq( 1991, 2009, 2))

# add timebar
f1.e <- 
  f.e + 
  geom_rect(data = rectdf, 
            aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf), 
            fill = "gray", alpha = 0.2, inherit.aes=FALSE) 

# color group
tredata.E                              <- tre.d
tredata.E[, ncol(tredata.E) + 1 ]      <- "black"
colnames(tredata.E)[ ncol(tredata.E) ] <- "colorr"

tredata.E$colorr[ Ia.E ] <- "#1f77b4"
tredata.E$colorr[ Ib.E ] <- "#2ca02c"
tredata.E$colorr[ II.E ] <- "#9467bd"

f1.e %<+% tredata.E + aes( color = I(colorr) ) + geom_tippoint()



## ORF tree  ---------------

# define groups
tre.d <- fortify( orf_mcc )
Ia.orf  <- c( getDes(82), 82)
Ib.orf  <- c( setdiff(getDes(47), c( getDes(55), 55)), 47)
II.orf  <- c( getDes(55), 55)

# color 
f.orf  <- ggtree( orf_mcc, right = TRUE, mrsd = "2002-12-22", size = 1) + theme_tree2() 

tredata.orf                                <- tre.d
tredata.orf[, ncol(tredata.orf) + 1 ]      <- "black"
colnames(tredata.orf)[ ncol(tredata.orf) ] <- "colorr"

tredata.orf$colorr[ Ia.orf ] <- "#1f77b4"
tredata.orf$colorr[ Ib.orf ] <- "#2ca02c"
tredata.orf$colorr[ II.orf ] <- "#9467bd"

f.orf <- 
f.orf %<+% tredata.orf + aes( color = I(colorr) ) + geom_tippoint()


# time bar
rectdf <- data.frame( xstart = seq(2000, 2002, 2), xend = seq(2001, 2003, 2))

f.orf.rec <-  
  f.orf %>%  scale_x_ggtree( breaks = seq(1999.5, 2002.5, by = 1), labels = seq(1999, 2002, by = 1) ) + 
  geom_rect(data = rectdf, 
          aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf), 
          fill = "gray", alpha = 0.2, inherit.aes=FALSE) +
  scale_y_continuous( expand = c(0,0.3) )  
  


# prep heatmap

orf_max <- as.data.frame( do.call( rbind, getSequence(orf_aa) ), stringsAsFactors = FALSE)
colnames( orf_max ) = c( "46", "271", "357" )
rownames( orf_max ) = attributes( orf_aa )$names

f.orf.rec  %>%
gheatmap( orf_max, width = 0.1, color = "black", colnames_position = "top", colnames_angle = 90, offset = 0.25, font.size = 1)+
#scale_x_ggtree( breaks = seq(1999.5, 2002.5, by = 1), 
#                labels = seq(1999, 2002, by = 1) ) +
theme_tree2( legend.position = c(0,0), legend.justification = c(0,0), 
             axis.text.x = element_text( size = 20 )) + 
  theme( axis.ticks  = element_blank() )

  


