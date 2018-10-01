#setwd("~/twDENV/method_ngs/")

require(pasillaBamSubset)
require(Rsamtools)
require(stringr)
require(ggplot2)
require(ggpubr)


ls_bam  <- paste0( "./raw/bam/", list.files("./raw/bam/") )
ls_read <- 
  lapply( as.list(ls_bam),
          function(x)
          {
            r <- scanBam( BamFile( x ) )
            v <- as.vector( coverage( IRanges( r[[1]]$pos, width = r[[1]]$qwidth ) ) )
            if( length(v) > 3000 )
            {
              v <- matrix( v[801:3000], nrow = 1, byrow = TRUE)
            }else
            {
              v <- matrix( c( v[ c( 801: length(v) ) ], rep( 0, 3000 - length(v) ) ), nrow = 1, byrow = TRUE)  
            }
          return(v)
          } 
          )
 
mx_read  <- do.call("rbind", ls_read)
df_depth <- data.frame( pos = seq(801, 3000), 
                        med = apply(mx_read, 2, quantile)[3,],
                        q1  = apply(mx_read, 2, quantile)[2,], 
                        q3  = apply(mx_read, 2, quantile)[4,] )
d1 <- 
ggplot( data = df_depth ) + 
  geom_line( aes(x = pos, y = q1/1000), color = "gray", size = 1.5) + 
  geom_line( aes(x = pos, y = q3/1000), color = "gray", size = 1.5) + 
  geom_line( aes(x = pos, y = med/1000), size = 1.8) +
  
  #1f + 1r
  geom_segment( x = 884, xend = 901, y = 41, yend = 41, color = "blue", size = 2) + 
  geom_segment( x = 1357, xend = 1376, y = 41, yend = 41, color = "red", size = 2) + 
  geom_segment( x = 884, xend = 1376, y = 41, yend = 41, color = "black", size = 0.25, linetype = "dashed") +
  #2f + 1r
  geom_segment( x = 1310, xend = 1327, y = 40, yend = 40, color = "blue", size = 2) + 
  geom_segment( x = 1752, xend = 1770, y = 40, yend = 40, color = "red", size = 2) + 
  geom_segment( x = 1310, xend = 1770, y = 40, yend = 40, color = "black", size = 0.25, linetype = "dashed") +
  
  #3f + 1r
  geom_segment( x = 1616, xend = 1635, y = 41, yend = 41, color = "blue", size = 2) + 
  geom_segment( x = 2035, xend = 2053, y = 41, yend = 41, color = "red", size = 2) + 
  geom_segment( x = 1616, xend = 2053, y = 41, yend = 41, color = "black", size = 0.25, linetype = "dashed") +
  
  #4f + 1r
  geom_segment( x = 2026, xend = 2044, y = 40, yend = 40, color = "blue", size = 2) + 
  geom_segment( x = 2504, xend = 2524, y = 40, yend = 40, color = "red", size = 2) + 
  geom_segment( x = 2026, xend = 2524, y = 40, yend = 40, color = "black", size = 0.25, linetype = "dashed") +
  
  
  theme_bw() + 
  ylab("Median coverage (1/1000)") + xlab("") +  
  theme( panel.grid.minor   = element_blank(), 
         panel.grid.major.y = element_blank(), 
         axis.title.x=element_blank() ) +
  scale_x_continuous( limits = c(800, 2550), breaks = seq(800, 2500, by = 100), labels = seq(800, 2500, by = 100) ) 
#  scale_y_continuous( breaks = seq(0, 40000, by = 5000), labels = seq(0, 40000, by = 5000) )


# check SNV call biased by depth 
# 1 position
variant              <- read.csv( "./variant.csv", stringsAsFactors = FALSE)
tab_variant.pos      <- data.frame( table(variant$pos), stringsAsFactors = FALSE)
tab_variant.pos$Var1 <- as.integer( as.character(tab_variant.pos$Var1) )

d2 <- 
ggplot( data = tab_variant.pos ) + 
  theme_bw() +
  xlab("Position") + ylab("Freq. of variant") +
  geom_line( aes( x = Var1, y = Freq) ) + 
  theme( panel.grid.minor   = element_blank(), 
         panel.grid.major.y = element_blank() ) + 
  scale_x_continuous( limits = c(800, 2550), breaks = seq(800, 2500, by = 100), labels = seq(800, 2500, by = 100) ) 

ggarrange( d1, d2, ncol = 1, nrow = 2) #a5


# 2 sample-wise and methods used
variantSum   <- read.csv( "./variantSum.csv", stringsAsFactors = FALSE)
avg.coverage <- sapply( ls_read, function(x){ median(x[which(x != 0)]) } )
snv_coverage <- data.frame( snv = variantSum$SNV_n, depth = avg.coverage, pi = variantSum$pi, SE = variantSum$SE)

summary( lm( depth ~ snv, data = snv_coverage) )
summary( lm( depth ~ SE, data = snv_coverage) )
summary( lm( depth ~ pi, data = snv_coverage) )

d3 <- 
ggplot( data = snv_coverage, aes( x = depth, y = snv ) ) + geom_point() + 
  theme( axis.title.x = element_blank() ) +
  ylab("No. of variants") +
  annotate( "text", label = "R^{2}==0.25 ", parse = TRUE, x = 7500, y = 70, size = 4, colour = "black") 
d4 <-
ggplot( data = snv_coverage, aes( x = depth, y = SE ) ) + geom_point() + 
  theme( axis.title.x = element_blank() ) +
  ylab("Shannon Entropy") + xlab("Median coverage") +
  annotate( "text", label = "R^{2}==0.07 ", parse = TRUE, x = 7500, y = 3.5, size = 4, colour = "black")
d5 <-
ggplot( data = snv_coverage, aes( x = depth, y = pi*10000 ) ) + geom_point() + 
  ylab("Genetic diverty (pi x10000)") + xlab("Median coverage") +
  annotate( "text", label = "R^{2}==0.04 ", parse = TRUE, x = 7500, y = 12.5, size = 4, colour = "black")

ggarrange( d3, d4, d5, ncol = 1, nrow = 3) #a6
