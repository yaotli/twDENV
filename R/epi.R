require(tidyverse)

# read-in
epi_csv    <- read.csv("./epi_data_csv.csv", stringsAsFactors = FALSE)
variantSum <- read.csv("./variantSum.csv", stringsAsFactors = FALSE)
sec_csv    <- read.csv("./20180329variant.csv", stringsAsFactors = FALSE)

ix      <- match( variantSum$sample, epi_csv$ID )
ix2     <- match( variantSum$sample, sec_csv$ID )

epi_var <- data.frame( sample = variantSum$sample,
                       snv_n  = variantSum$SNV_n, 
                       pi     = variantSum$pi,
                       se     = variantSum$SE, 
                       clade  = epi_csv$clade[ix],
                       time   = epi_csv$IsoD[ix],
                       onsetD = epi_csv$TpIllness[ix],
                       secInf = sec_csv$Secondary[ix2], stringsAsFactors = FALSE )

epi_var$geoix <- ifelse( variantSum$sample %in% sec_csv$ID, 1, 0)



# clade_variant

epi_var %>% 
  filter( geoix == 1 ) %>%
  
ggplot( aes(x = secInf , y = snv_n, group = secInf) ) + theme_bw() +
  geom_violin( position = position_dodge(0.8), draw_quantiles = c(0.5), size = 1) + 
  geom_jitter( aes( group = secInf), position = position_dodge(0.8), alpha = 0.4, size = 3) +  
  theme( axis.title.x = element_text(size = 18),
         axis.title.y = element_text(size = 16),
         axis.text = element_text(size = 15)) 

ggplot( data = epi_var ) + geom_point( aes( x = time, y = se) ) 
