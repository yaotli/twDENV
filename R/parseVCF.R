require(stringr)
require(ggplot2)

# primer and E region

E    <- seq( 937, 2421 )
P_f  <- c( seq(884, 901), seq(1310, 1327), seq(1616, 1635), seq(2026, 2044) )
P_r  <- c( seq(1357, 1376), seq(1752, 1770), seq(2035, 2053), seq(2504, 2524) )
P_fr <- unique( c( P_f, P_r ) )

# read-in

ls_vcf0  <- grep( "[0-9].vcf", paste0( "./vcf/", list.files( "./vcf/" ) ), value = TRUE )

# parse the vcf file ----

variant <- c()
for( i in 1: length(ls_vcf0) )
{
  vcf0 <- read.table( ls_vcf0[i] )
  vcfF <- read.table( gsub( "\\.vcf", "_f.vcf", ls_vcf0[i] ) )
  vcfR <- read.table( gsub( "\\.vcf", "_r.vcf", ls_vcf0[i] ) ) 
  
  # vatiants in E but not in primer
  i_0 <- setdiff( which( vcf0[,2] %in% E == TRUE), which( vcf0[,2] %in% P_fr == TRUE ) )
  
  # read variants in `reverse` primers region with `forward` reads
  i_r <- setdiff( intersect( which( vcfF[,2] %in% P_r == TRUE), which( vcfF[,2] %in% E == TRUE ) ), which( vcfF[,2] %in% P_f == TRUE) )
  
  # read variants in `forward` primers region with `reverse` reads
  i_f <- setdiff( intersect( which( vcfR[,2] %in% P_f == TRUE), which( vcfR[,2] %in% E == TRUE ) ), which( vcfR[,2] %in% P_r == TRUE) )
  
  # combine the eligible variants
  tab <- rbind( vcf0[i_0, ], vcfF[i_r, ], vcfR[i_f,] )
  tab <- tab[ order(tab[,2]), ]
  
  if ( TRUE %in% duplicated(tab[,2]) ){ print( ls_vcf0[i] ) }
  
  dp4 <- sapply( lapply( strsplit(as.character(tab$V8), ";"), "[", 4), 
                 function(x){ y = str_match(x, "=([0-9,]+)")[2] } )
  
  app <- sapply( strsplit( dp4, ","), 
                 function(x){ 
                   x      <- as.numeric(x)
                   x_sum  <- sum(x)
                   x_ref  <- sum(x[1:2])
                   x_v    <- sum(x[3:4])
                   v_freq <- x_v / x_sum
                   return( y = c( v_freq, x_v, x_ref, x_sum ) )
                 })
  
  y     <- cbind( tab, t(app) )[, c( 1, 2, 4, 5, 6, 9, 10, 11, 12 ) ]
  y[,1] <- str_match( ls_vcf0[i], "([0-9]+)." )[,2]
  
  variant = rbind( variant, y )
}

colnames(variant) <- c("sample", "pos", "ref", "alt", "qual", "v_freq", "v_n", "ref_n", "sum_n")
write.csv(variant, "./variant.csv", row.names = FALSE)

# summary vcf outcome ----

SNV_n <- c()
pi    <- c()
SE    <- c()

for( i in 1: length(unique(variant$sample)) )
{
  k     <- unique( unique(variant$sample) )[i]
  tab_k <- variant[ which( variant$sample == k), ]
  tab_k <- tab_k[ which( tab_k$sum_n >= 1000), ]
  
  if( TRUE %in% duplicated(tab_k$pos) )
  {
    
    dup_n <- length( which( duplicated(tab_k$pos) == TRUE) )
    dup_p <- tab_k$pos[ which( duplicated(tab_k$pos) ) ]

    # adjust for multiple SNV calls at one position
    
    pi.k.j <- c()
    se.k.j <- c()
    for( j in 1: length(dup_n) )
    {
      tab_t  <- tab_k[  which( tab_k$pos == dup_p[j] ), ]
      new_n  <- sum( tab_t$v_n ) + tab_t$ref_n[1] 
      
      pi.k.j <- c(pi.k.j, 
                  ( sum(tab_t$ref_n * tab_t$v_n) + prod( tab_t$v_n ) ) / ( ( new_n -1 )*( new_n )/2 ) )
      se.k.j <- c(se.k.j,
                  sum( (-1)*( tab_t$v_n/new_n )*log( tab_t$v_n/new_n ) ) + (-1)*( tab_t$ref_n[1]/new_n )*log( tab_t$ref_n[1]/new_n ) )
    }
    
    tab_k  <- tab_k[ -grep( paste(dup_p, collapse = "|"), tab_k$pos ), ]
  
    pi.k.1 <- sum( (tab_k$ref_n * tab_k$v_n) / ( (tab_k$sum_n - 1)*(tab_k$sum_n)/2 ) )
    pi.k   <- (pi.k.1 + sum(pi.k.j) )/ (2421-937+1) 
    
    se.k.1 <- sum(-( tab_k$v_freq*log(tab_k$v_freq) + (1-tab_k$v_freq)*log( (1-tab_k$v_freq) )  ))
    se.k   <- se.k.1 + sum(se.k.j)
    
    SNV.k  <- length( tab_k$pos ) + length( tab_t$pos )
      
  }else
  {
    pi.k  <- sum(  (tab_k$ref_n * tab_k$v_n) / ( (tab_k$sum_n - 1)*(tab_k$sum_n)/2 ) )/ (2421-937+1)
    se.k  <- sum( -( tab_k$v_freq*log(tab_k$v_freq) + (1-tab_k$v_freq)*log( (1-tab_k$v_freq) ) ) )
    SNV.k <- length( tab_k$pos )
    }
  
  SNV_n <- c(SNV_n, SNV.k)
  pi    <- c(pi, pi.k)
  SE    <- c(SE, se.k)
}

variantSum <- data.frame( sample = unique(variant$sample), SNV_n, pi, SE)
write.csv( variantSum, "./variantSum.csv", row.names = FALSE)




