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
  
  
  
  
  
  
  
  
  