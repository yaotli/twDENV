library(stringr)

  # primer and E protein regions

      E_nt <- seq(937, 2421)
  Primer_f <- c(seq(884, 901), seq(1310, 1327), seq(1616, 1635), seq(2026, 2044) )
  Primer_r <- c(seq(1357, 1376), seq(1752, 1770), seq(2035, 2053), seq(2504, 2524) )
    Primer <- unique(c(Primer_f, Primer_r))

  # read-in vcf files
    
    lsvcf <-  list.files("./vcfresult/")
  
  # list of whole, forward-only or reverse-only file lists
    
  lsvcf_f <- which(grepl("f.vcf", lsvcf) == TRUE)
  lsvcf_r <- which(grepl("r.vcf", lsvcf) == TRUE)
  lsvcf_w <-setdiff(  setdiff( seq(1, length(lsvcf)), lsvcf_f), lsvcf_r)
  
  vcf0 <- read.table(paste0("./vcfresult/", lsvcf[lsvcf_w[i]]) )
  vcf1 <- read.table(paste0("./vcfresult/", lsvcf[lsvcf_f[i]]) )
  vcf2 <- read.table(paste0("./vcfresult/", lsvcf[lsvcf_r[i]]) )
  
  # vatiants in E but not in primer
  gf1 <- setdiff(which(vcf0[,2] %in% E_nt == TRUE), 
                 which(vcf0[,2] %in% Primer == TRUE) )
  
  # forward in primer_r and in E 
  gf2 <- interaction(which(vcf1[,2] %in% Primer_r == TRUE),
                   which(vcf1[,2] %in% E_nt == TRUE) )
  
  # reverse in primer_f and in E
  gf3 <- intersect(which(vcf2[,2] %in% Primer_f == TRUE),
                   which(vcf2[,2] %in% E_nt == TRUE) )

  # combine the eligible variants
  gf_u <- rbind(vcf0[gf1, ], vcf1[gf2, ], vcf2[gf3,])
  
  
  x1 <- sapply(lapply(strsplit(as.character(gf_u$V8), ";"), "[", 4), 
         function(x){
           y = str_match(x, "=([0-9,]+)")[2]
           return(y)
         })
      
  x2 <- sapply(strsplit(x1, ","), function(x){
    
    x = as.numeric(x)
    x_sum <- sum(x)
    x_ref <- sum(x[1:2])
      x_v <- sum(x[3:4])
    v_feq <- x_v / x_sum
     pi_i <- ((x_v)*(x_sum - x_v))/ (((x_sum)*(x_sum - 1))/2)
      
    return(y = c(v_feq, pi_i, x_v, x_ref, x_sum))  
    })
         
  y <- cbind(gf_u, t(x2))
  
  
  
  
  
  
  
  