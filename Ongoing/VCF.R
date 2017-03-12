library(stringr)

      E_nt <- seq(937, 2421)
  Primer_f <- c(seq(884, 901), seq(1310, 1327), seq(1616, 1635), seq(2026, 2044) )
  Primer_r <- c(seq(1357, 1376), seq(1752, 1770), seq(2035, 2053), seq(2504, 2524) )
  Primer <- unique(c(Primer_f, Primer_r))
  
    lsvcf <-  list.files("./vcfresult/")
  
  lsvcf_f <- which(grepl("f.vcf", lsvcf) == TRUE)
  lsvcf_r <- which(grepl("r.vcf", lsvcf) == TRUE)
  lsvcf_w <-setdiff(  setdiff( seq(1, length(lsvcf)), lsvcf_f), lsvcf_r)
  
  vcf0 <- read.table(paste0("./vcfresult/", lsvcf[lsvcf_w[i]]) )
  vcf1 <- read.table(paste0("./vcfresult/", lsvcf[lsvcf_f[i]]) )
  vcf2 <- read.table(paste0("./vcfresult/", lsvcf[lsvcf_r[i]]) )
  
  # in E but not in primer
  gf1 <- setdiff(which(vcf0[,2] %in% E_nt == TRUE), 
                 which(vcf0[,2] %in% Primer == TRUE) )
  
  # forward in primer_r and in E 
  gf2 <- intersect(which(vcf1[,2] %in% Primer_r == TRUE),
                   which(vcf1[,2] %in% E_nt == TRUE) )
  
  # reverse in primer_f and in E
  gf3 <- intersect(which(vcf2[,2] %in% Primer_f == TRUE),
                   which(vcf2[,2] %in% E_nt == TRUE) )
  
  
  vcf0[gf1, ]
  
  
  df <- data.frame(drug1 = c("B+A", "B", "A+C"))
  df
  df$drug2 <- lapply( strsplit(as.character(df$drug1), "\\+"), "[", 2)
  df$drug1 <- lapply(strsplit(as.character(df$drug1), "\\+"), "[", 1)
  df
  
  
  
  
  
  
  
  