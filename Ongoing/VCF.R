# call variant number, and caculate pi

library(stringr)

  # primer and E protein regions

      E_nt <- seq(937, 2421)
  Primer_f <- c(seq(884, 901), seq(1310, 1327), seq(1616, 1635), seq(2026, 2044) )
  Primer_r <- c(seq(1357, 1376), seq(1752, 1770), seq(2035, 2053), seq(2504, 2524) )
    Primer <- unique(c(Primer_f, Primer_r))

  # read-in vcf files
    
      lsvcf <- list.files("./vcfresult/")
    lsvcf.o <- c(lsvcf[388:414], lsvcf[1:387])
  
  # list of whole, forward-only or reverse-only file lists
    
  lsvcf_f <- which(grepl("f.vcf", lsvcf.o) == TRUE)
  lsvcf_r <- which(grepl("r.vcf", lsvcf.o) == TRUE)
  lsvcf_w <-setdiff(  setdiff( seq(1, length(lsvcf.o)), lsvcf_f), lsvcf_r)
  
  vcf0 <- read.table( paste0("./vcfresult/", lsvcf.o[lsvcf_w[i]]) )
  vcf1 <- read.table( paste0("./vcfresult/", lsvcf.o[lsvcf_f[i]]) )
  vcf2 <- read.table( paste0("./vcfresult/", lsvcf.o[lsvcf_r[i]]) )
  
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
  
  
  # forward ref alleles / reverse ref / forward non-ref / reverse non-ref alleles
  
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
    v_freq <- x_v / x_sum
     pi_i <- ((x_v)*(x_sum - x_v))/ ( ((x_sum)*(x_sum - 1))/2 )
      
    return(y = c(v_freq, pi_i, x_v, x_ref, x_sum))  
    })
         
  y <- cbind(gf_u, t(x2))[, c(2, 4, 5, 6, 9, 10, 11, 12, 13)]
  colnames(y) <- c("pos", "ref", "alt", "qual","alt_freq", "pi", "alt_n", "ref_n", "n")
  
  
  
  
  
  
  
# create consensus seq for each vcf ####
  
   ls <- list.files("./allSUM/")
  ls0 <- gsub(".csv", "",gsub("SUM", "", ls))
  
  # aim: to generate .fas of consensus sequence for each sample 
  
  library(seqinr)
  
  # use 2001 consensus as backbone
  
  ff <- read.fasta("./2001Con.fas")
  ff.seq <- getSequence(ff)
  ff.seq.u <- unlist(ff.seq)
  
  for (i in 1:length(ls)){
    
    # readin summarized file (.csv) derived from .bam
    
    ff.i <- read.csv(paste0("./allSUM/", ls[i] ) )
    lth <- length(ff.i)
    cons <- c()
    
    # get consensus nucleotide for each position
    
    for(k in 3: lth){
      
      tcc = sum( ff.i[1,k] + ff.i[2,k] + ff.i[3,k] + ff.i[4,k])   
      
      if ( ( ff.i[1,k] | ff.i[2,k] | ff.i[3,k] | ff.i[4,k] | ff.i[5,k]) != 0 ){
        
        aamatrix <- c("a","t","c","g")
        
        nmax <-  which.max( c( ff.i[1,k], ff.i[2,k], ff.i[3,k], ff.i[4,k]) )
        
        cons[ length(cons) +1 ] <-  aamatrix[nmax]  
        
      }else{
        
        cons[length(cons) +1 ] = "-" }
    }
    
    # refill gap; original summarized file comprises 1-3000 in ORF
    
    cons <- c(cons, rep("-", length(ff.seq.u)-3000) )
    
    cons.u <- c()
    
    for (l in 1: length(ff.seq.u)){
      
      if (cons[l] == "-"){
        
        cons.u[l] <- ff.seq.u[l]
        
      }else{
        
        cons.u[l] <- cons[l]
      }
    }
    
    write.fasta(list(cons.u), file.out = paste0("./confas/", ls0[i], ".fasta"), 
                
                names = gsub("/", "_", attributes(ff)$names, ) )
    
  }
  