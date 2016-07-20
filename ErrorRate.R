#########################################
###### to calculate systemic error rate #

setwd("/Users/yaosmacbook/twDENV/allSUM")
ls = list.files(getwd())

vR = c() # variant/shadow reads
TR = c() # total reads of respective position
aR = c() # all reads
mTR = c() # total reads for all

for(k in 1: length(ls)){

  dddf <- read.csv(ls[k])
  lth = length(ls[k])
  
  for(i in 902: 2502){  #region actually followed by primers

    tcc = sum(dddf[1,i] + dddf[2,i] + dddf[3,i] + dddf[4,i])  #total counts
    
    aR = c(aR, dddf[1,i], dddf[2,i], dddf[3,i], dddf[4,i])  
    mTR = c(mTR, rep(tcc, 4))
    
    cc <- c(dddf[1,i], dddf[2,i], dddf[3,i], dddf[4,i])
    mv = which.max(cc)
    
    vR = c(vR, cc[-mv])
    TR = c(TR, rep(tcc, 3))
    
  }

  print((k/length(ls)*100))
  
}


