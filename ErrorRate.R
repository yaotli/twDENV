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
  
  for(i in 904: 2505){  #region actually followed by primers

    tcc = sum(dddf[1,i] + dddf[2,i] + dddf[3,i] + dddf[4,i])  #total counts
    
    if (tcc != 0){
      
      aR = c(aR, dddf[1,i], dddf[2,i], dddf[3,i], dddf[4,i])  
      mTR = c(mTR, rep(tcc, 4))
      
      cc <- c(dddf[1,i], dddf[2,i], dddf[3,i], dddf[4,i])
      mv = which.max(cc)
      
      vR = c(vR, cc[-mv])
      TR = c(TR, rep(tcc, 3))
      
    }

  }

  print((k/length(ls)*100))
  
}

rrR = aR/mTR 
plot(density(rrR))
polygon(density(rrR), col = "red", border = "gray") #distribution of all values

plot(TR, vR, pch=20)
mm = lm (vR ~ TR)
summary(mm) #beta = 0.001139 Error Rate = beta / ( 1 + beta ) = 0.0113
abline(mm)



