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

############################################################
###### Site-specific error rate ############################

# s.RR 
# er0728 (with aa position)

# plot(nn, s.RR, pch=20, xaxt='n',xlab="Position", ylab = "Error Rate", 
#     xlim=c(900, 2550), ylim=c(0, 0.07))
# lines(nn,s.RR)
# axis(1, at=bb,labels=bb)
# text(nn[which(s.RR > 0.01)], s.RR[which(s.RR > 0.01)], nn[which(s.RR > 0.01)], pos=3, cex=0.5, col="blue")


### slow method (and has bug) #####
for (i in 904: 2505){
  
  s.VR = c()
  s.TR = c()
  
  for(k in 1: length(ls)){
  
    
    dddf = read.csv(ls[k])
    
      t0 = which(dddf[,i] == 0)
      cc = which.max(dddf[,i][1:4])
      tcc = sum(dddf[,i])                       #total counts
      td = c(cc, t0)
      
      t = dddf[,i][1:4][-td]
      
      if (length(t) == 0){ s.RR[ length(s.RR) + 1] = 0 } else {
      
      s.VR = c(s.VR, t)  
      s.TR = c(s.TR, rep(tcc, length(t)))
      s.mm = lm (s.VR ~ s.TR)
      bt = s.mm$coefficients[2]
      s.RR[ length(s.RR) + 1] = bt / (1 + bt)}
  
}
  print(i)

  } ## too slow
  


#### modified version (site specific and pooled) ################################

cSUM = do.call("rbind", lapply(ls, read.csv)) # dim(cSUM)


#site-specific
ER1=c()           # mean(ER1) = 0.005243203
ER0=c()           # mean(ER0) = 0.006933332

#pooled
ERp1 # 0.008472968
  cERp1=c()
  tERp1=c()
ERp0 # 0.01125579
  cERp0=c()
  tERp0=c()


for(i in 904: 2505){                                  #real position: V902 - V2503
  
  LSi = split(cSUM[,i], rep(1:221, each = 5))
  
  VR1.c = unlist(lapply(LSi, ff1c), use.names = FALSE)
  VR1.t = unlist(lapply(LSi, ff1t), use.names = FALSE)
  
  VR0.c = unlist(lapply(LSi, ff0c), use.names = FALSE)
  VR0.t = unlist(lapply(LSi, ff0t), use.names = FALSE)
  
  if( length(which(VR1.c != 0)) == 0){

    ER1[length(ER1) + 1] = 0
    ER0[length(ER0) + 1] = 0
    
    
  }else{
    
    mm1 = lm (VR1.c ~ VR1.t)
    mm0 = lm (VR0.c ~ VR0.t)
    
    b1 = mm1$coefficients[2]
    ER1[length(ER1) + 1] = (b1 / (b1 + 1))
    
    b0 = mm0$coefficients[2]
    ER0[length(ER0) + 1] = (b0 / (b0 + 1))
    
  }
  
print(i)  
}

# plot 

nn = seq(902,2503)
plot(nn, ER0, pch=20, xaxt='n',xlab="Position", ylab = "Error Rate", 
          xlim=c(900, 2550), ylim=c(0, 0.07))
      lines(nn,ER0)
      axis(1, at=bb,labels=bb)
      text(nn[which(ER0 > 0.02)], ER0[which(ER0 > 0.02)], nn[which(ER0 > 0.02)], pos=3, cex=0.5, col="Blue")
     
      
### pooled re-accessment 
      
      for(i in 904: 2505){                             #real position: 902 - 2503
        
        LSi = split(cSUM[,i], rep(1:221, each = 5))
        
        VR1.c = unlist(lapply(LSi, ff1c), use.names = FALSE)
        VR1.t = unlist(lapply(LSi, ff1t), use.names = FALSE)
        
        VR0.c = unlist(lapply(LSi, ff0c), use.names = FALSE)
        VR0.t = unlist(lapply(LSi, ff0t), use.names = FALSE)
        
        if( length(which(VR1.c != 0)) == 0){
          
        }else{
            
          cERp1 = c(cERp1, VR1.c)
          tERp1 = c(tERp1, VR1.t)
          
          cERp0 = c(cERp0, VR0.c)
          tERp0 = c(tERp0, VR0.t)

        }
        print(i)  
      }
      
      mm1 = lm (cERp1 ~ tERp1)
      mm0 = lm (cERp0 ~ tERp0)
      b1 = mm1$coefficients[2]
      b0 = mm0$coefficients[2]
      
      ERp1 = (b1 / (b1 + 1))
      ERp0 = (b0 / (b0 + 1))
      
      
      
      
############################################################
###### to calculate systemic error rate ###### AA ##########
############################################################


vRa  = c() # variant/shadow reads
TRa  = c() # total reads of respective position
aRa  = c() # all reads
mTRa = c() # total reads for all

for(k in 1: length(ls)){
  
  dddf <- read.csv(ls[k])
  
  for(i in 303: 835){              #region actually followed by primers
    
    saa = c(1:23)[-(19:21)]        #19- 21 : "*", "-", "X"
    
    
    tcc = sum(dddf[,i][saa])         #total counts
    
    if (tcc != 0){
      
      aRa = c(aRa, dddf[,i][saa])  
      mTRa = c(mTRa, rep(tcc, 20))

          cc = which.max(dddf[,i][saa])
          
      vRa = c(vRa, dddf[,i][saa][-cc])
      TRa = c(TRa, rep(tcc, 19))
      
    }
    
  }
  
  print((k/length(ls)*100))
  
}

rrR = vRa/TRa 
plot(density(vRa))
polygon(density(rrR), col = "red", border = "gray") #distribution of all values
quantile(rrR, 0.99) # 0.08 = Q99.9; 0.0365 = Q99

plot(TRa, vRa, pch=20)
mm = lm (vR ~ TR)
summary(mm) #beta = 0.001139 Error Rate = beta / ( 1 + beta ) = 0.0113
abline(mm)

######################################################################################
# examine the distribution of the 2nd most variant % in each position in all samples #
######################################################################################


allsend.r0 = sndvariantp(lsfull, 1)           # max(allsend.r0)  
allsend.r500 = sndvariantp(lsfull, 2)         # max(allsend.r500) 

c=c(allsend.r0, allsend.r500)
cc=c(rep("0", length(allsend.r0)), rep("500", length(allsend.r500)))
a=data.frame(c, cc) 
names(a) = c("position", "reads")

aa=data.frame(allsend.r500, rep(500, length(allsend.r500)))
names(aa) = c("position", "reads")


#dis
p<-ggplot(a, aes(position, fill = reads)) + geom_density(alpha = 0.2)


#jitter
j=ggplot(a, aes(reads, position)) + geom_jitter(alpha=I(1/4), aes(color=position))

#pirate  
library(yarrr)
p=pirateplot(formula = position ~ reads , data = aa, theme.o = 3, 
             main="", inf = "ci",
             ylab="", gl.col = gray(.8),
             pal="basel")

box(which = "p")
#box  
b=ggplot(a, aes(reads, position, fill=reads)) + geom_boxplot(aes(fill=reads)) + ylab("Position") +
  theme(legend.position = "none") + ylim(0,0.05)

#t.test for CI

t.test(allsend.r0) # mean = 0.01565829 CI = 0.01555479, 0.01576180 
t.test(allsend.r500) # mean = 0.01418616 CI = 0.01410864, 0.01426368

sd(allsend.r0) 
0.01565829 + 2.58*(  sd(allsend.r0) )/sqrt(length(allsend.r0)) #set to 99% = ~ 0.0158
