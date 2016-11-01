########## General Data preparation ##############
#

ls0 = list.files(getwd())                                 # allSUM folder files name
  allSUM = do.call("rbind", lapply(ls0, read.csv))        # n = 221

ls.raw                                                    # ordered file name (n = 221)
  ls.clean                                                # ordered ID (n = 221, ID only)
  
no.c0                   # no. of c0 in ls.raw             # n = 151
  ls.raw.c0             # C0 files name in ls.raw
  ls.clean.c0 = ls.clean[no.c0]
  allSUM.c0 = do.call("rbind", lapply(ls.raw.c0, read.csv))
  
no.c1
  ls.raw.c1                                               # n = 70
  ls.clean.c1 = ls.clean[no.c1]
  allSUM.c1 = do.call("rbind", lapply(ls.raw.c1, read.csv))
 
##### Probe out all variant by Poisson method #####
#
# variantPT 
# lsvariantPT (df, pc, nr, adjER, listname)

adjER = rep(0.001, 3000)

  # out.lsv.0.01.1000r = lsvariantPT(allSUM, 0.01, 3, adjER, ls0)
  out.lsv.0.001.1000r = lsvariantPT(allSUM, 0.001, 3, adjER, ls0)  # n = 221, not in order
  
  # ID DV1461 (71), DV1464 (73): 0 variants
  
      df.out.lsv.0.01.1000r = do.call("rbind", out.lsv.0.01.1000r) 
      df.out.lsv.0.001.1000r = do.call("rbind", out.lsv.0.001.1000r)
  
  out.lsv.0.001.1000r.c0 = lsvariantPT(allSUM.c0, 0.001, 3, adjER, ls.clean.c0) # ordered 
  out.lsv.0.001.1000r.c1 = lsvariantPT(allSUM.c1, 0.001, 3, adjER, ls.clean.c1)
  
    #  id = sapply(strsplit(attributes(df.out.lsv.0.01.1000r)$row.names, split = ".", fixed = T), function(x) head(x,1))
    #  id = sapply(strsplit(attributes(df.out.lsv.0.001.1000r)$row.names, split = ".", fixed = T), function(x) head(x,1))

          df.out.lsv.0.01.1000r = cbind(id, df.out.lsv.0.01.1000r)
          df.out.lsv.0.001.1000r = cbind(id, df.out.lsv.0.001.1000r)

##### plot with position #####
#
      library(ggplot2)
      bk = seq(900, 2500, by=50)  
          
      which ( df.out.lsv.0.001.1000r$pVarN == 0 )
      
    
ggplot( df.out.lsv.0.001.1000r[-c(656,662),] , aes(sVarN, (pVarN), color= id)) + geom_point() + 
  
      theme_bw() + theme(legend.position="none") + scale_x_continuous(breaks=bk) + 
  
      xlab("Position") + ylab("Proportion") 
  

##### Variant disribution #####
#

out.v.m3e.3er = as.data.frame(table(df.out.lsv.0.001.1000r$sVarN))

      # out.v.m3e.3er.OD = out.v.m3e.3er[order(-out.v.m3e.3er$Freq),]
      
ggplot(out.v.m3e.3er[-1,], aes(x=as.numeric(as.character(Var1)), y=Freq)) +

      geom_point() + theme_bw() + scale_x_continuous(breaks=bk) + xlab("Position") +
  
      geom_line()
  

##### heatmap #####
#
# 1 list 

ls.raw = c(ls0[213:221], ls0[1:212])                # ordered ls0

ls.clean = gsub("SUM", replacement = "", ls.raw)
ls.clean = gsub(".csv", replacement = "", ls.clean) # ls.claen

c1 = sapply(strsplit(ls.clean, split = "_", fixed = T), function(x){length(x) >1})
no.c1 = which(c1 == "TRUE")  # n = 70
no.c0 = which(c1 == "FALSE") # n = 151, therefore

ls.raw.c0 = ls.raw[no.c0]
ls.raw.c1 = ls.raw[no.c1]

#
# 2 lspvariantPT (df, pc, nr, kd, b, adjER, filenames)

breaks = seq(930, 2430, by=5)


out.p.0.001.1000r.c0 = lspvariantPT(allSUM.c0, 0.001, 3, 0, breaks, adjER, ls.clean.c0)
out.p.0.001.1000r.c1 = lspvariantPT(allSUM.c1, 0.001, 3, 0, breaks, adjER, ls.clean.c1)
  
#
# 3 ggplotting 

breaks = seq(930, 2430, by=5)                 # length(breaks) = 301
ticks<-out.p.0.001.1000r.c0$region[1:300]     # originally what showed
ticksss<-ticks[seq(5,295, by=10)]             # the place we preferred to replace
ticklab<-seq(950, 2400, by=50)

# c0
  ggplot(out.p.0.001.1000r.c0, aes(region, no.sample)) + 
    
    geom_tile(aes(fill=Freq), colour="white") +  theme(legend.position="none") +
  
    scale_fill_gradient(low="white", high="Blue") + 
  
  theme(axis.text.x=element_text(angle=90, vjust = 0.5)) +
    
    scale_y_continuous(breaks = c(1:151), labels = ls.clean.c0) + 
    
    xlab("Position") + ylab("Sample") + scale_x_discrete(breaks=ticksss, labels=ticklab)

# c1 
  ggplot(out.p.0.001.1000r.c1, aes(region, no.sample)) + 
    
    geom_tile(aes(fill=Freq), colour="white") +  theme(legend.position="none") +
    
    scale_fill_gradient(low="white", high="green") + 
    
    theme(axis.text.x=element_text(angle=90, vjust = 0.5)) +
    
    scale_y_continuous(breaks = c(1:70), labels = ls.clean.c1) + 
    
    xlab("Position") + ylab("Sample") + scale_x_discrete(breaks=ticksss, labels=ticklab)
  



##### Close examination #####

library(dplyr)

 a=     df.out.lsv.0.001.1000r %>%
        
        select(pVarN)


##### Fit to Distribution #####
 
library(fitdistrplus)

plotdist(a$pVarN, histo = T, demp = T)
descdist(a$pVarN, boot = 1000)  

fw<-fitdist(a$pVarN, "weibull")  
fln<-fitdist(a$pVarN, "lnorm") 
fg<-fitdist(a$pVarN, "gamma")    # summary(fln): meanlog -4.642470; sdlog 1.561318
fp<-fitdist(a$pVarN, "pois", method = "mme")

image(as.matrix(leg),col=cx,axes=T)
par(mfrow = c(2,2))
plot.legend = c("Weibull", "lognormal", "Gamma")

denscomp(list(fw, fln, fg), legendtext = plot.legend)
qqcomp(list(fw, fln, fg), legendtext = plot.legend)
cdfcomp(list(fw, fln, fg), legendtext = plot.legend)
ppcomp(list(fw, fln, fg), legendtext = plot.legend)

# Gamma and Lognormal

  x = rgamma(100000, shape = 0.4315147, rate = 9.9349651)
  y = rlnorm(100000, meanlog = -4.642470, sdlog = 1.561318)
  plot(density(a$pVarN), lwd = 2) 
  lines(density(x), col="Red", lwd=2)
  
  qgamma(0.05, shape = 0.4315147, rate = 9.9349651, lower.tail = FALSE)
  # 0.175765
  
  qlnorm(0.05, meanlog = -4.642470, sdlog = 1.561318, lower.tail = FALSE)
  # 0.1256405
  

########## C0 sample vs C1 ##########
  
# out.lsv.0.001.1000r.c0 
# out.lsv.0.001.1000r.c1 
# out.n.c0 / ls.clean.c0
# out.n.c1 / ls.clean.c1
# out.n.c0.p
# out.n.c1.p
# df.c0c1
# df.c0c1.p

    
# out.n.c1 = as.numeric( sapply(out.lsv.0.001.1000r.c0, function(x){
#          length (which(x$sVarN != 0) )   })) 
# df.c0c1 = data.frame(ID=c(ls.clean.c0, ls.clean.c1), no=c(out.n.c0, out.n.c1),
#                     c0c1=c(rep("c0", length(ls.clean.c0)), rep("c1", length(ls.clean.c1))))  

  
# Pirate plot
#

install.packages("devtools")
library("devtools")
install_github("ndphillips/yarrr")
library(yarrr)  
  

p=pirateplot(formula = no ~ c0c1, data = df.c0c1, theme.o = 3, 
             main="All samples available", inf = "ci", inf.o = 0.3,
             xlab="", ylab="No. of Variants", gl.col = gray(.8),
             pal="basel")

box(which = "p")


# No contamination version (9/24)
#

# contami = as.character(c(5083, 5086, 5145, 5184, 5146, 5110, 5152, 5181, 5005, 5059, 5099, 5071,
#                         5111, 5183, 5102, 5177, 5118, 5196, 5009, 5151, 5180)) # n = 21
# no.contami = match(contami, ls.clean.c0)
# no.contami = sort(no.contami)
#
# df.c0c1 = data.frame(ID=c(ls.clean.c0, ls.clean.c1), no=c(out.n.c0, out.n.c1),
# c0c1=c(rep("c0", length(ls.clean.c0)), rep("c1", length(ls.clean.c1))))  



  ls.clean.c0.Dcon = ls.clean.c0[-no.contami]                       # n = 130 (151 - 21)
  out.n.c0.Dcon = out.n.c0[-no.contami]                             # note: DV5083
  
  df.c0c1.Dcon = data.frame(ID=c(ls.clean.c0.Dcon, ls.clean.c1), no=c(out.n.c0.Dcon, out.n.c1),
                       c0c1=c(rep("c0", length(ls.clean.c0.Dcon)), rep("c1", length(ls.clean.c1))))  
  
  
  p1=pirateplot(formula = no ~ c0c1, data = df.c0c1.Dcon, theme.o = 3, 
               main="All samples available", inf = "ci", inf.o = 0.6,
               xlab="", ylab="No. of Variants", gl.col = gray(.8),
               pal="basel")
  
  box(which = "p")
  

##### paired C0 vs C1 #####
#

# c1 = sapply(strsplit(ls.clean.c1, split = "_", fixed = T), 
#       function(x) head(x,1))

# out.n.c1.p = out.n.c1[-which(duplicated(c(ls.clean.c0, c1))[152:221] == "FALSE")]
#          a = c1[-which(duplicated(c(ls.clean.c0, c1))[152:221] == "FALSE")]
# out.n.c0.p = out.n.c0[which(duplicated(c(a, ls.clean.c0))[66:216] == "TRUE")]
# df.c0c1.p = data.frame(no=c(out.n.c0.p, out.n.c1.p),
#                       c0c1=c(rep("c0", length(out.n.c0.p)), rep("c1", length(out.n.c1.p))))  



# Pirate plot (paired sample)

## No contamination 
#

d = c(64, 129)
df.c0c1.p.Dcon = df.c0c1.p[-d,]  # 5083

  p = pirateplot(formula = no ~ c0c1, data = df.c0c1.p.Dcon, theme.o = 3, 
             main="Paired Samples (n = 64)", inf = "ci", inf.o = 0.3,
             xlab="", ylab="No. of Variants", gl.col = gray(.8),
             pal="basel")

  box(which = "p")


# High frequency variant 
#
# out.n.c0.p.hfv = out.n.c0.p.hfv[which(duplicated(c(a, ls.clean.c0))[66:216] == "TRUE")]
# out.n.c1.p.hfv = out.n.c1.p.hfv[-which(duplicated(c(ls.clean.c0, c1))[152:221] == "FALSE")]
# df.c0c1.p.hfv
# df.c0c1.p.hfv.pattern15

out.n.c0.p.hfv = as.numeric( sapply(out.lsv.0.001.1000r.c0, function(x){
  
  length (which(x$pVarN >= 0.175765) ) }))

df.c0c1.p.hfv = data.frame(no=c(out.n.c0.p.hfv, out.n.c1.p.hfv),
                       c0c1=c(rep("c0", length(out.n.c0.p.hfv)), rep("c1", length(out.n.c1.p.hfv))))  

p=pirateplot(formula = no ~ c0c1, data = df.c0c1.p.hfv, theme.o = 3, 
             main="All samples available", inf = "ci", inf.o = 0.6,
             xlab="", ylab="No. of Variants", gl.col = gray(.8),
             pal="basel")

box(which = "p")


# Freq change

    a=which(out.n.c0.p.hfv == 0)[-6]
    b=which(out.n.c1.p.hfv == 0)
    duplicated(c(b,a))

  df.c0c1.p.hfv.pattern15 = data.frame(no=c(out.n.c0.p.hfv[-a], out.n.c1.p.hfv[-a]), 
                                     sample=rep(1:16,2), 
                                     c0c1 = c(rep("C0",16), rep("C1", 16)))

## No Contamination 
#
# df.c0c1.p.Dcon


  p.c1mc0 = df.c0c1.p.Dcon$no[65:128] - df.c0c1.p.Dcon$no[1:64]
  df.p.c1mc0 = data.frame(p.c1mc0, id=1)

  p = pirateplot(formula = p.c1mc0 ~ id, data = df.p.c1mc0, theme.o = 3, 
           main="Variation change", inf = "ci", inf.o = 0.3,
           xlab="", ylab="No. of Variants", gl.col = gray(.8),
           pal="basel")

  box(which = "p")



########## Onset day and DHF ##########
#
# out.n.c0M
# c0.table


demo = read.csv(file.choose())

out.n.c0M = out.n.c0[-103]                    # to remove DENV5083 n = 150
c0.table = data.frame(demo, out.n.c0M)

# pirate


##### Eliminate Contamination (9/24)
#
# adjc0.table
#

# TpIllnessO = factor(adjc0.table$TpIllness, levels = c(1:14))
# c0.dhfonset.Dcon.2 = data.frame(adjc0.table, TpIllnessO)

# *** TpIllnessO is problematic; use f and ff instead (OnsetD) (10/22): 
# *** adjc0.table::DoSampling is incorrect; use p instead (IsoD): 
#
#

  f = as.character(adjc0.table$TpIllness)
  f[60] = 99
  f = as.numeric(f)
  
  ff = f
  ff[which(ff >= 5)] = 5
  ff[60] = 99
  OnsetD = ff
  
# 
  p = adjc0.table$DoSampling
  p[20:77] = p[20:77] + 365 
  p = p/365
  IsoD = p
#
  clade0 = read.csv(file.choose(), header = F)
  match(c0.dhfonset.Dcon.Oct$ID, clade0$V1)
  
  cl = c()
  cl[1:77] = as.character(clade0$V2[match(c0.dhfonset.Dcon.Oct$ID, clade0$V1)[1:77]])
  cl[78:130] = "III"
  clade = cl
   
#
  # c0.dhfonset.Dcon.Oct = data.frame(adjc0.table, OnsetD, IsoD, clade)

  c0.dhfonset.Dcon.Oct.0102 = 
    c0.dhfonset.Dcon.Oct %>%
    filter(DHF != ".") %>%
    filter(Outbreak == 1) %>%
    filter(OnsetD != 99) %>%
    select(ID, out.n.c0, DHF, OnsetD, IsoD, clade, Secondary)
  
p=pirateplot(formula = out.n.c0 ~ OnsetD + DHF, data = c0.dhfonset.Dcon.Oct.15, 
              theme  = 1, main="2015 ", inf.method = "hdi", inf.f.o = 0.2,
              ylab="No. variants", gl.col = gray(.8), 
              pal="basel")

box(which = "p")


##### Epitime and no.variant #####
#
# c0.table
# coEpitime

## Correction (10/22)
#
# c0.dhfonset.Dcon.Oct.0102
# c0.dhfonset.Dcon.Oct.15

a = ggplot(data = c0.dhfonset.Dcon.Oct.0102, 
           aes(x=IsoD, y=out.n.c0, color = clade)) + geom_point(size = 3) +
  
  ylab("No. Variant") + xlab("") + ggtitle("2001-2002") + 
  theme(axis.title=element_text(size=20), axis.text=element_text(size=18), 
        plot.title = element_text(size=25)) 

multiplot(a,b)


##### Secondary vs Primary Infection ##### 
#
# c0.dhfonset.Dcon.Oct.0102
# c0.dhfonset.Dcon.Oct.15

b=pirateplot(formula = out.n.c0 ~ Secondary, data = c0.dhfonset.Dcon.Oct.15, inf.method ="ci",
             theme  = 1, main="2015 ", inf.f.o = 0.1, inf.lwd = 0, bar.f.o = 0.1,
             point.cex = 2, avg.line.lwd = 10,
             ylab="No. variants", gl.col = gray(.8), 
             pal="basel")

box(which = "p")




########## Repeated Sample ##########


c0.rep = 
    c0.table %>%
      filter(Repeated != ".") %>%
      filter(Repeated != "Y") %>%
      select(out.n.c0M, Repeated, DHF, Outbreak, TpIllness)
  
ggplot(c0.rep, aes(x=TpIllness, y=out.n.c0M, color=DHF, group=Repeated)) + 
  geom_line() + facet_wrap(~ Repeated)


########## Region vs Onset day ##########
#
# out.p.0.001.1000r.c0 
# ls.clean.c0
# c0.table
# out.lsv.0.001.1000r.c0

contami = as.character(c(5083, 5086, 5145, 5184, 5146, 5110, 5152, 5181, 5005, 5059, 5099, 5071,
            5111, 5183, 5102, 5177, 5118, 5196, 5009, 5151, 5180))
no.contami = match(contami, ls.clean.c0)
no.contami = sort(no.contami)

    adjls.clean.c0 = ls.clean.c0[-no.contami]
    adjls.raw.c0 = ls.raw.c0[-no.contami]


co.site <-sapply(out.lsv.0.001.1000r.c0, function(x){ y = x$sVarN 
return(y) } )

  adjco.site=co.site[-no.contami]
  
# 1 Eliminate contamination sample

  demo151 = rbind(demo[1:102,], demo[103,], demo[103:150,])
  attributes(demo151)$row.names = c(1:151)
  demo151[103,1:9] = NA
  
  c0.table151 = data.frame(demo151, out.n.c0)
  adjc0.table = c0.table151[-no.contami,]
  attributes(adjc0.table)$row.names = c(1:130)
  adjc0.table = cbind(adjc0.table, no=c(1:130))
   
  
# 2 Base on the no to choose sample 
  
library(dplyr)  
  
 rg.15.ond = adjc0.table %>% 
    arrange(TpIllness) %>%
    filter(TpIllness != "." )  %>%
    filter(Outbreak == 2 )  %>%
    select(no, TpIllness)  
  
# 3 Heatmap preparation
#
# breaks = seq(930, 2430, by=5)       
 
input = rg.15.ond
output = data.frame()
  
  t.1 = as.data.frame(table(as.numeric(as.character(input$TpIllness))))
  
  for(i in 1: dim(t.1)[1]){
    
    n.1 = t.1$Freq[i]

    h = sum(t.1$Freq[1:i])
    k = h - n.1 + 1
    
    n.2 = input$no[h:k]
    t.2 = adjco.site[n.2]
    t.3 = as.numeric(do.call("c", t.2))
    
    t.4 = data.frame(table(cut(t.3, breaks, right=FALSE)))
    On.date = as.numeric(as.character(t.1$Var1[i]))
    AdjFreq = t.4$Freq/n.1
    t.4 = cbind(t.4, On.date, AdjFreq)
  
  output = rbind(output, t.4)  
  print(i)
}
 
# 4 Heatmap  
  
  
  ggplot(output, aes(Var1, On.date)) + 
    
    geom_tile(aes(fill=Freq), colour="white") +  theme(legend.position="none") +
    
    scale_fill_gradient(low="white", high="Green") + 
    
    theme(axis.text.x=element_text(angle=90, vjust = 0.5)) +
    
    scale_y_continuous(breaks = c(0:8), labels = c(0:8)) + 
    
    xlab("Position") + ylab("Onset Date") + scale_x_discrete(breaks=ticksss, labels=ticklab)

 
# breaks = seq(930, 2430, by=5)                 # length(breaks) = 301
# ticks<-out.p.0.001.1000r.c0$region[1:300]     # originally what showed
# ticksss<-ticks[seq(5,295, by=10)]             # the place we preferred to replace
# ticklab<-seq(950, 2400, by=50) 


  
  
  
    
  
  
  
  
  


  