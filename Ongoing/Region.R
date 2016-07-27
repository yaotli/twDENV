#this .R file focus on region recognition

# E protein: (nt) 937 - 2421 
# E protein: (aa) 313 - 807

# lsfullm = lsfull[-(pairedS.No[65])]
# lsfull.C0m = lsfull.C0[-103]
# pairedS.Nom

# lsfull.C0m.0102 = lsfull.C0m[1:77] #n = 77
# lsfull.C0m.15 = lsfull.C0m[78:150] #n = 73


# use pvariant -----
breaks = seq(930, 2430, by=5) # length(breaks) = 301
ticks<-outreg.df.r1000.0.02.vs$region[1:300]
ticksss<-ticks[seq(5,295, by=10)] #length(ticksss)
ticklab<-seq(950, 2400, by=50)

outreg.0102.r1000.0.05.vs = pvariant(lsfull.C0m.0102, breaks, 0.05, 3, 2)
outreg.0102.r1000.0.05.vns = pvariant(lsfull.C0m.0102, breaks, 0.05, 3, 4)
outreg.15.r1000.0.05.vs = pvariant(lsfull.C0m.15, breaks, 0.05, 3, 2)
outreg.15.r1000.0.05.vns = pvariant(lsfull.C0m.15, breaks, 0.05, 3, 4)

#heatmap
library(ggplot2)


hp = ggplot(outreg.15.r1000.0.05.vns, aes(region, k)) + geom_tile(aes(fill=Freq), colour="white") + 
            
            scale_fill_gradient(low="white", high="Blue") + 
  
            theme(axis.text.x=element_text(angle=90, vjust = 0.5)) + xlab("") + ylab("Sample") +
  
            scale_x_discrete(breaks=ticksss, labels=ticklab)


#########################################################################################
################ DHF ####################################################################

# t.0.025 from cdvariant

library(dplyr)  

no.DF.c0 = t.demo %>%
  filter(DHF == "DF") %>%
  select(c.1.150.)

ls.DF.c0 = lsfull.C0m[no.DF.c0[,1]] # n = 127
ls.DHF.c0 = lsfull.C0m[no.DHF.c0[,1]] # n = 22

outreg.df.r1000.0.02.vs = pvariant(ls.DF.c0, breaks, 0.02, 3, 1)
outreg.dhf.r1000.0.02.vs = pvariant(ls.DHF.c0, breaks, 0.02, 3, 1)

# outreg.df.r1000.0.02.vs
# outreg.dhf.r1000.0.02.vs

#heatmap

library(ggplot2)
hp = ggplot(outreg.df.r1000.0.02.vs, aes(region, k)) + geom_tile(aes(fill=Freq), colour="white") + 
  
  scale_fill_gradient(low="white", high="Red") + 
  
  theme(axis.text.x=element_text(angle=90, vjust = 0.5)) + xlab("") + ylab("Sample") +
  
  scale_x_discrete(breaks=ticksss, labels=ticklab)


multiplot(hp,hp2)


##########################
###distribution : variants 
# ouls.df.r1000.0.02.vs
# ouls.dhf.r1000.0.02.vs

ouls.df.r1000.0.02.vs = lsvariantNS(ls.DF.c0, 0.02,3)
ouls.dhf.r1000.0.02.vs = lsvariantNS(ls.DHF.c0, 0.02,3)

##hist

a=ouls.df.r1000.0.02.vs[,1]
a=data.frame(a)
b=ouls.dhf.r1000.0.02.vs[,1]
b=data.frame(b)

  p1=ggplot(a,aes(a)) + geom_histogram(binwidth = 5) + xlim(c(950,2500)) + xlab("")
  p2=ggplot(b,aes(b)) + geom_histogram(binwidth = 5) + xlim(c(950,2500)) + xlab("")
  
multiplot(p1,p2)
  
#normalized 

aa=as.data.frame(table(a)/127)
bb=as.data.frame(table(b)/22) 
names(bb) = c("a", "Freq")

df = rbind(aa,bb)
p2 = as.numeric(as.character(df$a))

df = data.frame(df, p2, c(rep("DF",772), rep("DHF", 453)))
names(df) = c("Position", "Freq", "Position2","DHF")

  ggplot(df, aes(x=Position2, y=Freq)) + geom_line(aes(color=DHF), size=1.5) +
    xlab("Position") + ylab("Normalized count") + 
     scale_color_manual(values=c("#56B4E9", "#CC79A7"))  
    
  
###########################
###distribution : position

# outls.df.r1000.0.02.ps
# outls.dhf.r1000.0.02.ps
    
outls.df.r1000.0.02.ps = pvariantM(ls.DF.c0, 0.02, 3, 1)
outls.dhf.r1000.0.02.ps = pvariantM(ls.DHF.c0, 0.02, 3, 1)

#his

a=outls.df.r1000.0.02.ps[,1]
a=data.frame(a)
b=outls.dhf.r1000.0.02.ps[,1]
b=data.frame(b)

p1=ggplot(a,aes(a)) + geom_histogram(binwidth = 5) + xlim(c(950,2500)) + xlab("")
p2=ggplot(b,aes(b)) + geom_histogram(binwidth = 5) + xlim(c(950,2500)) + xlab("")

multiplot(p1,p2)


#normalized 

aa=as.data.frame(table(a)/127)
bb=as.data.frame(table(b)/22) 
names(bb) = c("a", "Freq")

df = rbind(aa,bb)
p2 = as.numeric(as.character(df$a))

df = data.frame(df, p2, c(rep("DF",772), rep("DHF", 453)))
names(df) = c("Position", "Freq", "Position2","DHF")


ggplot(df, aes(x=Position2, y=Freq)) + geom_line(aes(color=DHF), size=1.5) +
  xlab("Position") + ylab("Normalized count") + 
  scale_color_manual(values=c("#56B4E9", "#CC79A7"))  


################################################################################
### AA position #################################################### DF vs DHF #


# lsfull.aa <- c(list.files(getwd())[213:221], list.files(getwd())[1:212])
# lsfull.c.aa <- strsplit(lsfull, split="-", fixed = T) #cut by "-"
# cstatus.aa <- sapply(lsfull.c.aa, function(x) length(x) == 1 )
# cstatus.aa <- gsub("TRUE", replacement = "C0", cstatus.aa)
# cstatus.aa <- gsub("FALSE", replacement = "C1", cstatus.aa)         #all C status

# lsfull.C0.aa <- lsfull.aa[which(cstatus.aa == "C0")]    # n = 151
# lsfull.C1.aa <- lsfull.aa[which(cstatus.aa == "C1")]    # n = 70
# lsfull.C0m.aa = lsfull.C0.aa[-103]                      # eliminate DENV5083

# ls.DF.c0.aa = lsfull.C0m.aa[no.DF.c0[,1]]               # n = 127
# ls.DHF.c0.aa = lsfull.C0m.aa[no.DHF.c0[,1]]             # n = 22; eliminate DV2355 


##### Use pvariantAA.E #####

#position
 outreg.df.r1000.0.02.p = pvariantAA.E(ls.DF.c0.aa, 0.02, 3, 1)
outreg.dhf.r1000.0.02.p = pvariantAA.E(ls.DHF.c0.aa, 0.02, 3, 1)
 outreg.df.r1000.0.05.p = pvariantAA.E(ls.DF.c0.aa, 0.05, 3, 1)
outreg.dhf.r1000.0.05.p = pvariantAA.E(ls.DHF.c0.aa, 0.05, 3, 1)

#variant
 outreg.df.r1000.0.02.v = pvariantAA.E(ls.DF.c0.aa, 0.02, 3, 2)
outreg.dhf.r1000.0.02.v = pvariantAA.E(ls.DHF.c0.aa, 0.02, 3, 2)
 outreg.df.r1000.0.05.v = pvariantAA.E(ls.DF.c0.aa, 0.05, 3, 2)
outreg.dhf.r1000.0.05.v = pvariantAA.E(ls.DHF.c0.aa, 0.05, 3, 2)


test = rbind(outreg.df.r1000.0.02.p, outreg.dhf.r1000.0.02.p)

library(ggplot2)

ticks.aa<-outreg.df.r1000.0.02.v$position[1:495]
tt=seq(0, 495, by=5)[-1]
ticklab.aa<-seq(0, 495, by=5)[-1]

# from Region.R

test = data.frame(0, outreg.df.r1000.0.05.v$position[1:495], g.E.poly-1, 0)
names(test) = c("k","position", "variantaa", "tcc")   
test = rbind(outreg.df.r1000.0.08.v, test)


vh = ggplot(test, aes(position, k)) + geom_tile(aes(fill=variantaa), colour="white") + 
  
  scale_fill_gradient(low="white", high="Blue") + theme_bw()  +
  
  theme(axis.text.x=element_text(angle=90, vjust = 0.5)) + xlab("") + ylab("Sample") +
  
  scale_x_continuous(breaks=tt, labels=tt)

ch = ggplot(outreg.dhf.r1000.0.08.v, aes(position, k)) + geom_tile(aes(fill=variantaa), colour="white") + 
  
  scale_fill_gradient(low="white", high="Red") + theme_bw()  +
  
  theme(axis.text.x=element_text(angle=90, vjust = 0.5)) + xlab("") + ylab("Sample") +
  
  scale_x_continuous(breaks=tt, labels=tt)


multiplot(vh, ch)







#########################################################################################
################ Repeated Sapmle ########################################################


library(dplyr)

no.repeated = 
t.demo %>% 
  filter(Repeated != ".") %>% 
    select(c.1.150.)

ls.rp.c0 = lsfull.C0m[no.repeated[,1]] # n = 58 ! should remove DV2132 k=5

outls.rp.r1000.0.016.ps = pvariantM(ls.rp.c0, 0.016, 3, 1)
outls.rp.r1000.0.02.ps = pvariantM(ls.rp.c0, 0.02, 3, 1)

# modify to remove DV2132 and add ID and T

outls.rp.r1000.0.016.psM = 
outls.rp.r1000.0.016.ps %>%
  filter(mxk != 5) %>%
    select(mx, mxk)

#call Variability.R:: t.0.016r

cc = t.0.016r$Repeated
t = c(duplicated(t.0.016r$Repeated))
t = recode(as.character(t), "FALSE" = "T1","TRUE" = "T2")
tt = c(17, 33, 36, 39, 44, 45, 48, 51, 54, 57)
t[tt] = "T3"
t[52] = "T4"

d0 = sapply(outls.rp.r1000.0.016.psM$mxk, function(x){if(x >4){ x -1} else{x} })
d1 = sapply(d0, function(x) cc[x])
d2 = sapply(d0, function(x) t[x])
pp = as.numeric(as.character(outls.rp.r1000.0.016.psM$mx))
outls.rp.r1000.0.016.psM$mx

outls.rp.r1000.0.016.psM = data.frame(outls.rp.r1000.0.016.psM, d1, d2, pp)

psm=ggplot(outls.rp.r1000.0.016.psM,aes(mx, color=d2)) + 
  geom_histogram(binwidth = 5, fill="white") + facet_wrap(~ d1) +  theme_bw() +
  scale_color_manual(values=c("Green", "Blue", "Red", "Purple")) + xlab("Position") +
  ylab("Cumulative counts")
