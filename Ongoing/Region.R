#this .R file focus on region recognition

#lsfullm = lsfull[-(pairedS.No[65])]
#lsfull.C0m = lsfull.C0[-103]
#pairedS.Nom

lsfull.C0m.0102 = lsfull.C0m[1:77] #n = 77
lsfull.C0m.15 = lsfull.C0m[78:150] #n = 73


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

#################################
# outreg.0102.r1000.0.05.ps
# outreg.15.r1000.0.05.ps
# outreg.0102.r1000.0.05.pns
# outreg.15.r1000.0.05.pns
# 
# outreg.0102.r1000.0.01.ps
# outreg.15.r1000.0.01.ps
# outreg.0102.r1000.0.01.pns
# outreg.15.r1000.0.01.pns
#
# outreg.0102.r1000.0.05.vs
# outreg.15.r1000.0.05.vs
# outreg.0102.r1000.0.05.vns
# outreg.15.r1000.0.05.vns


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





  
  
  
  