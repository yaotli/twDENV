#this .R file focus on region recognition

#lsfullm = lsfull[-(pairedS.No[65])]
#lsfull.C0m = lsfull.C0[-103]
#pairedS.Nom

lsfull.C0m.0102 = lsfull.C0m[1:77] #n = 77
lsfull.C0m.15 = lsfull.C0m[78:150] #n = 73


# use pvariant -----
breaks = seq(930, 2430, by=5) # length(breaks) = 301
ticks<-outreg.15.r1000.0.05.ps$region[1:300]
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



