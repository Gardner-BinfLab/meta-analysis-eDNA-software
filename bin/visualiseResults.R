#!/usr/bin/Rscript
#cd data
#R CMD BATCH --no-save ../bin/visualiseResults.R

library(RColorBrewer)
meths<-read.table("metagenome-meta-analysis-F-measures.tsv", header=T, sep = "\t")
methods<-sort(unique(meths$Method), decreasing = TRUE)
papers <-sort(unique(meths$Paper ), decreasing = FALSE)
cols   <-brewer.pal(length(papers),"Dark2")

######################################################################
###########prescence/absence papers vs methods plot (point-size = median(F)):
#FIGURE 3A

pdf(file="../manuscript/figures/fig3a-matrix-methods-benchmarks.pdf",height=22,width=10)
par(las=2,cex=2,mar=c(5.1+2, 4.1, 4.1+1, 8.1), xpd=TRUE)   
plot(0,xlim=c(1,length(papers)),ylim=c(1,length(methods)),xlab="",ylab="",yaxt='n',xaxt='n',bty='n',pch='')
axis(4, at=1:length(methods), labels=methods  )
axis(1, at=1:length(papers),  labels=papers   )
for(j in length(papers):1){
      lines(c(j,j), c(1,length(methods)), lwd=0.1, col=grey(0.75))
}
for(i in 1:length(methods)){
      lines(c(1,length(papers)), c(i,i), lwd=0.1, col=grey(0.75))
      for(j in length(papers):1){
      	    p  <- papers[j]
	    fp <- meths$F1.measure[methods[i]==meths$Method & p==meths$Paper]
      	    if(length(fp) > 0){
	    	points(j, i, pch=19, col = cols[j], cex=log10(15*(length(fp)-0.9) )+0.5 )                   #cex=2.25*median(fp)+0.75)
		#cat(paste(p, methods[i], length(fp), "\n" ))
		}
      }
}
text(1.05-0.5,length(methods)+2.75, 'Number of F-measure estimates (N)', pos=4)
cnts <- c(2,4,6,8,10)
for(j in 1:length(cnts)){
      points(j-0.5, length(methods)+1.5, pch=19, col=grey(0.5), cex=log10(15*(cnts[j]-0.9))+0.5)
      text(  j+0.05-0.5, length(methods)+1.5, cnts[j], pos=4, cex=0.75)
}
dev.off()
system("convert -flatten ../manuscript/figures/fig3a-matrix-methods-benchmarks.pdf ../manuscript/figures/fig3a-matrix-methods-benchmarks.png")


######################################################################
#bootstrap confidence intervals:
library(bootstrap)
theta <- function(x){median(x)}

######################################################################
#####sort on median normalised F-measure:

fnorm<-matrix(0, nrow = length(meths$Paper), ncol = 1)
for(j in 1:length(papers)){
      p <- papers[j]
      fnorm[p==meths$Paper] <- (meths$F1.measure[p==meths$Paper]-median(meths$F1.measure[p==meths$Paper]) )/mad(meths$F1.measure[p==meths$Paper])
}
medianF <- matrix(0, nrow = length(methods), ncol = 1)
for(i in 1:length(methods)){
      fp <- fnorm[methods[i]==meths$Method]
      medianF[i] <- median(fp)
}
medianF<-sort(medianF, index.return=TRUE)
#write.table(cbind(as.character(methods[medianF$ix]), as.numeric(medianF$x)), file = "medianValsTot.tsv", sep = "\t", append = FALSE, quote = FALSE, col.names =FALSE)  #, names = FALSE
#cat("", file = "medianVals.tsv", sep = "", append = FALSE)

scatterFCI <- matrix(0, nrow = length(methods), ncol = 3)
rownames(scatterFCI) <- methods[medianF$ix]

#FIGURE 3B
pdf(file="../manuscript/figures/fig3b-vizFnorm-dist-sorted.pdf",height=16,width=9)
par(las=1,cex=2,mar=c(5.1, 4.1, 4.1+1, 8.1), xpd=TRUE)
plot(0,xlim=c(-2.5,2.5),ylim=c(1,length(methods)),xlab="Robust Z-score (F-measure)",ylab="",yaxt='n',bty='n',pch='')
axis(4, at=1:length(methods), labels=methods[ medianF$ix ]  )
for(i in 1:length(methods)){
      lines(c(-2.5,2.5), c(i,i), lwd=0.1, col=grey(0.75))
      for(j in length(papers):1){
      	    p <- papers[j]
	    m <- methods[medianF$ix[i]]
	    fp <- fnorm[m==meths$Method & p==meths$Paper]
      	    if(length(fp) > 0){
      	    	points(fp, i*fp/fp, pch=19, col=cols[j])
	    	points(median(fp), i, pch='|', col = cols[j], cex=2, lwd=7)
		#cat(paste(c(as.character(p), as.character(m), median(fp), "\n"), sep="\t"), file = "medianVals.tsv", sep = "\t", append = TRUE)
	    }
      }
      points(medianF$x[i], i, pch='|', col = 'black', cex=1.5, lwd=5,font=2)
      scatterFCI[i,1] <- medianF$x[i]
      f <- fnorm[m==meths$Method ]
      minFnorm <- -1.75
      maxFnorm <-  1.75
      results <- bootstrap(c(rep(minFnorm,5), f, rep(maxFnorm,5)),1000,theta)     
      q       <- quantile(results$thetastar,probs = c(0.025,0.975))
      upper<- q[2]
      lower<- q[1]
      scatterFCI[i,2] <- lower
      scatterFCI[i,3] <- upper
      points(upper, i, pch='|', col = 'black', cex=1, lwd=5,font=2)
      points(lower, i, pch='|', col = 'black', cex=1, lwd=5,font=2)
      lines(c(lower,upper), c(i,i), lwd=0.2, col='black')
}
lines(c(0,0), c(1,length(methods)), lty=2)
for(j in 1:length(papers)){
      p <- papers[j]
      xm <- (j %% 2)
      points(5*xm-3, length(methods)-j*0.5 -xm*0.5+4, pch=19, col=cols[j])
      text(  5*xm-3, length(methods)-j*0.5 -xm*0.5+4, p, pos=4)
}
dev.off()
system("convert -flatten ../manuscript/figures/fig3b-vizFnorm-dist-sorted.pdf ../manuscript/figures/fig3b-vizFnorm-dist-sorted.png")

######################################################################
######################################################################
#Supplementary Figure 4

#####sort on median sensitivity:
medianSens <- matrix(0, nrow = length(methods), ncol = 1)
for(i in 1:length(methods)){
      fp <- meths$Sensitivity[methods[i]==meths$Method]
      medianSens[i] <- median(fp)
}
medianSens<-sort(medianSens, index.return=TRUE)
#write.table(cbind(as.character(methods[medianSens$ix]), as.numeric(medianSens$x)), file = "medianValsTot.tsv", sep = "\t", append = FALSE, quote = FALSE, col.names =FALSE) #names = FALSE, 
cat("", file = "medianValsSens.tsv", sep = "", append = FALSE)
#SUPFIGURE 4A
pdf(file="../manuscript/figures/suppfig4a-vizSens-dist-sorted.pdf",height=22,width=10)
par(las=1,cex=2,mar=c(5.1, 4.1, 4.1+1, 8.1), xpd=TRUE)
plot(0,xlim=c(0,1),ylim=c(1,length(methods)),xlab="Sensitivity",ylab="",yaxt='n',bty='n',pch='')
axis(4, at=1:length(methods), labels=methods[ medianSens$ix ]  )
for(i in 1:length(methods)){
      lines(c(0,1), c(i,i), lwd=0.1, col=grey(0.75))
      for(j in length(papers):1){
      	    p <- papers[j]
	    m <- methods[medianSens$ix[i]]
	    fp <- meths$Sensitivity[m==meths$Method & p==meths$Paper]
      	    if(length(fp) > 0){
      	    	points(fp, i*fp/fp, pch=19, col=cols[j])
	    	points(median(fp), i, pch='|', col = cols[j], cex=2, lwd=7)
		#cat(paste(c(as.character(p), as.character(m), median(fp), "\n"), sep="\t"), file = "medianVals.tsv", sep = "\t", append = TRUE)
	    }
      }
      points(medianSens$x[i], i, pch='|', col = 'black', cex=1.5, lwd=5,font=2)
      f <- meths$Sensitivity[m==meths$Method ]
      results <- bootstrap(c(0,0,0, f, 1,1,1),1000,theta)     
      q       <- quantile(results$thetastar,probs = c(0.025,0.975))
      upper<- q[2]
      lower<- q[1]
      points(upper, i, pch='|', col = 'black', cex=1, lwd=5,font=2)
      points(lower, i, pch='|', col = 'black', cex=1, lwd=5,font=2)
      lines(c(lower,upper), c(i,i), lwd=0.2, col='black')
}
for(j in 1:length(papers)){
      p <- papers[j]
      xm <- (j %% 2)
      points(0.75*xm, length(methods)-j*0.5 -xm*0.5+4, pch=19, col=cols[j])
      text(  0.75*xm, length(methods)-j*0.5 -xm*0.5+4, p, pos=4)
}
dev.off()
system("convert -flatten ../manuscript/figures/suppfig4a-vizSens-dist-sorted.pdf ../manuscript/figures/suppfig4a-vizSens-dist-sorted.png")

######################################################################
#####sort on median PPV:
medianPPV <- matrix(0, nrow = length(methods), ncol = 1)
for(i in 1:length(methods)){
      fp <- meths$PPV[methods[i]==meths$Method]
      medianPPV[i] <- median(fp)
}
medianPPV<-sort(medianPPV, index.return=TRUE)
#write.table(cbind(as.character(methods[medianPPV$ix]), as.numeric(medianPPV$x)), file = "medianValsTot.tsv", sep = "\t", append = FALSE, quote = FALSE, col.names =FALSE) #names = FALSE, 
cat("", file = "medianValsPPV.tsv", sep = "", append = FALSE)
#SUPFIGURE 4B
pdf(file="../manuscript/figures/suppfig4b-vizPPV-dist-sorted.pdf",height=22,width=10)
par(las=1,cex=2,mar=c(5.1, 4.1, 4.1+1, 8.1), xpd=TRUE)
plot(0,xlim=c(0,1),ylim=c(1,length(methods)),xlab="PPV",ylab="",yaxt='n',bty='n',pch='')
axis(4, at=1:length(methods), labels=methods[ medianPPV$ix ]  )
for(i in 1:length(methods)){
      lines(c(0,1), c(i,i), lwd=0.1, col=grey(0.75))
      for(j in length(papers):1){
      	    p <- papers[j]
	    m <- methods[medianPPV$ix[i]]
	    fp <- meths$PPV[m==meths$Method & p==meths$Paper]
      	    if(length(fp) > 0){
      	    	points(fp, i*fp/fp, pch=19, col=cols[j])
	    	points(median(fp), i, pch='|', col = cols[j], cex=2, lwd=7)
		#cat(paste(c(as.character(p), as.character(m), median(fp), "\n"), sep="\t"), file = "medianVals.tsv", sep = "\t", append = TRUE)
	    }
      }
      points(medianPPV$x[i], i, pch='|', col = 'black', cex=1.5, lwd=5,font=2)
      f <- meths$PPV[m==meths$Method ]
      results <- bootstrap(c(0,0,0, f, 1,1,1),1000,theta)     
      q       <- quantile(results$thetastar,probs = c(0.025,0.975))
      upper<- q[2]
      lower<- q[1]
      points(upper, i, pch='|', col = 'black', cex=1, lwd=5,font=2)
      points(lower, i, pch='|', col = 'black', cex=1, lwd=5,font=2)
      lines(c(lower,upper), c(i,i), lwd=0.2, col='black')
}
for(j in 1:length(papers)){
      p <- papers[j]
      xm <- (j %% 2)
      points(0.75*xm, length(methods)-j*0.5 -xm*0.5+4, pch=19, col=cols[j])
      text(  0.75*xm, length(methods)-j*0.5 -xm*0.5+4, p, pos=4)
}
dev.off()
system("convert -flatten ../manuscript/figures/suppfig4b-vizPPV-dist-sorted.pdf ../manuscript/figures/suppfig4b-vizPPV-dist-sorted.png")


######################################################################
#####sort on median F-measure:
medianF <- matrix(0, nrow = length(methods), ncol = 1)
for(i in 1:length(methods)){
      fp <- meths$F1.measure[methods[i]==meths$Method]
      medianF[i] <- median(fp)
}
medianF<-sort(medianF, index.return=TRUE)
write.table(cbind(as.character(methods[medianF$ix]), as.numeric(medianF$x)), file = "medianValsTot.tsv", sep = "\t", append = FALSE, quote = FALSE, col.names =FALSE)  #, names = FALSE
cat("", file = "medianVals.tsv", sep = "", append = FALSE)
#SUPFIGURE 4C
pdf(file="../manuscript/figures/suppfig4c-vizF-dist-sorted.pdf",height=22,width=10)
par(las=1,cex=2,mar=c(5.1, 4.1, 4.1+1, 8.1), xpd=TRUE)
plot(0,xlim=c(0,1),ylim=c(1,length(methods)),xlab="F-measure",ylab="",yaxt='n',bty='n',pch='')
axis(4, at=1:length(methods), labels=methods[ medianF$ix ]  )
for(i in 1:length(methods)){
      lines(c(0,1), c(i,i), lwd=0.1, col=grey(0.75))
      for(j in length(papers):1){
      	    p <- papers[j]
	    m <- methods[medianF$ix[i]]
	    fp <- meths$F1.measure[m==meths$Method & p==meths$Paper]
      	    if(length(fp) > 0){
      	    	points(fp, i*fp/fp, pch=19, col=cols[j])
	    	points(median(fp), i, pch='|', col = cols[j], cex=2, lwd=7)
		cat(paste(c(as.character(p), as.character(m), median(fp), "\n"), sep="\t"), file = "medianVals.tsv", sep = "\t", append = TRUE)
	    }
      }
      points(medianF$x[i], i, pch='|', col = 'black', cex=1.5, lwd=5,font=2)
      f <- meths$F1.measure[m==meths$Method ]
      results <- bootstrap(c(0,0,0, f, 1,1,1),1000,theta)     
      q       <- quantile(results$thetastar,probs = c(0.025,0.975))
      upper<- q[2]
      lower<- q[1]
      points(upper, i, pch='|', col = 'black', cex=1, lwd=5,font=2)
      points(lower, i, pch='|', col = 'black', cex=1, lwd=5,font=2)
      lines(c(lower,upper), c(i,i), lwd=0.2, col='black')
}
for(j in 1:length(papers)){
      p <- papers[j]
      xm <- (j %% 2)
      points(0.75*xm, length(methods)-j*0.5 -xm*0.5+4, pch=19, col=cols[j])
      text(  0.75*xm, length(methods)-j*0.5 -xm*0.5+4, p, pos=4)
}
dev.off()
system("convert -flatten ../manuscript/figures/suppfig4c-vizF-dist-sorted.pdf ../manuscript/figures/suppfig4c-vizF-dist-sorted.png")

######################################################################
#F, Sensitivity and PPV distributions for each paper:

#SUPFIGURE 2A
pdf(file="../manuscript/figures/suppfig2a-hists-Sensitivity-per-paper.pdf",height=21,width=21)
#par(las=1,cex=4,mfrow=c(3,2),mar=c(5.1+0, 4.1, 4.1+1, 2.1))
par(las=1,cex=3.0,mfrow=c(2,2),mar=c(5.1+3, 4.1+3, 4.1+3, 2.1)) #mar: c(bottom, left, top, right); c(5, 4, 4, 2) + 0.1
xrange<-seq(0,1,length.out=11)
for(j in 1:length(papers)){
      p <- papers[j]
      fp<-meths$Sensitivity[p==meths$Paper]
      hv1 <- hist(fp,breaks=xrange,plot=F)$counts
      barplot(hv1,names.arg=format(round(xrange[-1], digits = 2), nsmall = 2),space=0,las=2,col=cols[j],main=paste(p, "Sensitivity", sep="\n"),xlab='', cex.axis = 4, cex.names = 4, cex.main = 4)
}
dev.off()
system("convert -flatten ../manuscript/figures/suppfig2a-hists-Sensitivity-per-paper.pdf ../manuscript/figures/suppfig2a-hists-Sensitivity-per-paper.png")

#SUPFIGURE 2B
pdf(file="../manuscript/figures/suppfig2b-hists-PPV-per-paper.pdf",height=21,width=21)
par(las=1,cex=3.0,mfrow=c(2,2),mar=c(5.1+3, 4.1+3, 4.1+3, 2.1)) #mar: c(bottom, left, top, right); c(5, 4, 4, 2) + 0.1
xrange<-seq(0,1,length.out=11)
for(j in 1:length(papers)){
      p <- papers[j]
      fp<-meths$PPV[p==meths$Paper]
      hv1 <- hist(fp,breaks=xrange,plot=F)$counts
      barplot(hv1,names.arg=format(round(xrange[-1], digits = 2), nsmall = 2),space=0,las=2,col=cols[j],main=paste(p, "PPV", sep="\n"),xlab='', cex.axis = 4, cex.names = 4, cex.main = 4)
}
dev.off()
system("convert -flatten ../manuscript/figures/suppfig2b-hists-PPV-per-paper.pdf ../manuscript/figures/suppfig2b-hists-PPV-per-paper.png")

#SUPFIGURE 3A
pdf(file="../manuscript/figures/suppfig3a-hists-F-per-paper.pdf",height=21,width=21)
par(las=1,cex=3.0,mfrow=c(2,2),mar=c(5.1+3, 4.1+3, 4.1+3, 2.1)) #mar: c(bottom, left, top, right); c(5, 4, 4, 2) + 0.1
xrange<-seq(0,1,length.out=11)
for(j in 1:length(papers)){
      p <- papers[j]
      fp<-meths$F1.measure[p==meths$Paper]
      hv1 <- hist(fp,breaks=xrange,plot=F)$counts
      barplot(hv1,names.arg=format(round(xrange[-1], digits = 2), nsmall = 2),space=0,las=2,col=cols[j],main=paste(p, "F-measure", sep="\n"),xlab='', cex.axis = 4, cex.names = 4, cex.main = 4)
}
dev.off()
system("convert -flatten ../manuscript/figures/suppfig3a-hists-F-per-paper.pdf ../manuscript/figures/suppfig3a-hists-F-per-paper.png")

#SUPFIGURE 3B
pdf(file="../manuscript/figures/suppfig3b-hists-F-norm-per-paper.pdf",height=21,width=21)
par(las=1,cex=3.0,mfrow=c(2,2),mar=c(5.1+3, 4.1+3, 4.1+3, 2.1)) #mar:c(bottom, left, top, right); c(5, 4, 4, 2) + 0.1
xrange<-seq(-2.5,2.5,length.out=11)
#CAN THROW AN ERROR, AUTOMATE FINDING THE LIMITS????
for(j in 1:length(papers)){
      p <- papers[j]
      fp<-(meths$F1.measure[p==meths$Paper]-median(meths$F1.measure[p==meths$Paper]))/mad(meths$F1.measure[p==meths$Paper])
      cat(paste(p, min(fp), max(fp),"\n" ))
      hv1 <- hist(fp,breaks=xrange,plot=F)$counts
      barplot(hv1,names.arg=format(round(xrange[-1], digits = 1), nsmall =1),space=0,las=2,col=cols[j],main=paste(p, "Robust Z-scores (F)",sep="\n"),xlab='', cex.axis = 4, cex.names = 4, cex.main = 4)
}
dev.off()
system("convert -flatten ../manuscript/figures/suppfig3b-hists-F-norm-per-paper.pdf ../manuscript/figures/suppfig3b-hists-F-norm-per-paper.png")


#############################

#FIGURE 3C
pdf(file="../manuscript/figures/fig2c-boxplots-F-per-paper.pdf",height=8,width=8)
par(las=2,cex=2,mar=c(5.1+2, 4.1, 4.1-2, 1.1))  #c(bottom, left, top, right) c(5, 4, 4, 2) + 0.1
boxplot(meths$F1.measure ~ meths$Paper, ylab="F-measure")
for(j in 1:length(papers)){
      p <- papers[j]
      fp<-meths$F1.measure[p==meths$Paper]
      x<-rnorm(length(fp), mean = j, sd = 0.15)
      x[x>j+0.35]<-j+0.35
      x[x<j-0.35]<-j-0.35
      points(x,fp, col="black", bg=cols[j], pch=21, cex=0.35)
}
dev.off()

######################################################################
#Citations
cites<-read.table("citations.tsv", header=T, sep = "\t")
#HOOP JUMPING TO COMBINE PAPER COLOURS:
citeCols = matrix('black', length(cites$Method))
tfCitesMethods <- cites$Method %in%  as.character(methods)
#true if method has been benchmarked!
for(i in 1:length(cites$Method)){
      if(tfCitesMethods[i]){
	citePapers <- as.character(sort(unique(meths$Paper[ meths$Method %in% cites$Method[i] ])))
	cat(citePapers, "\n")
	cp <- cols[papers %in% citePapers]
	citeCols[i] <- cp[1]
	if(length(citePapers)>1){
	   for(j in 2:length(citePapers)){
	      citeCols[i] <- colorRampPalette( c(citeCols[i], cp[j]), space = "Lab")(3)[2]
	   }
	}
      }
}

#FIGURE 2A-B
pdf(file="../manuscript/figures/fig2ab-citations2-year.pdf",height=17,width=10)
par(las=3,cex=3,mfrow=c(2,1), mar=c(5, 4+2, 4, 2))
xrange <- 2007:2019-0.5
hv1 <- hist(cites$Year[cites$Evaluated == 'y'], breaks=xrange,plot=F)
hv2 <- hist(cites$Year[cites$Evaluated == 'n'], breaks=xrange,plot=F)
cnts<-rbind(hv1$counts,hv2$counts)
colnames(cnts) <- hv1$mids
rownames(cnts) <- c('Evaluated','Unevaluated')
barplot(cnts,space=0,las=3,xlab='', ylab='Number of eDNA analysis methods', main='', col=c('lightskyblue','peachpuff'), cex.main = 2, cex.lab = 2, cex.axis = 2, cex.names = 2) 
legend("topleft", c('Evaluated','Unevaluated'), col=c('lightskyblue','peachpuff'), fil=c('lightskyblue','peachpuff'), cex=2)
plot(0, ylab="Number of citations",yaxt='n', xlab="",yaxt='n',xaxt='n', xlim=c(2007, 2018.5), ylim=c(0,4),bty='n'  , cex.axis = 2, cex.main = 2, cex.lab = 2)  
par(las=3)
axis(1, at=2007:2018, labels=c(2007,'',2009,'',2011,'',2013,'',2015,'',2017,'')  , cex.axis = 2, cex.main = 2) 
axis(2, at=0:4, labels=10^(0:4)                                               , cex.axis = 2, cex.main = 2) 
for(i in 1:length(cites$Method)){
      if(cites$Evaluated[i] == 'y'){
	points(cites$Year[i], log10(cites$NumberCitations.GS.2017.07.07[i]), col=citeCols[i], pch=19, cex=1.5)
      	text(cites$Year[i]-0.0, log10(cites$NumberCitations.GS.2017.07.07[i]), cites$Method[i], cex=1.5, col=citeCols[i], pos=4, font=2)
      }
      else{
	points(cites$Year[i], log10(cites$NumberCitations.GS.2017.07.07[i]), col='black', pch=20, cex=1.2)
	if(cites$NumberCitations.GS.2017.07.07[i]>100){
	     text(cites$Year[i]-0.0, log10(cites$NumberCitations.GS.2017.07.07[i]), cites$Method[i], cex=1.0, col='black', pos=4)
	}
      }
}
dev.off()

#FIGURE 2D
#Citations vs median Robust Z:
xx<-matrix(NA, nrow = length(scatterFCI[,1]), ncol = 1)
yy<-matrix(NA, nrow = length(scatterFCI[,1]), ncol = 1)
nn<-matrix(NA, nrow = length(scatterFCI[,1]), ncol = 1)
j<-1
for(i in 1:length(cites$Method)){

      if(cites$Evaluated[i] == 'y' && length(scatterFCI[as.character(cites$Method[i])==rownames(scatterFCI),1])){
	    xx[j] <-  scatterFCI[as.character(cites$Method[i])==rownames(scatterFCI),1]
            yy[j] <- log10( cites$NumberCitations.GS.2017.07.07[i]/(2019 - cites$Year[i]) )
	    nn[j] <- as.character(cites$Method[i])
	    j<-j+1
	}
}

pdf(file="../manuscript/figures/fig2d-citations2-year-vs-zF.pdf",height=8,width=8)
par(cex=1.5,mar=c(5.1, 4.1+2, 4.1, 2.1)) #mar: c(bottom, left, top, right); c(5, 4, 4, 2) + 0.1
plot(xx,yy,pch=19, col="blue", ylab="Citations/Year",yaxt='n',yaxt='n', xlab="Robust Z-score (F-measure)",xlim=c(-2, 1), ylim=c(0,3),bty='n'  , cex.axis = 2, cex.main = 2, cex.lab = 2)
axis(2, at=0:3, labels=10^(0:3)                                               , cex.axis = 2, cex.main = 2) 
abline( lm(yy ~ xx), col="red", lwd=5)
fitted <- predict( lm(yy ~ xx) )
for( i in 1:length(xx)){
    lines( c(xx[i], xx[i]), c(yy[i], fitted[i]), lwd=3, col="purple4")
    text(xx[i], yy[i], nn[i], pos=4, cex=0.75)
}
ct<-cor.test(xx,yy)
text(-2,3.00, paste("Pearson's correlation:", signif(as.numeric(ct$estimate), digits = 2), sep=""), pos=4,cex=0.75)
text(-2,2.85, paste("P-value:", signif(ct$p.value, digits = 2), sep=" "), pos=4,cex=0.75)
dev.off()

######################################################################
######################################################################
#NETWORK DIAGRAM:

library(igraph)
system("../bin/makeGephiFiles.sh")

edges<-read.csv("medianVals-gephiEdges.tsv", header=T)
nodes<-read.table("medianValsTot.tsv", header=F, row.names=1) #medianVals-gephiNodes.tsv
net <- graph_from_data_frame(d=edges, vertices=nodes, directed=T) 

#FIGURE 4B  -- uses a random seed, so different each time
pdf(file="../manuscript/figures/fig4b-networkDAG-accuracy-ranks.pdf",height=10,width=10)
par(mar=c(.1, .1, 5, .1))
plot(net, layout=layout_with_lgl, 
edge.arrow.size=1.0/edges$weight, edge.width=5/edges$weight, edge.curved=0.2, edge.color=cols[as.factor(edges$type)],
vertex.color="#FB9A99", vertex.frame.color="lightgrey", vertex.size=10*nodes$V3+1,
vertex.label.color="black", vertex.label.cex=1.5, vertex.label.dist=0.0)
legend("topright", as.character(papers), col=cols, fil=cols, cex=1.2)
#
text(0,1.2, 'Median F-measure')
for(j in seq(0,1.0,length.out=6)){
      points(j-0.5,     1.25, pch=21, col="lightgrey", bg="#FB9A99", cex=4.5*j+0.65)
      text(  j-0.5+0.02, 1.25, j, pos=4)
}
dev.off()
system("convert -flatten ../manuscript/figures/fig4b-networkDAG-accuracy-ranks.pdf ../manuscript/figures/fig4b-networkDAG-accuracy-ranks.png")

######################################################################
######################################################################
#Network meta-analysis

library(netmeta)
#Example:
#Senn2013

#INITIALIZE:
#treatment, n, mean, sd, studlab
studlab.1<-character() #as.character(papers[1])
methlab.1<-character() #as.character(methods[1])
te.1<-numeric() #mc$TE
sete.1<-numeric() #mc$seTE
ctrl<-character() #'Control'

#rbind(data frameA, data frameB)
for(i in length(papers):1){
      for(j in 1:(length(methods)-1)){
      	    p <- papers[i]
	    m <- methods[j]
	    fp <- c(0,meths$F1.measure[m==meths$Method & p==meths$Paper], 1)
      	    if(length(fp) > 2){
	        fp2<-fp
		#for(k in (j+1):(length(methods)-1)){
		for(k in (j):(length(methods))){
		      m3 <- methods[k]
		      print(paste(m, " vs ", m3))
	    	      fp3 <- c(0,meths$F1.measure[m3==meths$Method & p==meths$Paper], 1)
		      if(length(fp3) > 2 & length(m)>0 & length(m3)>0 & m3 != m){
		      		  fp4<-fp3
				  mc2<-metacont(length(fp3), mean(fp3), sd(fp3), length(fp2), mean(fp2), sd(fp2), studlab=p,
				                comb.random=TRUE, method.bias='rank', label.e=as.character(m3), label.c=as.character(m))
				  studlab.1<-rbind(studlab.1, as.character(p))
				  methlab.1<-rbind(methlab.1, as.character(m3))
				  te.1     <-rbind(te.1,    mc2$TE)
				  sete.1   <-rbind(sete.1,  mc2$seTE)
				  ctrl     <-rbind(ctrl,    as.character(m))  
		      }
		}
	    }
      }
}

df.mc<-data.frame(TE=te.1, seTE=round(sete.1, 4), studlab=studlab.1, treat1=methlab.1, treat2=ctrl, row.names=1:length(te.1))
head(df.mc)
df.mc

#Function ‘pairwise’ automatically calculates all pairwise comparisons for multi-arm studies.
#netconnection(treat1, treat2, studlab, data = df.mc)

mn1 <- netmeta(TE, seTE, treat1, treat2, studlab, data=df.mc, sm="OR", comb.random=TRUE)
#"RD"’, ‘"RR"’, ‘"OR"’, ‘"ASD"’, ‘"HR"’, ‘"MD"’, ‘"SMD"’, or ‘"ROM"
netgraph(mn1)

#FIGURE 4A 
pdf(file="../manuscript/figures/fig4a-forest-plot-F-measure.pdf",height=10,width=5) #font
fr.out<-forest.netmeta(mn1, sortvar=-TE, layout="RevMan5", ref="EBI-mg", prediction=TRUE)
dev.off()

metabias(df.mc$TE, seTE=df.mc$seTE, method.bias='rank', plotit=TRUE, correct=TRUE)

#download.packages(pkgs = "netmeta", destdir = ".",type = "source")
#edit forest.netmeta function to extract ORs & CIs 
print(fr.out)

#cat forest-vals.tsv | tr -d '[;]()%\42' | egrep -v 'Wrandom|EBI-mg|V1' | perl -lane '$a=join("\t", @F[1..$#F]); print $a' > forest-vals2.tsv
fr.out<-read.table("forest-vals2.tsv", header=F, sep = "\t")

#Supplementary Figure 5
#Effect size vs confidence width plots
pdf(file="../manuscript/figures/suppfig5-effect-vs-confidence.pdf",height=5,width=10) #font
par(las=1,cex=1.5,mfrow=c(1,2),mar=c(5.1+0, 4.1, 4.1+1, 2.1))
plot(scatterFCI[,1], abs(scatterFCI[,3]-scatterFCI[,2]), xlab="Robust Z-score (F-measure)", ylab="Width of 95% confidence interval", col='black', pch=19, ylim=c(0,4.5))
for(i in 1:length(methods)){
      text(scatterFCI[i,1], abs(scatterFCI[i,3]-scatterFCI[i,2]),rownames(scatterFCI)[i], pos=3, cex=0.5)
}
lines(c(0,1,1.25),c(0,1,1.25), lwd=0.1, col=grey(0.75))
#
plot(fr.out$V2, abs(fr.out$V4-fr.out$V3), xlab="Odds-ratio", ylab="Width of 95% confidence interval", col='black', pch=19, xlim=c(1,2), ylim=c(1,2.7))
for(i in 1:length(fr.out$V2)){
      text(fr.out$V2[i], abs(fr.out$V4[i]-fr.out$V3[i]),fr.out$V1[i], pos=3, cex=0.5)
}
dev.off()
#Methods with CI that excludes 1:
fr.out[fr.out$V3 > 1, ]

#removing EBI-mg (ground state for network meta-analysis: 
ixFCI<-sort(row.names(scatterFCI)[ row.names(scatterFCI) != "EBI-mg"], index.return = T)$ix
ixFR<-sort(as.vector(fr.out$V1), index.return = T)$ix

x<-scatterFCI[ixFCI+1,1]
y<-matrix(fr.out$V2[ixFR], nrow = 1, ncol = length(ixFR))
colnames(y) <- fr.out$V1[ixFR]

pdf(file="../manuscript/figures/suppfig6-robustZ-vs-oddsratios.pdf",height=5,width=5) 
par(cex=1.0)
plot(x,y,pch=19, col="blue",xlab="Robust Z-score", 
                 ylab = "Odds ratio",  
                 main = "Comparison of robust Z\n & network meta-analysis values")
abline( lm(y[1,] ~ x), col="red", lwd=5)
fitted <- predict( lm(y[1,] ~ x) )
for( i in 1:length(x)){
    lines( c(x[i], x[i]), c(y[1,i], fitted[i]), lwd=3, col="purple4")
    text(x[i], y[1,i], names(x)[i], pos=4, cex=0.3)
}
ct<-cor.test(x,y[1,])
text(0.5,1.2, paste("Pearson's correlation:", signif(as.numeric(ct$estimate), digits = 2), sep=""), pos=2,cex=0.5)
text(0.5,1.175, paste("P-value:", signif(ct$p.value, digits = 2), sep=" "), pos=2,cex=0.5)
dev.off()

######################################################################

# Bazinet.2012 	  MEGAN (0.80)      > MetaPhyler (0.01) 
# Lindgreen.2016 	  CLARK (0.98)      > Kraken (0.95) > MEGAN (0.70) > MetaPhyler (0.01)
# McIntyre.2017 	  CLARK (0.93)      > Kraken (0.90) > MEGAN (0.87)
# Peabody.2015	  MetaPhyler (0.44) > CLARK (0.28) > Kraken (0.28)
# Sczyrba.2017	  MetaPhyler (0.33) > Kraken (0.31) > CLARK (0.20) > MEGAN (0.06)

# Bazinet.2012 	  MEGAN (0.80)      > MetaPhyler (0.01) 
# Lindgreen.2016 	  CLARK (0.98)      > MEGAN (0.70) > MetaPhyler (0.01)
# McIntyre.2017 	  CLARK (0.93)      > MEGAN (0.87)

# Peabody.2015	  MetaPhyler (0.44) > CLARK (0.28)
# Sczyrba.2017	  MetaPhyler (0.33) > CLARK (0.20) > MEGAN (0.06)

# Versions:       CLARK      MEGAN      MetaPhyler
# Bazinet.2012 	  -          4.61.5     1.13
# Lindgreen.2016  1.1.3      5.7.0      1.25
# McIntyre.2017   1.2.2-beta 5.10.6     -
# Peabody.2015	  1.1.3      4.70.4     1.25
# Sczyrba.2017	  1.1.3      6.4.9      1.25



















