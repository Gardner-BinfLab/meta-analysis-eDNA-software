#!/usr/bin/Rscript
#R CMD BATCH --no-save vizF-dist.R

#TODO:
#1. split taxonomic binners and profilers.


library(RColorBrewer)
meths<-read.table("metagenome-meta-analysis-F-measures2.tsv", header=T, sep = "\t")
methods<-sort(unique(meths$Method), decreasing = TRUE)
papers <-sort(unique(meths$Paper ), decreasing = FALSE)
cols   <-brewer.pal(length(papers),"Dark2")


######################################################################
###########prescence/absence papers vs methods plot (point-size = median(F)):

pdf(file="vizF-matrix-F.pdf",height=22,width=10)
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
	    	points(j, i, pch=19, col = cols[j], cex=log10(length(fp))+0.75 )                   #cex=2.25*median(fp)+0.75)
	    }
      }
}
text(2.5+1,length(methods)+4, 'Number of F-measure estimates (N)')
cnts <- c(1,10,50,100,150,300)
for(j in 1:length(cnts)){
      points(j, length(methods)+2, pch=19, col=grey(0.5), cex=log10(cnts[j])+0.75)
      text(  j+0.05, length(methods)+2, cnts[j], pos=4, cex=0.75)
      #points(j*5+1, length(methods)+2, pch=19, col=grey(0.5), cex=2.25*j+0.75)
      #text(  j*5+1.05, length(methods)+2, j, pos=4)
}
dev.off()


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
pdf(file="vizFnorm-dist-sorted.pdf",height=22,width=10)
par(las=1,cex=2,mar=c(5.1, 4.1, 4.1+1, 8.1), xpd=TRUE)
plot(0,xlim=c(-5.5,3.5),ylim=c(1,length(methods)),xlab="Robust Z-score (F-measure)",ylab="",yaxt='n',bty='n',pch='')
axis(4, at=1:length(methods), labels=methods[ medianF$ix ]  )
for(i in 1:length(methods)){
      lines(c(-5.5,3.5), c(i,i), lwd=0.1, col=grey(0.75))
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
      f <- fnorm[m==meths$Method ]
      results <- bootstrap(c(-5,-5,-5, f, 3,3,3),1000,theta)     
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
      points(6*xm-4, length(methods)-j*0.5 -xm*0.5+6, pch=19, col=cols[j])
      text(  6*xm-4, length(methods)-j*0.5 -xm*0.5+6, p, pos=4)
}
dev.off()



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
pdf(file="vizF-dist-sorted.pdf",height=22,width=10)
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
      points(0.75*xm, length(methods)-j*0.5 -xm*0.5+6, pch=19, col=cols[j])
      text(  0.75*xm, length(methods)-j*0.5 -xm*0.5+6, p, pos=4)
}
dev.off()



######################################################################
#####sort on median sensitivity:
medianSens <- matrix(0, nrow = length(methods), ncol = 1)
for(i in 1:length(methods)){
      fp <- meths$Sensitivity[methods[i]==meths$Method]
      medianSens[i] <- median(fp)
}
medianSens<-sort(medianSens, index.return=TRUE)
#write.table(cbind(as.character(methods[medianSens$ix]), as.numeric(medianSens$x)), file = "medianValsTot.tsv", sep = "\t", append = FALSE, quote = FALSE, col.names =FALSE) #names = FALSE, 
cat("", file = "medianValsSens.tsv", sep = "", append = FALSE)
pdf(file="vizSens-dist-sorted.pdf",height=22,width=10)
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
      points(0.75*xm, length(methods)-j*0.5 -xm*0.5+6, pch=19, col=cols[j])
      text(  0.75*xm, length(methods)-j*0.5 -xm*0.5+6, p, pos=4)
}
dev.off()

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
pdf(file="vizPPV-dist-sorted.pdf",height=22,width=10)
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
      points(0.75*xm, length(methods)-j*0.5 -xm*0.5+6, pch=19, col=cols[j])
      text(  0.75*xm, length(methods)-j*0.5 -xm*0.5+6, p, pos=4)
}
dev.off()


######################################################################
#NETWORK DIAGRAM:


library(igraph)
library(extrafont)
#loadfonts() #run once

#echo "from,to,type,weight" > medianVals-gephiEdges.tsv
#cat medianVals.tsv | sort -k1,1d -k3,3nr | perl -lane 'if(defined($prevM) && ($prevP eq $F[0]) ){$cnt++; print "$F[1],$prevM,$F[0],$cnt"}else{$cnt=0;} ($prevP, $prevM)=($F[0], $F[1]); ' >> medianVals-gephiEdges.tsv

system("makeGephiFiles.sh")

edges<-read.csv("medianVals-gephiEdges.tsv", header=T)
nodes<-read.table("medianValsTot.tsv", header=F, row.names=1) #medianVals-gephiNodes.tsv
net <- graph_from_data_frame(d=edges, vertices=nodes, directed=T) 

pdf(file="networkDAG-accuracy-ranks.pdf",family="FreeSans",height=10,width=10)
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

######################################################################
#F-score distributions for Kraken, Clark & MetaPhyler: 
pdf(file="kraken-clark-metaphyler-hists.pdf",height=17,width=10)
par(las=1,cex=4,mfrow=c(4,1),mar=c(5.1+0, 4.1, 4.1+1, 2.1))
xrange<-seq(0,1,length.out=21)
remove(cnts)
for(j in 1:length(papers)){
      p <- papers[j]
      fp<-meths$F1.measure['Kraken'==meths$Method & p==meths$Paper]
      hv1 <- hist(fp,breaks=xrange,plot=F)$counts
      if(exists("cnts")){
	cnts<-rbind(cnts,hv1)
      }
      else{
	cnts<-hv1
      }
}
barplot(cnts,names.arg=format(round(xrange[-1], digits = 2), nsmall = 2),space=0,las=2,col=cols,main='Kraken',xlab='', cex.axis = 2, cex.names = 2, cex.main = 2)
remove(cnts)
for(j in 1:length(papers)){
      p <- papers[j]
      fp<-meths$F1.measure['CLARK'==meths$Method & p==meths$Paper]
      hv1 <- hist(fp,breaks=xrange,plot=F)$counts
      if(exists("cnts")){
	cnts<-rbind(cnts,hv1)
      }
      else{
	cnts<-hv1
      }
}
barplot(cnts,names.arg=format(round(xrange[-1], digits = 2), nsmall = 2),space=0,las=2,col=cols,main='CLARK',xlab='', cex.axis = 2, cex.names = 2, cex.main = 2)
remove(cnts)
for(j in 1:length(papers)){
      p <- papers[j]
      fp<-meths$F1.measure['MEGAN'==meths$Method & p==meths$Paper]
      hv1 <- hist(fp,breaks=xrange,plot=F)$counts
      if(exists("cnts")){
	cnts<-rbind(cnts,hv1)
      }
      else{
	cnts<-hv1
      }
}
barplot(cnts,names.arg=format(round(xrange[-1], digits = 2), nsmall = 2),space=0,las=2,col=cols,main='MEGAN',xlab='', cex.axis = 2, cex.names = 2, cex.main = 2)
remove(cnts)
for(j in 1:length(papers)){
      p <- papers[j]
      fp<-meths$F1.measure['MetaPhyler'==meths$Method & p==meths$Paper]
      hv1 <- hist(fp,breaks=xrange,plot=F)$counts
      if(exists("cnts")){
	cnts<-rbind(cnts,hv1)
      }
      else{
	cnts<-hv1
      }
}
par(mar=c(5.1+2, 4.1, 4.1+1, 2.1))
barplot(cnts,names.arg=format(round(xrange[-1], digits = 2), nsmall = 2),space=0,las=2,col=cols,main='MetaPhyler',xlab='', cex.axis = 2, cex.names = 2, cex.main = 2)
mtext("F-measure", side = 1, line = 5.5, cex = 1.5, font=2)
legend("topright", as.character(papers), col=cols, fil=cols, cex=2)
text(0.5, -1.0, 'F-measure', cex=2)
dev.off()

######################################################################
#F, Sensitivity and PPV distributions for each paper:

pdf(file="hists-F-per-paper.pdf",height=21,width=10)
par(las=1,cex=4,mfrow=c(3,2),mar=c(5.1+0, 4.1, 4.1+1, 2.1))
xrange<-seq(0,1,length.out=21)
for(j in 1:length(papers)){
      p <- papers[j]
      fp<-meths$F1.measure[p==meths$Paper]
      hv1 <- hist(fp,breaks=xrange,plot=F)$counts
      barplot(hv1,names.arg=format(round(xrange[-1], digits = 2), nsmall = 2),space=0,las=2,col=cols[j],main=paste(p, "F-measure", sep="\n"),xlab='', cex.axis = 2, cex.names = 2, cex.main = 2)
}
dev.off()

pdf(file="hists-Sensitivity-per-paper.pdf",height=21,width=10)
par(las=1,cex=4,mfrow=c(3,2),mar=c(5.1+0, 4.1, 4.1+1, 2.1))
xrange<-seq(0,1,length.out=21)
for(j in 1:length(papers)){
      p <- papers[j]
      fp<-meths$Sensitivity[p==meths$Paper]
      hv1 <- hist(fp,breaks=xrange,plot=F)$counts
      barplot(hv1,names.arg=format(round(xrange[-1], digits = 2), nsmall = 2),space=0,las=2,col=cols[j],main=paste(p, "Sensitivity", sep="\n"),xlab='', cex.axis = 2, cex.names = 2, cex.main = 2)
}
dev.off()

pdf(file="hists-PPV-per-paper.pdf",height=21,width=10)
par(las=1,cex=4,mfrow=c(3,2),mar=c(5.1+0, 4.1, 4.1+1, 2.1))
xrange<-seq(0,1,length.out=21)
for(j in 1:length(papers)){
      p <- papers[j]
      fp<-meths$PPV[p==meths$Paper]
      hv1 <- hist(fp,breaks=xrange,plot=F)$counts
      barplot(hv1,names.arg=format(round(xrange[-1], digits = 2), nsmall = 2),space=0,las=2,col=cols[j],main=paste(p, "PPV", sep="\n"),xlab='', cex.axis = 2, cex.names = 2, cex.main = 2)
}
dev.off()


pdf(file="hists-F-norm-per-paper.pdf",height=21,width=10)
par(las=1,cex=4,mfrow=c(3,2),mar=c(5.1+0, 4.1, 4.1+1, 2.1))
xrange<-seq(-5.4,2.5,length.out=20)   #CAN THROW AN ERROR, AUTOMATE FINDING THE LIMITS????
for(j in 1:length(papers)){
      p <- papers[j]
      fp<- (meths$F1.measure[p==meths$Paper]-median(meths$F1.measure[p==meths$Paper]) )/mad(meths$F1.measure[p==meths$Paper])
      cat(paste(p, min(fp), max(fp), "\n" ))
      hv1 <- hist(fp,breaks=xrange,plot=F)$counts
      barplot(hv1,names.arg=format(round(xrange[-1], digits = 1), nsmall = 1),space=0,las=2,col=cols[j],main=paste(p, "Normalised F-measure", sep="\n"),xlab='', cex.axis = 2, cex.names = 2, cex.main = 2)
}
dev.off()


######################################################################
#Sensitivity distributions for Kraken, Clark & MetaPhyler: 
pdf(file="kraken-clark-metaphyler-hists-Sensitivity.pdf",height=17,width=10)
par(las=1,cex=4,mfrow=c(4,1),mar=c(5.1+0, 4.1, 4.1+1, 2.1))
xrange<-seq(0,1,length.out=21)
remove(cnts)
for(j in 1:length(papers)){
      p <- papers[j]
      fp<-meths$Sensitivity['Kraken'==meths$Method & p==meths$Paper]
      hv1 <- hist(fp,breaks=xrange,plot=F)$counts
      if(exists("cnts")){
	cnts<-rbind(cnts,hv1)
      }
      else{
	cnts<-hv1
      }
}
barplot(cnts,names.arg=format(round(xrange[-1], digits = 2), nsmall = 2),space=0,las=2,col=cols,main='Kraken',xlab='', cex.axis = 2, cex.names = 2, cex.main = 2)
remove(cnts)
for(j in 1:length(papers)){
      p <- papers[j]
      fp<-meths$Sensitivity['CLARK'==meths$Method & p==meths$Paper]
      hv1 <- hist(fp,breaks=xrange,plot=F)$counts
      if(exists("cnts")){
	cnts<-rbind(cnts,hv1)
      }
      else{
	cnts<-hv1
      }
}
barplot(cnts,names.arg=format(round(xrange[-1], digits = 2), nsmall = 2),space=0,las=2,col=cols,main='CLARK',xlab='', cex.axis = 2, cex.names = 2, cex.main = 2)
remove(cnts)
for(j in 1:length(papers)){
      p <- papers[j]
      fp<-meths$Sensitivity['MEGAN'==meths$Method & p==meths$Paper]
      hv1 <- hist(fp,breaks=xrange,plot=F)$counts
      if(exists("cnts")){
	cnts<-rbind(cnts,hv1)
      }
      else{
	cnts<-hv1
      }
}
barplot(cnts,names.arg=format(round(xrange[-1], digits = 2), nsmall = 2),space=0,las=2,col=cols,main='MEGAN',xlab='', cex.axis = 2, cex.names = 2, cex.main = 2)
remove(cnts)
for(j in 1:length(papers)){
      p <- papers[j]
      fp<-meths$Sensitivity['MetaPhyler'==meths$Method & p==meths$Paper]
      hv1 <- hist(fp,breaks=xrange,plot=F)$counts
      if(exists("cnts")){
	cnts<-rbind(cnts,hv1)
      }
      else{
	cnts<-hv1
      }
}
par(mar=c(5.1+2, 4.1, 4.1+1, 2.1))
barplot(cnts,names.arg=format(round(xrange[-1], digits = 2), nsmall = 2),space=0,las=2,col=cols,main='MetaPhyler',xlab='', cex.axis = 2, cex.names = 2, cex.main = 2)
mtext("Sensitivity", side = 1, line = 5.5, cex = 1.5, font=2)
legend("topright", as.character(papers), col=cols, fil=cols, cex=2)
text(0.5, -1.0, 'Sensitivity', cex=2)
dev.off()

######################################################################
#PPV distributions for Kraken, Clark & MetaPhyler: 
pdf(file="kraken-clark-metaphyler-hists-PPV.pdf",height=17,width=10)
par(las=1,cex=4,mfrow=c(4,1),mar=c(5.1+0, 4.1, 4.1+1, 2.1))
xrange<-seq(0,1,length.out=21)
remove(cnts)
for(j in 1:length(papers)){
      p <- papers[j]
      fp<-meths$PPV['Kraken'==meths$Method & p==meths$Paper]
      hv1 <- hist(fp,breaks=xrange,plot=F)$counts
      if(exists("cnts")){
	cnts<-rbind(cnts,hv1)
      }
      else{
	cnts<-hv1
      }
}
barplot(cnts,names.arg=format(round(xrange[-1], digits = 2), nsmall = 2),space=0,las=2,col=cols,main='Kraken',xlab='', cex.axis = 2, cex.names = 2, cex.main = 2)
remove(cnts)
for(j in 1:length(papers)){
      p <- papers[j]
      fp<-meths$PPV['CLARK'==meths$Method & p==meths$Paper]
      hv1 <- hist(fp,breaks=xrange,plot=F)$counts
      if(exists("cnts")){
	cnts<-rbind(cnts,hv1)
      }
      else{
	cnts<-hv1
      }
}
barplot(cnts,names.arg=format(round(xrange[-1], digits = 2), nsmall = 2),space=0,las=2,col=cols,main='CLARK',xlab='', cex.axis = 2, cex.names = 2, cex.main = 2)
remove(cnts)
for(j in 1:length(papers)){
      p <- papers[j]
      fp<-meths$PPV['MEGAN'==meths$Method & p==meths$Paper]
      hv1 <- hist(fp,breaks=xrange,plot=F)$counts
      if(exists("cnts")){
	cnts<-rbind(cnts,hv1)
      }
      else{
	cnts<-hv1
      }
}
barplot(cnts,names.arg=format(round(xrange[-1], digits = 2), nsmall = 2),space=0,las=2,col=cols,main='MEGAN',xlab='', cex.axis = 2, cex.names = 2, cex.main = 2)
remove(cnts)
for(j in 1:length(papers)){
      p <- papers[j]
      fp<-meths$PPV['MetaPhyler'==meths$Method & p==meths$Paper]
      hv1 <- hist(fp,breaks=xrange,plot=F)$counts
      if(exists("cnts")){
	cnts<-rbind(cnts,hv1)
      }
      else{
	cnts<-hv1
      }
}
par(mar=c(5.1+2, 4.1, 4.1+1, 2.1))
barplot(cnts,names.arg=format(round(xrange[-1], digits = 2), nsmall = 2),space=0,las=2,col=cols,main='MetaPhyler',xlab='', cex.axis = 2, cex.names = 2, cex.main = 2)
mtext("PPV", side = 1, line = 5.5, cex = 1.5, font=2)
legend("topright", as.character(papers), col=cols, fil=cols, cex=2)
text(0.5, -1.0, 'PPV', cex=2)
dev.off()



######################################################################
#Citations 
cites<-read.table("citations2.tsv", header=T, sep = "\t")

#EVIL HOOP JUMPING TO COMBINE PAPER COLOURS:
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


library(extrafont)

pdf(file="citations2-year.pdf",family="FreeSans",height=17,width=10)
par(las=3,cex=3,mfrow=c(2,1), mar=c(5, 4+2, 4, 2))
xrange <- 2007:2018-0.5
hv1 <- hist(cites$Year[cites$Evaluated == 'y'], breaks=xrange,plot=F)
hv2 <- hist(cites$Year[cites$Evaluated == 'n'], breaks=xrange,plot=F)
cnts<-rbind(hv1$counts,hv2$counts)
colnames(cnts) <- hv1$mids
rownames(cnts) <- c('Evaluated','Unevaluated')
barplot(cnts,space=0,las=3,xlab='', ylab='Number of eDNA analysis methods', main='', col=c('lightskyblue','peachpuff'), cex.main = 2, cex.lab = 2, cex.axis = 2, cex.names = 2) #cex.names = 2, 
legend("topleft", c('Evaluated','Unevaluated'), col=c('lightskyblue','peachpuff'), fil=c('lightskyblue','peachpuff'), cex=2)
plot(0, ylab="Number of citations",yaxt='n', xlab="",yaxt='n',xaxt='n', xlim=c(2007, 2017.5), ylim=c(0,4),bty='n'  , cex.axis = 2, cex.main = 2, cex.lab = 2)  #cex.names = 2, 
par(las=3)
axis(1, at=2007:2017, labels=c(2007,'',2009,'',2011,'',2013,'',2015,'',2017)  , cex.axis = 2, cex.main = 2) #cex.names = 2, 
axis(2, at=0:4, labels=10^(0:4)                                               , cex.axis = 2, cex.main = 2) #cex.names = 2, 
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

######################################################################
#Sensitivity vs PPV

c25 <- c(
"#E31A1C", # red
"green4",
"#6A3D9A", # purple
"#FF7F00", # orange
"gold1",
"skyblue2","#FB9A99", # lt pink
"palegreen2",
"#CAB2D6", # lt purple
"#FDBF6F", # lt orange
"gray70",
"khaki2",
"deeppink1",
"blue1",
"green1",
"yellow4",
"yellow3",
"darkorange4",
"brown",
"#023fa5", "#7d87b9", "#bb7784", "#8e063b", "#4a6fe3", "#e6afb9", "#e07b91", "#d33f6a", "#11c638", "#8dd593", "#c6dec7", "#f0b98d", "#ef9708", "#0fcfc0", "#f79cd4")

pch<-c(3, 4, 7, 8, 9, 10, 12, 13)
pdf(file="sens-vs-ppv.pdf",family="FreeSans",height=10,width=11.5)
par(las=1,cex=3,mfrow=c(1,1), mar=c(5+3, 4+4, 4, 2))
plot(0, xlab="Positive predictive value",ylab="Sensitivity",xlim=c(0, 1.2), ylim=c(0,1.2),xaxt='n',yaxt='n', cex.axis = 2, bty='n')
axis(1, at=(0:5)/5, labels=c((0:5)/5)  , cex.axis = 2, cex.names = 2, cex.main = 2)
axis(2, at=(0:5)/5, labels=c((0:5)/5)  , cex.axis = 2, cex.names = 2, cex.main = 2)
cnt <- 1;
for(i in 1:length(cites$Method)){
      if(cites$Evaluated[i] == 'y'){
         cat(as.character(cites$Method[i]), c25[cnt], '\n')
	 for(j in 1:length(papers)){
      	     p <- papers[j]
             sens <- meths$Sensitivity[as.character(cites$Method[i])==as.character(meths$Method)  & p==meths$Paper]
             ppv  <- meths$PPV[        as.character(cites$Method[i])==as.character(meths$Method)  & p==meths$Paper]
             points(ppv, sens, col=c25[cnt], cex=0.75, pch=pch[j])
        }
        sens <- meths$Sensitivity[as.character(cites$Method[i])==as.character(meths$Method)]
        ppv  <- meths$PPV[        as.character(cites$Method[i])==as.character(meths$Method)]
	points(median(ppv), median(sens), col=c25[cnt], pch='+', cex=2.5)
	text(median(ppv), median(sens), as.character(cites$Method[i]), cex=1.5, col=c25[cnt], pos=4)
        cnt <- cnt+1
	cnt <- (cnt %% length(c25))+1
      }
}
for(j in 1:length(papers)){
      p <- papers[j]
      xm <- (j %% 2)
      points(0.3*xm, 1.2 -((j+1)%%3)*0.05, pch=pch[j], col='black')
      text(  0.3*xm, 1.2 -((j+1)%%3)*0.05, p, pos=4)
}
dev.off()




##MISSING DATA, CONTINGENCY TABLE:
#library("Skillings.Mack")
#Ski.Mack(meths$F1.measure,groups = meths$Method,blocks = meths$Paper)

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
		for(k in (j+1):(length(methods)-1)){
		      m3 <- methods[k]
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

pdf(file="forest-plot-F-measure.pdf",height=10,width=5,family="FreeSans") #font
forest(mn1, sortvar=-TE, layout="RevMan5", ref="EBI-mg", prediction=TRUE)
dev.off()

metabias(df.mc$TE, seTE=df.mc$seTE, method.bias='rank', plotit=TRUE, correct=TRUE)

######################################################################

system("convert -flatten citations2-year.pdf citations2-year.png")
system("convert -flatten kraken-clark-metaphyler-hists.pdf kraken-clark-metaphyler-hists.png")
system("convert -flatten vizF-matrix-F.pdf vizF-matrix-F.png")
system("convert -flatten vizF-dist-sorted.pdf vizF-dist-sorted.png")
system("convert -flatten vizSens-dist-sorted.pdf vizSens-dist-sorted.png")
system("convert -flatten vizPPV-dist-sorted.pdf  vizPPV-dist-sorted.png")
system("convert -flatten networkDAG-accuracy-ranks.pdf networkDAG-accuracy-ranks.png")
system("convert -flatten forest-plot-F-measure.pdf forest-plot-F-measure.png")
system("convert -flatten vizFnorm-dist-sorted.pdf vizFnorm-dist-sorted.png")
system("convert -flatten hists-F-per-paper.pdf hists-F-per-paper.png")
system("convert -flatten hists-PPV-per-paper.pdf hists-PPV-per-paper.png")
system("convert -flatten hists-Sensitivity-per-paper.pdf hists-Sensitivity-per-paper.png")
system("convert -flatten hists-F-norm-per-paper.pdf hists-F-norm-per-paper.png")

#egrep 'CLARK|MetaPhyler|MEGAN|Kraken' metagenome-meta-analysis-F-measures2.tsv
######################################################################

methodsSm <- c('CLARK', 'MetaPhyler', 'MEGAN', 'Kraken')
for(i in 1:length(methodsSm)){
      for(j in length(papers):1){
      	    p  <- papers[j]
	    fp <- meths$F1.measure[methodsSm[i]==meths$Method & p==meths$Paper]
      	    if(length(fp) > 0){
	    	cat(paste(methodsSm[i], p, median(fp), "\n"))
	    }
      }
}


#cat metagenome-meta-analysis-F-measures2-small.tsv | sort -k2,2d -k3,3nr

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



















