#analysis of Edze data
setwd("/Users/stevep11/Dropbox/bioinf/edze")

targets.bed <- read.table("target_seqs1.bed",sep="\t",col.names=c("clone","start0","stop"))
allcnts.list <- vector(mode="list",length=9)
for(i in 1:9) allcnts.list[[i]] <- read.table(dir()[grep('^Sample[0-9]_pileup',dir())][i],sep="\t",col.names=c("chr","pos","ref","depth","refcnt","altcnt","ncnt","indelcnt"))
#all the same size

targets.bed.unique <- targets.bed[match(unique(targets.bed$start0),targets.bed$start0),]

bedcnts.list <- vector(mode="list",length=9)
for(i in 1:9){
	tmp <- allcnts.list[[i]][0,]
	for(j in 1:nrow(targets.bed.unique)) tmp <- rbind(tmp,allcnts.list[[i]][allcnts.list[[i]]$pos > targets.bed.unique[j,"start0"] & allcnts.list[[i]]$pos <= targets.bed.unique[j,"stop"],])
	bedcnts.list[[i]] <- tmp
}

targets.maxfreq <- targets.bed.unique

findMaxFreq <- function(data,start0,stop){
	tmp <- data[data$pos>start0 & data$pos <= stop,]
	tmp$freq <- (tmp$depth - tmp$refcnt)/tmp$depth
	max(tmp$freq,na.rm=TRUE)
}

for(i in 1:9){
	sample <- paste("Sample",i,sep="")
	targets.maxfreq[,sample] <- 0
	for(j in 1:nrow(targets.maxfreq)) targets.maxfreq[j,sample] <- findMaxFreq(bedcnts.list[[i]],start0=targets.maxfreq[j,"start0"],stop=targets.maxfreq[j,"stop"])	
}

write.table(targets.maxfreq,file="targets.maxfreq.csv",sep=",",row.names=FALSE)

#need to sum the frequencies

sumFreq <- function(data,start0,stop){
	tmp <- data[data$pos>start0 & data$pos <= stop,]	
	tmp$freq <- tmp$depth - tmp$refcnt
	nrow(tmp)*sum(tmp$freq)/sum(tmp$depth)
}

targets.sumfreq <- targets.bed.unique
for(i in 1:9){
	sample <- paste("Sample",i,sep="")
	targets.sumfreq[,sample] <- 0
	for(j in 1:nrow(targets.sumfreq)) targets.sumfreq[j,sample] <- sumFreq(bedcnts.list[[i]],start0=targets.sumfreq[j,"start0"],stop=targets.sumfreq[j,"stop"])	
}
#sample 8 seems quite high, but a smattering of mutations all over

#do what Edze suggests and just look at PAM and 8bp flanking

targets.bed10 <- read.table("target_seqs_10.bed",sep="\t",col.names=c("clone","start0","stop"))
targets.bed10.unique <- targets.bed10[match(unique(targets.bed10$start0),targets.bed10$start0),]
bedcnts10.list <- vector(mode="list",length=9)
for(i in 1:9){
	tmp <- allcnts.list[[i]][0,]
	for(j in 1:nrow(targets.bed10.unique)) tmp <- rbind(tmp,allcnts.list[[i]][allcnts.list[[i]]$pos > targets.bed10.unique[j,"start0"] & allcnts.list[[i]]$pos <= targets.bed10.unique[j,"stop"],])
	bedcnts10.list[[i]] <- tmp
}

targets.maxfreq10 <- targets.bed10.unique
targets.sumfreq10 <- targets.bed10.unique
for(i in 1:9){
	sample <- paste("Sample",i,sep="")
	targets.maxfreq10[,sample] <- 0
	for(j in 1:nrow(targets.maxfreq10)) targets.maxfreq10[j,sample] <- findMaxFreq(bedcnts10.list[[i]],start0=targets.maxfreq10[j,"start0"],stop=targets.maxfreq10[j,"stop"])	
}
for(i in 1:9){
	sample <- paste("Sample",i,sep="")
	targets.sumfreq10[,sample] <- 0
	for(j in 1:nrow(targets.sumfreq10)) targets.sumfreq10[j,sample] <- sumFreq(bedcnts10.list[[i]],start0=targets.sumfreq10[j,"start0"],stop=targets.sumfreq10[j,"stop"])	
}

#reshape this to plot
targets.sumfreq10.long <- cbind(data.frame(clone=rep(targets.sumfreq10$clone,9)),stack(targets.sumfreq10[,grep('Sample',names(targets.sumfreq10))]))
plot(values~ind,data=targets.sumfreq10.long[targets.sumfreq10.long$ind %in% paste("Sample",c(1,5:9),sep=""),])

#genetic distance
#sum(((bedcnts10.list[[5]]$depth - bedcnts10.list[[5]]$refcnt)/bedcnts10.list[[5]]$depth)^2)
#heat map with cut offs at 0.01 and 0.02
heat.mat <- matrix(0,nrow=nrow(targets.sumfreq10),ncol=5)
#heat.mat[targets.sumfreq10[,8:12]>0.01] <- 0.01
heat.mat[targets.sumfreq10[,8:12]>0.02] <- 0.02
rownames(heat.mat) <- targets.sumfreq10$clone
colnames(heat.mat) <- c(1,6,12,24,48)
heatmap3(heat.mat[seq(nrow(heat.mat),1),],Rowv=NA,Colv=NA,scale="none",col=c("darkgreen","red"),margins=c(4,9),xlab="No of clones")

write.table(targets.sumfreq10,file="targets.sumfreq10.csv",sep=",",row.names=FALSE)