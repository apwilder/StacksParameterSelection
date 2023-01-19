# This R script will plot summary statistics output from StacksDeNovo.sh run in parameter exploration mode
# Run from main stacks output directory with the command, eg: Rscript PlotStacksSummary.R BFI 3 2 0 1 8 1 8
args <- commandArgs(trailingOnly = TRUE)

species=args[1]

mdef=args[2]
Mdef=args[3]
ndef=args[4]
mrng=args[5]:args[6]
Mrng=args[7]:args[8]

#
if(species=='Waterbuck'){
	nsamps=nrow(read.table(paste(c(species,'_m',mdef,'M',Mdef,'n',ndef,'_ustack.stats'),collapse=''),header=T))
	} else {
	nsamps=nrow(read.table(paste(c(species,'_m',mdef,'M',Mdef,'_ustack.stats'),collapse=''),header=T))
	}

####Varying m
covg=matrix(nrow=nsamps,ncol=max(mrng))
totreads=matrix(nrow=nsamps,ncol=max(mrng))
totstacks=matrix(nrow=nsamps,ncol=max(mrng))
for (m in mrng){
if(species=='Waterbuck'){
	ustack=read.table(paste(c(species,'_m',m,'M',Mdef,'n0_ustack.stats'),collapse=''),header=T)
	}else{
	ustack=read.table(paste(c(species,'_m',m,'M',Mdef,'_ustack.stats'),collapse=''),header=T)}

covg[,m]=ustack$final_cvg_mean
totreads[,m]=ustack$final_reads_in_stacks
totstacks[,m]=ustack$assmbld_stacks
}

shrd=read.table(paste0(species,'_SharedPolymorphisms.txt'),header=T)
png(file=paste0(species,"_meanCovg_m.png"),res=300,width=6,height=6,unit='in')
boxplot(covg,xlab='m',ylab='Mean Coverage after ustacks',main=paste0(species))
dev.off()
png(file=paste0(species,"_ReadsInStacks_m.png"),res=300,width=6,height=6,unit='in')
boxplot(totreads,xlab='m',ylab='Reads in stacks after ustacks',main=paste0(species))
dev.off()

idx=which(shrd$n==ndef & shrd$M==Mdef)
shrdi=shrd[idx,]
shrdi=shrdi[order(shrdi$m),]
png(file=paste0(species,"_Nloci_m.png"),res=300,width=6,height=6,unit='in')
boxplot(totstacks,xlab='m',ylab='Number of loci',main=paste(c(species,': optimal m=',shrdi$m[which(shrdi$Loci_80pctShared==max(shrdi$Loci_80pctShared))],': ',max(shrdi$Loci_80pctShared),' loci'),collapse=''))
points(shrdi$m,shrdi$Loci_50pctShared,col='darkorange',pch=16)
points(shrdi$m,shrdi$Loci_80pctShared,col='dodgerblue',pch=16)
abline(h=max(shrdi$Loci_80pctShared),lty=2)
legend('topright',legend=c('50% of samples','80% of samples'),col=c('darkorange','dodgerblue'),pch=16,bty='n')
dev.off()

png(file=paste0(species,"_SharedPolymorphicLoci_m.png"),res=300,width=6,height=6,unit='in')
barplot(height=as.matrix(t(cbind(shrdi$PolymLoci_50pctShared,shrdi$PolymLoci_80pctShared))),beside=T,names=as.factor(shrdi$m),xlab='m',ylab='Polymorphic Loci',col=c('darkorange','dodgerblue'),legend.text=c('50% of samples','80% of samples'),args.legend=list(bty='n'),main=paste(c(species,': optimal m=',shrdi$m[which(shrdi$PolymLoci_80pctShared==max(shrdi$PolymLoci_80pctShared))],': ',max(shrdi$PolymLoci_80pctShared),' Polymorphic loci'),collapse=''))
abline(h=max(shrdi$PolymLoci_80pctShared),lty=2)
dev.off()

png(file=paste0(species,"_SharedSNPs_m.png"),res=300,width=6,height=6,unit='in')
barplot(height=as.matrix(t(cbind(shrdi$SNPs_50pctShared,shrdi$SNPs_80pctShared))),beside=T,names=as.factor(shrdi$m),xlab='m',ylab='SNPs',col=c('darkorange','dodgerblue'),legend.text=c('50% of samples','80% of samples'),args.legend=list(bty='n'),main=paste(c(species,': optimal m=',shrdi$m[which(shrdi$SNPs_80pctShared==max(shrdi$SNPs_80pctShared))],': ',max(shrdi$SNPs_80pctShared),' SNPs'),collapse=''))
abline(h=max(shrdi$SNPs_80pctShared),lty=2)
dev.off()

####Varying M
covg=matrix(nrow=nsamps,ncol=max(mrng))
totreads=matrix(nrow=nsamps,ncol=max(mrng))
totstacks=matrix(nrow=nsamps,ncol=max(mrng))
for (M in Mrng){
	if(species=='Waterbuck'){
ustack=read.table(paste(c(species,'_m',mdef,'M',M,'n0_ustack.stats'),collapse=''),header=T)}
	else{
ustack=read.table(paste(c(species,'_m',mdef,'M',M,'_ustack.stats'),collapse=''),header=T)}
covg[,M]=ustack$final_cvg_mean
totreads[,M]=ustack$final_reads_in_stacks
totstacks[,M]=ustack$assmbld_stacks
}

png(file=paste0(species,"_meanCovg_bigM.png"),res=300,width=6,height=6,unit='in')
boxplot(covg,xlab='M',ylab='Mean Coverage after ustacks',main=paste0(species))
dev.off()
png(file=paste0(species,"_ReadsInStacks_bigM.png"),res=300,width=6,height=6,unit='in')
boxplot(totreads,xlab='M',ylab='Reads in stacks after ustacks',main=paste0(species))
dev.off()

idx=which(shrd$n==ndef & shrd$m==mdef)
shrdi=shrd[idx,]
shrdi=shrdi[order(shrdi$M),]
png(file=paste0(species,"_Nloci_bigM.png"),res=300,width=6,height=6,unit='in')
boxplot(totstacks,xlab='M',ylab='Number of loci',main=paste(c(species,': optimal M=',shrdi$M[which(shrdi$Loci_80pctShared==max(shrdi$Loci_80pctShared))],': ',max(shrdi$Loci_80pctShared),' Loci'),collapse=''),ylim=c(0,max(totstacks,na.rm=T)+50000))
points(shrdi$M,shrdi$Loci_50pctShared,col='darkorange',pch=16)
points(shrdi$M,shrdi$Loci_80pctShared,col='dodgerblue',pch=16)
legend('topright',legend=c('50% of samples','80% of samples'),col=c('darkorange','dodgerblue'),pch=16,bty='n')
abline(h=max(shrdi$Loci_80pctShared),lty=2)
dev.off()

png(file=paste0(species,"_SharedSNPs_bigM.png"),res=300,width=6,height=6,unit='in')
barplot(height=as.matrix(t(cbind(shrdi$SNPs_50pctShared,shrdi$SNPs_80pctShared))),beside=T,names=as.factor(shrdi$M),xlab='M',ylab='SNPs',col=c('darkorange','dodgerblue'),legend.text=c('50% of samples','80% of samples'),args.legend=list(bty='n'),main=paste(c(species,': optimal M=',shrdi$M[which(shrdi$SNPs_80pctShared==max(shrdi$SNPs_80pctShared))],': ',max(shrdi$SNPs_80pctShared),' SNPs'),collapse=''),ylim=c(0,100000))
abline(h=max(shrdi$SNPs_80pctShared),lty=2)
dev.off()

png(file=paste0(species,"_SharedPolymorphicLoci_bigM.png"),res=300,width=6,height=6,unit='in')
barplot(height=as.matrix(t(cbind(shrdi$PolymLoci_50pctShared,shrdi$PolymLoci_80pctShared))),beside=T,names=as.factor(shrdi$M),xlab='M',ylab='Polymorphic Loci',col=c('darkorange','dodgerblue'),legend.text=c('50% of samples','80% of samples'),args.legend=list(bty='n'),main=paste(c(species,': optimal M=',shrdi$M[which(shrdi$PolymLoci_80pctShared==max(shrdi$PolymLoci_80pctShared))],': ',max(shrdi$PolymLoci_80pctShared),' Polymorphic loci'),collapse=''),ylim=c(0,100000))
abline(h=max(shrdi$PolymLoci_80pctShared),lty=2)
dev.off()

####Varying n

idx=which(shrd$m==mdef & shrd$M==Mdef)
shrdi=shrd[idx,]
shrdi=shrdi[order(shrdi$n),]
png(file=paste0(species,"_SharedSNPs_n.png"),res=300,width=6,height=6,unit='in')
barplot(height=as.matrix(t(cbind(shrdi$SNPs_50pctShared,shrdi$SNPs_80pctShared))),beside=T,names=as.factor(shrdi$n),xlab='n',ylab='SNPs',col=c('darkorange','dodgerblue'),legend.text=c('50% of samples','80% of samples'),args.legend=list(bty='n'),main=paste(c(species,': optimal n=',shrdi$n[which(shrdi$SNPs_80pctShared==max(shrdi$SNPs_80pctShared))],': ',max(shrdi$SNPs_80pctShared),' SNPs'),collapse=''),ylim=c(0,120000))
abline(h=max(shrdi$SNPs_80pctShared),lty=2)
dev.off()

png(file=paste0(species,"_SharedPolymorphicLoci_n.png"),res=300,width=6,height=6,unit='in')
barplot(height=as.matrix(t(cbind(shrdi$PolymLoci_50pctShared,shrdi$PolymLoci_80pctShared))),beside=T,names=as.factor(shrdi$n),xlab='n',ylab='Polymorphic Loci',col=c('darkorange','dodgerblue'),legend.text=c('50% of samples','80% of samples'),args.legend=list(bty='n'),main=paste(c(species,': optimal n=',shrdi$n[which(shrdi$PolymLoci_80pctShared==max(shrdi$PolymLoci_80pctShared))],': ',max(shrdi$PolymLoci_80pctShared),' Polymorphic loci'),collapse=''),ylim=c(0,120000))
abline(h=max(shrdi$PolymLoci_80pctShared),lty=2)
dev.off()

####Including Reference-based
idx=which(shrd$n==ndef & shrd$M==Mdef | shrd$m=='Ref')
shrdi=shrd[idx,]
shrdi=shrdi[order(shrdi$m),]

png(file=paste0(species,"_SharedPolymorphicLoci_m_vRef.png"),res=300,width=6,height=6,unit='in')
barplot(height=as.matrix(t(cbind(shrdi$PolymLoci_50pctShared,shrdi$PolymLoci_80pctShared))),beside=T,names=as.factor(shrdi$m),xlab='m',ylab='Polymorphic Loci',col=c('darkorange','dodgerblue'),legend.text=c('50% of samples','80% of samples'),args.legend=list(bty='n'),main=paste(c(species,': optimal m=',as.character(shrdi$m[which(shrdi$PolymLoci_80pctShared==max(shrdi$PolymLoci_80pctShared))]),': ',max(shrdi$PolymLoci_80pctShared),' Polymorphic loci'),collapse=''))
abline(h=max(shrdi$PolymLoci_80pctShared),lty=2)
dev.off()

png(file=paste0(species,"_SharedSNPs_m_vRef.png"),res=300,width=6,height=6,unit='in')
barplot(height=as.matrix(t(cbind(shrdi$SNPs_50pctShared,shrdi$SNPs_80pctShared))),beside=T,names=as.factor(shrdi$m),xlab='m',ylab='SNPs',col=c('darkorange','dodgerblue'),legend.text=c('50% of samples','80% of samples'),args.legend=list(bty='n'),main=paste(c(species,': optimal m=',as.character(shrdi$m[which(shrdi$SNPs_80pctShared==max(shrdi$SNPs_80pctShared))]),': ',max(shrdi$SNPs_80pctShared),' SNPs'),collapse=''))
abline(h=max(shrdi$SNPs_80pctShared),lty=2)
dev.off()