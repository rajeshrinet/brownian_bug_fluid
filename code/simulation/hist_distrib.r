rm(list=ls())
graphics.off()

library(expint)

Delta=10^(-7)
N_0=20000
area=10
C_0=N_0/area

#Ticks for log scale from https://stackoverflow.com/questions/47890742/logarithmic-scale-plot-in-r
log10Tck <- function(side, type){
   lim <- switch(side,
     x = par('usr')[1:2],
     y = par('usr')[3:4],
     stop("side argument must be 'x' or 'y'"))
   at <- floor(lim[1]) : ceiling(lim[2])
   return(switch(type,
     minor = outer(1:9, 10^(min(at):max(at))),
     major = 10^at,
     stop("type argument must be 'major' or 'minor'")
   ))
}


#Theoretical value of G
G_theoretical=function(gamma,rDelta,C_0,tmax=NA){ #U=0 in the absence of advection
        tau=1
        lambda=0.5
        if(gamma>0){
        tmp=-1*lambda*tau/(2*pi*C_0*Delta^2)*(log(gamma)-log(1/(tau*rDelta^2)+gamma)) #Assuming that the advection U > 0
        }else{ #This corresponds to the case U=0
        D=Delta^2/(2*tau)
        tmp=2*lambda/C_0*(expint((rDelta*Delta)^2/(8*tmax*D))/(8*pi*D)) #We are multiplying rDelta by Delta because the input file contains r/Delta, no just r
        }
        return(tmp)
}


pdf("pcf_and_distribution_pairdist_1000part.pdf",height=4,width=10)
par(mfrow=c(1,3))

U_unique=c("tmax1000_U0p0","tmax5000_U0p0","tmax1000_U2p5")
val_U=c(0,0,2.5)
tmax=c(1000,5000,1000)

U_unique=c("tmax2000_U0p0")
val_U=c(0)
tmax=c(2000)


for (u in 1:length(U_unique)){

	f=read.table(paste('pcf_dpow0p25_',U_unique[u],'_table_distance.txt',sep=""),sep=";",header=F,dec=".")
	colnames(f)=c("rDelta","U","pcf")

	tt=read.table(paste("nb_individuals_dpow0p25_intervalle_",U_unique[u],"_table_distance.txt",sep=""),sep=";",header=F,nrows=1)
	nb_ind=tt[1,2]
	C_0=nb_ind/area

	plot(0.1,0.1,t="n",xlim=range(f$rDelta),ylim=c(0.01,max(f$pcf)),log="xy",xlab="r",ylab="g(r,t)-1",axes=F,main=paste("Utau/2=",val_U[u],", tmax=",tmax[u],sep=""),cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
	axis(1, at=log10Tck('x','major'), tcl= 0.2,cex.axis=1.5) # bottom
	axis(2, at=log10Tck('y','major'), tcl= 0.2,cex.axis=1.5) # left
	axis(1, at=log10Tck('x','minor'), tcl= 0.1, labels=NA) # bottom
	axis(2, at=log10Tck('y','minor'), tcl= 0.1, labels=NA) # left
	#plot(0.1,0.1,t="n",xlim=range(log10(f$rDelta)),ylim=c(0.01,max(f$pcf)),log="y",xlab="log10(r)",ylab="g(r,t)-1",axes=T,main=paste("Utau/2=",val_U[u],", tmax=",tmax[u],sep=""),cex.axis=1.5,cex.lab=1.5,cex.main=1.5)

	box()

	#lines(log10(Delta*f$rDelta),f$pcf-1,lty=1,col="black",lwd=2,t="o")
	#lines(log10(Delta*f$rDelta),G_theoretical(val_U[u],f$rDelta,C_0,tmax[u]),col="black",lty=2,lwd=2)
	lines(f$rDelta,f$pcf-1,lty=1,col="black",lwd=2,t="o")
	lines(f$rDelta,G_theoretical(val_U[u],f$rDelta/Delta,C_0,tmax[u]),col="black",lty=2,lwd=2)

	tt=read.table(paste("nb_individuals_dpow0p25_intervalle_",U_unique[u],"_table_distance.txt",sep=""),sep=";",header=F,skip=1)
	colnames(tt)=c("Utot","pow","nb")

	plot(-1*tt$pow,tt$nb/2,log="y",xlab="dist ij",ylab="Nb pairs",xaxt="n",t="h",lwd=2,cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
	axis(1,label=c("<=10^(-10)",paste(10^(seq(-8,0,by=2)))),at=seq(-10,0,by=2),cex.axis=1.5,cex.lab=1.5)
	#xlim=c(10^(-10),0.01)

tt=read.table("distance_table.txt",sep=";",header=T)
colnames(tt)=c("i","j","dist")
tt=unique(tt)

hist(log10(tt$dist))
axis(1,seq(-15,0,by=1))
}

dev.off()
