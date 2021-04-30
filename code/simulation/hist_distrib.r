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
G_theoretical=function(gamma,rDelta,C_0,U="1",tmax=NA){ #U=0 in the absence of advection
        tau=1
        lambda=0.5
        if(U!="0"){
        tmp=-1*lambda*tau/(2*pi*C_0*Delta^2)*(log(gamma)-log(1/(tau*rDelta^2)+gamma)) #Assuming that the advection U > 0
        }else{ #This corresponds to the case U=0
        D=Delta^2/(2*tau)
        tmp=2*lambda/C_0*(expint((rDelta*Delta)^2/(8*tmax*D))/(8*pi*D)) #We are multiplying rDelta by Delta because the input file contains r/Delta, no just r
        }
        return(tmp)
}


par(mfrow=c(1,2))
f=read.table('pcf_dpow0p25_tmax1000_U0p1_test.txt',sep=";",header=F,dec=".")
colnames(f)=c("rDelta","U","pcf")

f=unique(f)

plot(0.1,0.1,t="n",xlim=range(Delta*f$rDelta),ylim=c(0.01,max(f$pcf)),log="xy",xlab="r",ylab="g(r,t)-1",axes=F,main="Utau/2=0.1")

axis(1, at=log10Tck('x','major'), tcl= 0.2) # bottom
axis(2, at=log10Tck('y','major'), tcl= 0.2) # left
axis(1, at=log10Tck('x','minor'), tcl= 0.1, labels=NA) # bottom
axis(2, at=log10Tck('y','minor'), tcl= 0.1, labels=NA) # left
box()

f0=subset(f,U==0.1)
lines(Delta*f0$rDelta,f0$pcf-1,lty=1,col="black",lwd=2,t="o")
lines(Delta*f0$rDelta,G_theoretical(0.026,f0$rDelta,19450/area,"1",1000),col="black",lty=2,lwd=2)

tt=read.table("nb_individuals_dpow0p25_intervalle_tmax1000_U0p1_test.txt",sep=";",header=F,skip=1)
colnames(tt)=c("Utot","pow","nb")

U_unique=unique(tt$Utot)

for(U in 1:length(U_unique)){
	plot(-1*tt$pow,tt$nb/2,log="y",main=paste("Utau/2=",U_unique[U]),xlab="dist ij",ylab="Nb",xaxis="n")
	#xlim=c(10^(-10),0.01)
}
