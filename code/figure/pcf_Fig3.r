# 27/01/21 CP: uses outputs from main_Fig3.cpp to show pcf=f(r/Delta) for different U_tot

rm(list=ls())
graphics.off()

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


library(expint)

Delta=10^(-7)

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


pdf("pcf_test_Utot_modif_dx_dp.pdf",width=10,height=7.5)
par(mfrow=c(1,2),mar=c(4,4.5,1,1))
f0=read.table('../simulation/pcf_dpow0p25_area10_tmax1000_N200000_U0p0.txt',sep=";",header=T,dec=".")
f0$rDelta=f0$r/Delta
f_ab0=read.table("../simulation/nb_individuals_dpow0p25_area10_tmax1000_N200000_U0p0.txt",sep=";",header=T,dec=".")

plot(0.1,0.1,t="n",xlim=range(f0$rDelta),ylim=c(0.01,max(f0$pcf_dx)),log="xy",xlab=expression(r/Delta),ylab="g(r,t)-1",axes=F)

axis(1, at=log10Tck('x','major'), tcl= 0.2) # bottom
axis(2, at=log10Tck('y','major'), tcl= 0.2) # left
axis(1, at=log10Tck('x','minor'), tcl= 0.1, labels=NA) # bottom
axis(2, at=log10Tck('y','minor'), tcl= 0.1, labels=NA) # left
box()

#U=0
lines(f0$rDelta,f0$pcf_dx-1,lty=1,col="black",lwd=2)
#lines(f0$rDelta,f0$pcf_dp-1,lty=1,col="grey",lwd=2,t="o")
lines(f0$rDelta,G_theoretical(0.0,f0$rDelta,f_ab0$Nb_ind/f_ab0$area,"0",1000),col="black",lty=2,lwd=2)

f1=read.table('../simulation/pcf_dpow0p25_area10_tmax1000_N200000_U0p1.txt',sep=";",header=T,dec=".")
f1$rDelta=f1$r/Delta
f_ab1=read.table("../simulation/nb_individuals_dpow0p25_area10_tmax1000_N200000_U0p1.txt",sep=";",header=T,dec=".")

lines(f1$rDelta,f1$pcf_dx-1,lty=1,col="blue",lwd=2)
#lines(f1$rDelta,f1$pcf_dp-1,lty=1,col="cyan",lwd=2)
lines(f1$rDelta,G_theoretical(0.026,f1$rDelta,f_ab1$Nb_ind/f_ab1$area,"1",1000),col="blue",lty=2,lwd=2)

f2=read.table('../simulation/pcf_dpow0p25_area10_tmax1000_N200000_U0p5.txt',sep=";",header=T,dec=".")
f2$rDelta=f2$r/Delta
f_ab2=read.table("../simulation/nb_individuals_dpow0p25_area10_tmax1000_N200000_U0p5.txt",sep=";",header=T,dec=".")

lines(f2$rDelta,f2$pcf_dx-1,lty=1,col="red",lwd=2)
#lines(f2$rDelta,f2$pcf_dp-1,lty=1,col="pink",lwd=2)
lines(f2$rDelta,G_theoretical(0.51,f2$rDelta,f_ab2$Nb_ind/f_ab2$area,"1",1000),col="red",lty=2,lwd=2)

f3=read.table('../simulation/pcf_dpow0p25_area10_tmax1000_N200000_U2p5.txt',sep=";",header=T,dec=".")
f3$rDelta=f3$r/Delta
f_ab3=read.table("../simulation/nb_individuals_dpow0p25_area10_tmax1000_N200000_U2p5.txt",sep=";",header=T,dec=".")

lines(f3$rDelta,f3$pcf_dx-1,lty=1,col="darkgreen",lwd=2)
#lines(f3$rDelta,f3$pcf_dp-1,lty=1,col="green",lwd=2)
lines(f3$rDelta,G_theoretical(2.4,f3$rDelta,f_ab3$Nb_ind/f_ab3$area,"1",1000),col="darkgreen",lty=2,lwd=2)

legend("bottomleft",c(expression(paste("U",tau,"/2=0",sep="")),expression(paste("U",tau,"/2=0.1",sep="")),expression(paste("U",tau,"/2=0.5",sep="")),expression(paste("U",tau,"/2=2.5",sep="")),"Analytical solution"),col=c("black","blue","red","darkgreen","black"),lty=c(1,1,1,1,2,2),bty="n",lwd=2)

#Drawing the -2 line
y1=10^7
x1=f0$rDelta[15]
x2=f0$rDelta[21]
y2=x2^(-2)*(y1/x1^(-2))
lines(c(x1,x2),c(y1,y2),lwd=3)
text((x1+x2)*0.15,(y1+y2)/2,"-2")

##Zoom in
plot(0,0,t="n",xlim=c(0,20),ylim=c(0.,3*10^9),xlab="r/Delta",ylab=expression(paste("g(r,t)-1 (x 10"^"9",")")),yaxt="n")
axis(2,at=c(0,10^9,2*10^9,3*10^9),labels=c(0,1,2,3))

#U=0
lines(f0$rDelta,f0$pcf_dx-1,lty=1,col="black",lwd=2)
lines(f0$rDelta,G_theoretical(0.0,f0$rDelta,f_ab0$Nb_ind/f_ab0$area,"0",1000),col="black",lty=2,lwd=2)

#U=0.1
lines(f1$rDelta,f1$pcf_dx-1,lty=1,col="blue",lwd=2)
lines(f1$rDelta,G_theoretical(0.026,f1$rDelta,f_ab1$Nb_ind/f_ab1$area,"1"),col="blue",lty=2,lwd=2)

#U=0.5
lines(f2$rDelta,f2$pcf_dx-1,lty=1,col="red",lwd=2)
lines(f2$rDelta,G_theoretical(0.51,f2$rDelta,f_ab2$Nb_ind/f_ab2$area,"1"),col="red",lty=2,lwd=2)

#U=2.5
lines(f3$rDelta,f3$pcf_dx-1,lty=1,col="darkgreen",lwd=2)
lines(f3$rDelta,G_theoretical(2.4,f3$rDelta,f_ab3$Nb_ind/f_ab3$area,"1"),col="darkgreen",lty=2,lwd=2)
dev.off()
