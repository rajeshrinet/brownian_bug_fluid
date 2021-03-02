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
N_0=200000
area=10
C_0=N_0/area

G_theoretical=function(gamma,rDelta,C_0,U="1"){
	tau=1
	lambda=0.5
	if(U!="0"){
	tmp=-1*lambda*tau/(2*pi*C_0*Delta^2)*(log(gamma)-log(1/(tau*rDelta^2)+gamma)) #Assuming that the advection U > 0
	}else{ #This corresponds to the case U=0
	D=Delta^2/(2*tau)
	tmp=2*lambda/C_0*(expint((rDelta*Delta)^2/(8*1000*D))/(8*pi*D)) #We are multiplying rDelta by Delta because the input file contains r/Delta, no just r
	}
	return(tmp)
}


pdf("pcf_per_Utot_dx10m8.pdf",width=10,height=7.5)
par(mfrow=c(1,2),mar=c(4,4.5,1,1))
f=read.table('../simulation/pcf_particle_Fig3_dx10m8_big_area_0.txt',sep=";",header=F,dec=".")
colnames(f)=c("rDelta","U","pcf")

f=unique(f)

plot(0.1,0.1,t="n",xlim=range(f$rDelta),ylim=c(0.01,max(f$pcf)),log="xy",xlab="r/Delta",ylab="g(r,t)-1",axes=F)

axis(1, at=log10Tck('x','major'), tcl= 0.2) # bottom
axis(2, at=log10Tck('y','major'), tcl= 0.2) # left
axis(1, at=log10Tck('x','minor'), tcl= 0.1, labels=NA) # bottom
axis(2, at=log10Tck('y','minor'), tcl= 0.1, labels=NA) # left
box()

print("1")
#U=0
f0=subset(f,U==0)
lines(f0$rDelta,f0$pcf-1,lty=1,col="black",lwd=2)
lines(f0$rDelta,G_theoretical(0.0,f0$rDelta,C_0,"0"),col="black",lty=2,lwd=2)


f=read.table('../simulation/pcf_particle_Fig3_dx10m8_big_area_0p1.txt',sep=";",header=F,dec=".")
colnames(f)=c("rDelta","U","pcf")
print("2")
#U=0.1
f1=subset(f,U==0.1)
lines(f1$rDelta,f1$pcf-1,lty=1,col="blue",lwd=2)
lines(f1$rDelta,G_theoretical(0.026,f1$rDelta,C_0,"1"),col="blue",lty=2,lwd=2)

print("3")
f=read.table('../simulation/pcf_particle_Fig3_dx10m8_big_area_0p5.txt',sep=";",header=F,dec=".")
colnames(f)=c("rDelta","U","pcf")
#U=0.5
f2=subset(f,U==0.5)
lines(f2$rDelta,f2$pcf-1,lty=1,col="red",lwd=2)
lines(f2$rDelta,G_theoretical(0.51,f2$rDelta,C_0,"2"),col="red",lty=2,lwd=2)

print("4")
#U=2.5
f=read.table('../simulation/pcf_particle_Fig3_dx10m8_big_area_2p5.txt',sep=";",header=F,dec=".")
colnames(f)=c("rDelta","U","pcf")
f3=subset(f,U==2.5)
lines(f3$rDelta,f3$pcf-1,lty=1,col="darkgreen",lwd=2)
lines(f3$rDelta,G_theoretical(2.4,f3$rDelta,C_0,"3"),col="darkgreen",lty=2,lwd=2)

legend("bottomleft",c('U=0','U=0.1','U=0.5','U=2.5',"Theory"),col=c("black","blue","red","darkgreen","black"),lty=c(1,1,1,1,2,2),bty="n")

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
lines(f0$rDelta,f0$pcf-1,lty=1,col="black",lwd=2)
lines(f0$rDelta,G_theoretical(0.0,f0$rDelta,C_0,"0"),col="black",lty=2,lwd=2)

#U=0.1
lines(f1$rDelta,f1$pcf-1,lty=1,col="blue",lwd=2)
lines(f1$rDelta,G_theoretical(0.026,f1$rDelta,C_0,"1"),col="blue",lty=2,lwd=2)

#U=0.5
lines(f2$rDelta,f2$pcf-1,lty=1,col="red",lwd=2)
lines(f2$rDelta,G_theoretical(0.51,f2$rDelta,C_0,"1"),col="red",lty=2,lwd=2)

#U=2.5
lines(f3$rDelta,f3$pcf-1,lty=1,col="darkgreen",lwd=2)
lines(f3$rDelta,G_theoretical(2.4,f3$rDelta,C_0,"1"),col="darkgreen",lty=2,lwd=2)
dev.off()

