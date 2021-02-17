# 27/01/21 CP: uses outputs from main_Fig3.cpp to show pcf=f(r/Delta) for different U_tot

rm(list=ls())
graphics.off()

library(expint)

Delta=10^(-7)
N_0=20000
G_theoretical=function(gamma,rDelta,C_0,B="0"){
	tau=1
	lambda=0.5
	if(B=="0"){
	tmp=-1*lambda*tau/(2*pi*C_0*Delta^2)*(log(gamma)-log(1/(tau*rDelta^2)+gamma)) #Assuming B=0
	}else{ #This corresponds to the case U=0
	D=Delta^2/(2*tau)
	tmp=2*lambda/C_0*(expint((rDelta*Delta)^2/(8*1000*D))/(8*pi*D)) #We are multiplying rDelta by Delta because the input file contains r/Delta, no just r
	}
	return(tmp)
}


pdf("pcf_per_Utot_dx10m8.pdf",width=10,height=5)
par(mfrow=c(1,2))
f=read.table('pcf_particle_Fig3_dx10m8_small_area_high_U.txt',sep=";",header=F,dec=".")
colnames(f)=c("rDelta","U","pcf")
area=0.1
C_0=N_0/area

f=unique(f)

#plot(0,0,t="n",xlim=range(f$rDelta),ylim=range(f$pcf),log="xy") #pcf=0 for U=0 and r/Delta>10^2
plot(0.1,0.1,t="n",xlim=range(f$rDelta),ylim=c(0.01,max(f$pcf)),log="xy",xlab="r/Delta",ylab="g(r,t)-1")

print("1")
#U=0
f0=subset(f,U==0)
lines(f0$rDelta,f0$pcf-1,lty=1,col="black",lwd=2)
lines(f0$rDelta,G_theoretical(0.0,f0$rDelta,C_0,"1"),col="black",lty=2,lwd=2)

print("2")
#U=0.1
f1=subset(f,U==0.1)
lines(f1$rDelta,f1$pcf-1,lty=1,col="blue",lwd=2)
lines(f1$rDelta,G_theoretical(0.026,f1$rDelta,C_0*10,"0"),col="blue",lty=2,lwd=2)

print("3")
#U=0.5
f2=subset(f,U==0.5)
lines(f2$rDelta,f2$pcf-1,lty=1,col="red",lwd=2)
lines(f2$rDelta,G_theoretical(0.51,f2$rDelta,C_0*10,"0"),col="red",lty=2,lwd=2)

print("4")
#U=2.5
f3=subset(f,U==2.5)
lines(f3$rDelta,f3$pcf-1,lty=1,col="darkgreen",lwd=2)
lines(f2$rDelta,G_theoretical(2.4,f2$rDelta,C_0*10,"0"),col="darkgreen",lty=2,lwd=2)

abline(a=5,b=-2,col="black",lty=3,lwd=3)
legend("bottomleft",c('U=0','U=0.1','U=0.5','U=2.5',"Theory"),col=c("black","blue","red","darkgreen","black"),lty=c(1,1,1,1,2,2))



plot(0,0,t="n",xlim=c(0,20),ylim=range(f$pcf),xlab="r/Delta",ylab="g(r,t)-1")

#U=0
f0=subset(f,U==0)
lines(f0$rDelta,f0$pcf-1,lty=1,col="black")
lines(f0$rDelta,G_theoretical(0.0,f0$rDelta,C_0,"1"),col="black",lty=2,lwd=2)

#U=0.1
f1=subset(f,U==0.1)
lines(f1$rDelta,f1$pcf-1,lty=1,col="blue",lwd=2)
lines(f1$rDelta,G_theoretical(0.026,f1$rDelta,C_0,"0"),col="blue",lty=2,lwd=2)

#U=0.5
f2=subset(f,U==0.5)
lines(f2$rDelta,f2$pcf-1,lty=1,col="red",lwd=2)
lines(f2$rDelta,G_theoretical(0.51,f2$rDelta,C_0,"0"),col="red",lty=2,lwd=2)

#U=2.5
f2=subset(f,U==2.5)
lines(f2$rDelta,f2$pcf-1,lty=1,col="darkgreen",lwd=2)
lines(f2$rDelta,G_theoretical(2.4,f2$rDelta,C_0,"0"),col="darkgreen",lty=2,lwd=2)
dev.off()

pdf("Comparison_dx.pdf")
#Comparison with dx=10^(-9)
fm8=read.table('pcf_particle_Fig3_dx10m9.txt',sep=";",header=F,dec=".")
colnames(fm8)=c("rDelta","U","pcf")

plot(0,0,t="n",xlim=range(fm8$rDelta),ylim=c(0.1,max(f$pcf)),xlab="r/Delta",ylab="g(r,t)-1",log="xy")

f0=subset(f,U==0)
lines(f0$rDelta,f0$pcf-1,lty=1,col="black",lwd=2)
f0=subset(fm8,U==0)
points(f0$rDelta,f0$pcf-1,pch=16,col="grey",lwd=2,cex=1.5)

print("2")
#U=0.1
f1=subset(f,U==0.1)
lines(f1$rDelta,f1$pcf-1,lty=1,col="blue",lwd=2)
f1=subset(fm8,U==0.1)
points(f1$rDelta,f1$pcf-1,pch=16,col="cyan",lwd=2,cex=1.5)

print("3")
#U=0.5
f2=subset(f,U==0.5)
lines(f2$rDelta,f2$pcf-1,lty=1,col="red",lwd=2)
f2=subset(fm8,U==0.5)
points(f2$rDelta,f2$pcf-1,col="pink",lwd=2,pch=16,cex=1.5)

print("4")
#U=2.5
f3=subset(f,U==2.5)
lines(f3$rDelta,f3$pcf-1,lty=1,col="darkgreen",lwd=2)
f3=subset(fm8,U==2.5)
points(f3$rDelta,f3$pcf-1,pch=16,col="green",lwd=2,cex=1.5)

legend("bottomleft",c("dx=10^(-8)","dx=10^(-9)"),col=c("black","grey"),pch=c(NA,16),lty=c(1,NA),bty="n")
dev.off()

pdf("G_function_U0.pdf",width=10,height=5)
par(mfrow=c(1,2))
plot(0.1,0.1,t="n",xlim=range(f$rDelta),ylim=c(0.01,max(f$pcf)),xlab="r/Delta",ylab="g(r,t)-1",log="xy")

print("1")
#U=0
f0=subset(f,U==0)
lines(f0$rDelta,f0$pcf-1,lty=1,col="black",lwd=2)
plot(f0$rDelta,G_theoretical(0.0,f0$rDelta,C_0,"1"),col="black",lty=2,lwd=2,xlab="r/Delta",log="xy")
dev.off()

pdf("Comparison_area.pdf")
#Comparison with dx=10^(-9)
f=read.table('pcf_particle_Fig3_dx10m8.txt',sep=";",header=F,dec=".")
f_small=read.table('pcf_particle_Fig3_dx10m8_small_area.txt',sep=";",header=F,dec=".")
colnames(f_small)=c("rDelta","U","pcf")

plot(0,0,t="n",xlim=range(fm8$rDelta),ylim=c(0.1,max(f$pcf)),xlab="r/Delta",ylab="g(r,t)-1",log="xy")

f0=subset(f,U==0)
lines(f0$rDelta,f0$pcf-1,lty=1,col="black",lwd=2)
pcf_tmp=f0$pcf
f0=subset(f_small,U==0)
points(f0$rDelta,f0$pcf-1,pch=16,col="grey",lwd=2,cex=1.5)

print("2")
#U=0.1
f1=subset(f,U==0.1)
lines(f1$rDelta,f1$pcf-1,lty=1,col="blue",lwd=2)
f1=subset(f_small,U==0.1)
points(f1$rDelta,f1$pcf-1,pch=16,col="cyan",lwd=2,cex=1.5)

f_small=read.table('pcf_particle_Fig3_dx10m8_small_area_high_U.txt',sep=";",header=F,dec=".")
colnames(f_small)=c("rDelta","U","pcf")
print("3")
#U=0.5
f2=subset(f,U==0.5)
lines(f2$rDelta,f2$pcf-1,lty=1,col="red",lwd=2)
f2=subset(f_small,U==0.5)
points(f2$rDelta,f2$pcf-1,col="pink",lwd=2,pch=16,cex=1.5)

print("4")
U=2.5
f3=subset(f,U==2.5)
lines(f3$rDelta,f3$pcf-1,lty=1,col="darkgreen",lwd=2)
f3=subset(f_small,U==2.5)
points(f3$rDelta,f3$pcf-1,pch=16,col="green",lwd=2,cex=1.5)

legend("bottomleft",c("Area=1","Area=0.1"),col=c("black","grey"),pch=c(NA,16),lty=c(1,NA),bty="n")
dev.off()

