#23/01/2021 CP: This script uses the output from  compute_gamma to compute the distance r(t) between pairs of particles with time, and then compute gamma as the slope of (1/d)<ln(r(t))> where <ln(r(t))> is the average obtained with all the 800 pairs of particles. We do that for 4 different values of Utot. This will be used for eq(2) for Fig.3
#28/01/21: Comapring with the analytical version of Birch et al. 2007

rm(list=ls())
graphics.off()

analytical_Birch=function(U,x){ #U is the advection parameter which corresponds to U\tau/2 in Young et al. 2001, x is the multiplying factor for U to find the values we obtain with the simulation
	tau_u=U*x*2*pi #I'm taking k_m=2\pi et \tau=1
	tmp_log=log(1+tau_u^2/10+tau_u^4/67)/2
	lambda=tmp_log
	return(lambda)
}

x=seq(-2*pi,2*pi,0.1)

#List estimated coef with Young et al. 2001
coef_estimated=c()

pdf("gamma_for_different_Utot.pdf")
par(mfrow=c(2,2))

f=read.table("Spatial_distribution_particle_compute_gamma.txt",header=F,sep=";",dec='.')
colnames(f)=c("Utot","t","x","y","yfirst","type","first_parent")

Utot_list=unique(f$Utot)
#Utot_list="0.1"
for (u in 1:length(Utot_list)){
Utot_tmp=Utot_list[u]
print(paste("Utot",Utot_tmp))
timestep=unique(f$t)
unique_pair=unique(f$first_parent)

if (Utot_tmp=="0.5"){ #Stabilizes quickly: 15 time steps seems too be a good proxy
subset=1:15
}else if (Utot_tmp=="2.5"){ ##Stabilizes very quickly: 4 time steps and that's done
subset=1:4
}else{subset=1:length(timestep)}
dist=rep(0,length(subset))
dist[1]=log(10^(-7))*0.5
#Compute distance between parent and children
for (p in 1:length(unique_pair)){
#for (p in 1:1){
pair=unique_pair[p]
f_pair=subset(f,first_parent==pair&Utot==Utot_tmp)
for(t in 2:subset[length(subset)]){
	xP=f_pair$x[f_pair$t==timestep[t]&f_pair$type=="P"]
	yP=f_pair$y[f_pair$t==timestep[t]&f_pair$type=="P"]
	xC=f_pair$x[f_pair$t==timestep[t]&f_pair$type=="C"]
	yC=f_pair$y[f_pair$t==timestep[t]&f_pair$type=="C"]
	dist_tmp=sqrt((xP-xC)^2+(yP-yC)^2)
#	print(dist_tmp)
	if(dist_tmp==0.0){dist_tmp=10^(-9)}
	if(is.nan(log(dist_tmp))){stop()}
	dist[t]=dist[t]+log(dist_tmp)
}
}
dist[2:length(dist)]=0.5*dist[2:length(dist)]/length(unique_pair)

l=lm(dist~timestep[subset])
slope=l$coefficients[2]
coef_estimated=c(coef_estimated,slope)
plot(timestep[subset],dist,main=paste("U=",Utot_tmp,"\nslope=",format(slope,digits=3),sep=""),xlab="t",ylab="log(r)")

}

dev.off()

pdf("gamma_Birch_Young.pdf")
x=seq(0,2,0.01)
par(mfrow=c(2,2))
for (u in 1:length(Utot_list)){
	val_Birch=analytical_Birch(Utot_list[u],x)
	#find the closest valuei, knowing already that this is an increasing function
	first_dist=abs(val_Birch[1]-coef_estimated[u])
	j=1
	while(first_dist>=abs(val_Birch[j]-coef_estimated[u]) & j<length(x)){
		first_dist=abs(val_Birch[j]-coef_estimated[u])
		j=j+1
	}
	best_val=x[j]
	plot(x,val_Birch,main=paste("Best factor=",best_val,sep=""),xlab="x",ylab="gamma",t="l")
	abline(h=coef_estimated[u],col="red")
	if(u==1){
	legend("bottomright",c("Birch et al. 2007","Young et al. 2001"),col=c("black","red"),lty=1,bty="n")
	}else{
	abline(v=best_val)
	}
}
dev.off()
