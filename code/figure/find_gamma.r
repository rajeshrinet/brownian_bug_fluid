#23/01/2021 CP: This script uses the output from  compute_gamma to compute the distance r(t) between pairs of particles with time, and then compute gamma as the slope of (1/d)<ln(r(t))> where <ln(r(t))> is the average obtained with all the 800 pairs of particles. We do that for 4 different values of Utot. This will be used for eq(2) for Fig.3

rm(list=ls())
graphics.off()

x=seq(-2*pi,2*pi,0.1)

#List of estimated coefficients
coef_estimated=c()

pdf("gamma_for_different_Utot.pdf")
par(mfrow=c(2,2))

f=read.table("../simulation/Spatial_distribution_particle_compute_gamma.txt",header=F,sep=";",dec='.')
colnames(f)=c("Utot","t","x","y","yfirst","type","first_parent")

Utot_list=unique(f$Utot)

#Each subset corresponds to a value of U_tot
for (u in 1:length(Utot_list)){
	Utot_tmp=Utot_list[u]
	print(paste("Utot",Utot_tmp))
	timestep=unique(f$t)
	unique_pair=unique(f$first_parent)

#To compute the linear coefficient, we need to only examine the separation as a function of time before steady state. We therefore only use a subset of timestep
	if (Utot_tmp=="0.5"){ #Stabilizes quickly: 15 time steps seems too be a good proxy
		subset=1:15
	}else if (Utot_tmp=="2.5"){ ##Stabilizes very quickly: 4 time steps and stable state
		subset=1:4
	}else{subset=1:length(timestep)}

#Compute distance between parent and children
	dist=rep(0,length(subset))
	dist[1]=log(10^(-7))*0.5 #The first distance between particle is  10-7
	for (p in 1:length(unique_pair)){ #A pair is defined by a unique id, pair, for particles identified as "P" or "C"
		pair=unique_pair[p] 
		f_pair=subset(f,first_parent==pair&Utot==Utot_tmp)
		for(t in 2:subset[length(subset)]){ #For each time step, we compute the mean distance between particles
			xP=f_pair$x[f_pair$t==timestep[t]&f_pair$type=="P"]
			yP=f_pair$y[f_pair$t==timestep[t]&f_pair$type=="P"]
			xC=f_pair$x[f_pair$t==timestep[t]&f_pair$type=="C"]
			yC=f_pair$y[f_pair$t==timestep[t]&f_pair$type=="C"]
			dist_tmp=sqrt((xP-xC)^2+(yP-yC)^2)
			if(dist_tmp==0.0){dist_tmp=10^(-9)} #If the distance is too small, we set an arbitrary low value
			if(is.nan(log(dist_tmp))){stop()}
			dist[t]=dist[t]+log(dist_tmp) #We sum all distances and will then divide by the number of particles
		}
	}
	dist[2:length(dist)]=0.5*dist[2:length(dist)]/length(unique_pair)

	l=lm(dist~timestep[subset]) #Finally, we estimate the value of the slope
	slope=l$coefficients[2]
	coef_estimated=c(coef_estimated,slope)
	plot(timestep[subset],dist,main=eval(substitute(expression(paste("U",tau,"/2=",Utot_tmp,", ",gamma,"=",slope,sep="")),list(Utot_tmp=Utot_tmp,slope=format(slope,digits=3)))),xlab="t",ylab="log(r)")

}

dev.off()
