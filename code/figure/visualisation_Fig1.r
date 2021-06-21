#23/01/2021 CP: This script uses the output from main_Fig1.cpp to reproduce Fig1 in Young et al. 2001

rm(list=ls())
graphics.off()

f=read.table("../simulation/Spatial_distribution_particle_Fig1.txt",header=F,sep=";",dec='.')
colnames(f)=c("t","x","y","yfirst","first_parent")
 
png("spatial_distribution_Fig1.png",width=400,height=1150)
par(mfrow=c(3,1),cex=1.5,mar=c(4.,4,0.,0.2),oma=c(0.,.5,3.,1))

##Create color and parent correspondance, to keep the first y position. The color of a particle is a function of the y position of the first parent.
select_t=0

yfirst=f$yfirst[f$t==select_t]
first_parent=f$first_parent[f$t==select_t]
table_parent=matrix(NA,length(yfirst),3) #Y position, name of the parent, color
order_yfirst = order(yfirst)
cols=rainbow(length(yfirst))

table_parent[,1]=yfirst[order_yfirst]
table_parent[,2]=first_parent[order_yfirst]
table_parent[,3]=cols

#First plot, at t=0
print(select_t)
x=f$x[f$t==select_t]
print(length(x))
y=f$y[f$t==select_t]
plot(0,0,xlim=c(0,1),ylim=c(0,1),t="n",xlab="",ylab="y",xaxt="n")
for(i in 1:length(y)){ #Plotting one point after the other
        find_col=table_parent[table_parent[,2]==first_parent[i],3] #Identifying the ancestor of the point, and the corresponding color of such parent
        points(x[i],y[i],pch=19,col=find_col,cex=0.1)
}
text(-0.25,1.0,"a",xpd=T,font=2)

#Now, we want the spatial distribution of each point at t=100
par(mar=c(4.,4.,0.,0.2))
select_t=100
print(select_t)
x=f$x[f$t==select_t]
y=f$y[f$t==select_t]
yfirst=f$yfirst[f$t==select_t]
first_parent=f$first_parent[f$t==select_t]

print(length(x))
plot(0,0,xlim=c(0,1),ylim=c(0,1),t="n",xlab="",ylab="y",xaxt="n")
for(i in 1:length(y)){
        find_col=table_parent[table_parent[,2]==first_parent[i],3] #Identifying the ancestor of the point, and the corresponding color of such parent
        points(x[i],y[i],pch=19,col=find_col,cex=0.1)
}
text(-0.25,1.0,"b",xpd=T,font=2)

#Now, we want the spatial distribution of each point at t=1000
par(mar=c(4.,4.,0.,0.2))
select_t=max(f$t)
print(select_t)
x=f$x[f$t==select_t]
y=f$y[f$t==select_t]
yfirst=f$yfirst[f$t==select_t]
first_parent=f$first_parent[f$t==select_t]

print(length(x))

plot(0,0,xlim=c(0,1),ylim=c(0,1),t="n",xlab="x",ylab="y")
for(i in 1:length(y)){
	find_col=table_parent[table_parent[,2]==first_parent[i],3] #Identifying the ancestor of the point, and the corresponding color of such parent
	points(x[i],y[i],pch=19,col=find_col,cex=0.1)
}
text(-0.25,1.0,"c",xpd=T,font=2)
dev.off()
