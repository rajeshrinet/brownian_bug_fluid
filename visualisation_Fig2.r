#23/01/2021 CP: This script uses the output from main_Fig2a.cpp and main_Fig2b.cpp to reproduce Fig2 in Young et al. 2001

rm(list=ls())
graphics.off()

f=read.table("Spatial_distribution_particle_Fig2a.txt",header=F,sep=";",dec='.')
colnames(f)=c("t","x","y","yfirst","first_parent")
 
png("spatial_distribution_Fig2.png",width=450,height=900)
par(mfrow=c(2,1),cex=1.5,oma=c(0.,.5,.25,0.25))

##Create color and parent correspondance, to keep the first y position
#pal = colorRampPalette(c("blue", "red"))
select_t=0

yfirst=f$yfirst[f$t==select_t]
first_parent=f$first_parent[f$t==select_t]
table_parent=matrix(NA,length(yfirst),3) #Y, name of the parent, color
order_yfirst = order(yfirst)
#cols=pal(length(yfirst))
cols=rainbow(length(yfirst))

table_parent[,1]=yfirst[order_yfirst]
table_parent[,2]=first_parent[order_yfirst]
table_parent[,3]=cols

#Now, we want the spatial distribution of each point at t=select_t
par(mar=c(2.,4.,2.,0.2))
select_t=30
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
text(-0.25,1.0,"a",xpd=T,font=2)

f=read.table("Spatial_distribution_particle_Fig2b.txt",header=F,sep=";",dec='.')
colnames(f)=c("t","x","y","yfirst","first_parent")
##Create color and parent correspondance, to keep the first y position
#pal = colorRampPalette(c("blue", "red"))
select_t=0

yfirst=f$yfirst[f$t==select_t]
first_parent=f$first_parent[f$t==select_t]
table_parent=matrix(NA,length(yfirst),3) #Y, name of the parent, color
order_yfirst = order(yfirst)
#cols=pal(length(yfirst))
cols=rainbow(length(yfirst))

table_parent[,1]=yfirst[order_yfirst]
table_parent[,2]=first_parent[order_yfirst]
table_parent[,3]=cols

#Now, we want the spatial distribution of each point at t=select_t
par(mar=c(2.,4.,2.,0.2))
select_t=1000
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
dev.off()
