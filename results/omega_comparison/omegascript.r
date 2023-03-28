####please provide the command as Rscript omegascript.r filewithLPS filewithoutLPS purity
args = commandArgs(trailingOnly=TRUE)
f1=read.table(args[1],header=F)
f2=read.table(args[2],header=F)
pur=as.character(args[3])
f1=f1[,c(1,3)]
f2=f2[,c(1,3)]
f1$stat="With_LPSs"
f2$stat="Without_LPSs"
data=rbind(f1,f2)
colnames(data)=c("clades","omega","status")
data_new <- data
data_new$clades <- factor(data_new$clades,c("Afrotheria","Primates","Rodentia","Chiroptera","Artiodactyla","Perissodactyla","Carnivora","Marsupialia","Testudines","Aves","Squamata","Amphibia"))


##new_order <- with(data, reorder(clades , omega, mean , na.rm=T))
image_name=paste('omega_cladewise_LPSs_nonLPSs_',pur,'_comparison.jpeg',sep="")
jpeg(image_name,width=28,height=18,units="in",res=300)
par(mar=c(5,4,4,1))
myplot <- boxplot(omega ~ status*clades, data=data_new,ylim=c(-0.75,2), 
        boxwex=0.4,outline=F,
        ann=F,
        col=c("seagreen4" , "coral3") , 
        border=c("seagreen3","coral2"),cex.axis=1.25,
        xaxt="n")
maintext=paste("Omega variation in proteins with LPSs and proteins without LPSs at ",pur," purity cutoff",sep="")
title(main=maintext ,line=1,cex.main=2.3)
title(ylab="Omega (dN/dS)", line=2, cex.lab=2.5)        
title(xlab="Clades", line=3, cex.lab=2.5)
my_names <- sapply(strsplit(myplot$names , '\\.') , function(x) x[[2]] )
my_names <- my_names[seq(1 , length(my_names) , 2)]
##change the number in seq below according to your number of clades. It should be 2*(number of clades)
axis(1, 
     at = seq(1.5 , 24 , 2), 
     labels = my_names , 
     tick=FALSE , cex.axis=1.75)
     
# Add the grey vertical lines
##change the seq range (in our case it is 30) below according to your number of clades 
for(i in seq(0.5 , 30 , 2)){ 
  abline(v=i,lty=3, col="grey",lwd=1.75)
  }
 
# Add a legend ##Adjust the position of legend by changing x and y below
legend(x=21,y=-0.6, legend = c("Proteins with LPSs", "Proteins without LPSs"), 
       col=c("seagreen4" , "coral3"),
       pch = 15, pt.cex = 6, cex = 2,  horiz = F, inset = c(0.1, 0.1),box.lwd = 0,box.col = "white",bg = "white")

###adding the significant values for each clade paired comparison
##in yax, I am adding +4 to keep it slightly above the 3rd quantile range of boxplot, you can adjust it according to your plot
xpos=1
for(clade in unique(my_names)){
a1=data[data$clades==clade & data$status=="With_LPSs",]
b1=data[data$clades==clade & data$status=="Without_LPSs",]
b1=b1[!is.na(as.numeric(b1$omega)), ]
a1=a1[!is.na(as.numeric(a1$omega)), ]
yax=max(quantile(a1$omega,probs=0.75,na.rm=T),quantile(b1$omega,probs=0.75,na.rm=T))+0.2
pval=signif(as.numeric(wilcox.test(a1$omega,b1$omega,alternative="less")[3]),digits=3)
segments(xpos,yax,xpos+1,yax,lty=3,col="purple")
text(xpos+0.5,yax+0.6,labels=pval,cex=1.5)
a1mean=round(mean(a1$omega),digits=2)
b1mean=round(mean(b1$omega),digits=2)
yminax=min(quantile(a1$omega,probs=0.25,na.rm=T),quantile(b1$omega,probs=0.25,na.rm=T))-0.2
text(xpos,yminax,bquote(bar("x") == .(a1mean)),cex=1.1)
text(xpos+1,yminax,bquote(bar("x") == .(b1mean)),cex=1.1)
text(xpos,yminax-0.25,bquote("n" == .(dim(a1)[1])),cex=1.1)
text(xpos+1,yminax-0.25,bquote("n" == .(dim(b1)[1])),cex=1.1)
xpos=xpos+2
}
dev.off()

data_new <- data
data_new=data_new[data_new$omega<=2,]
data_new$clades <- factor(data_new$clades,c("Afrotheria","Primates","Rodentia","Chiroptera","Artiodactyla","Perissodactyla","Carnivora","Marsupialia","Testudines","Aves","Squamata","Amphibia"))


##new_order <- with(data, reorder(clades , omega, mean , na.rm=T))
image_name=paste('omega_two_cutoff_cladewise_LPSs_nonLPSs_',pur,'_comparison.jpeg',sep="")
jpeg(image_name,width=28,height=18,units="in",res=300)
par(mar=c(5,4,4,1))
myplot <- boxplot(omega ~ status*clades, data=data_new,ylim=c(-0.75,2), 
        boxwex=0.4,outline=F,
        ann=F,
        col=c("seagreen4" , "coral3") , 
        border=c("seagreen3","coral2"),cex.axis=1.25,
        xaxt="n")
maintext=paste("Omega (<=2) variation in proteins with LPSs and proteins without LPSs at ",pur," purity cutoff",sep="")
title(main=maintext ,line=1,cex.main=2.3)
title(ylab="Omega (dN/dS)", line=2, cex.lab=2.5)        
title(xlab="Clades", line=3, cex.lab=2.5)
my_names <- sapply(strsplit(myplot$names , '\\.') , function(x) x[[2]] )
my_names <- my_names[seq(1 , length(my_names) , 2)]
##change the number in seq below according to your number of clades. It should be 2*(number of clades)
axis(1, 
     at = seq(1.5 , 24 , 2), 
     labels = my_names , 
     tick=FALSE , cex.axis=1.75)
     
# Add the grey vertical lines
##change the seq range (in our case it is 30) below according to your number of clades 
for(i in seq(0.5 , 30 , 2)){ 
  abline(v=i,lty=3, col="grey",lwd=1.75)
  }

# Add a legend ##Adjust the position of legend by changing x and y below
legend(x=21,y=-0.6, legend = c("Proteins with LPSs", "Proteins without LPSs"), 
       col=c("seagreen4" , "coral3"),
       pch = 15, pt.cex = 6, cex = 2,  horiz = F, inset = c(0.1, 0.1),box.lwd = 0,box.col = "white",bg = "white")

###adding the significant values for each clade paired comparison
##in yax, I am adding +4 to keep it slightly above the 3rd quantile range of boxplot, you can adjust it according to your plot
xpos=1
for(clade in unique(my_names)){
a1=data[data$clades==clade & data$status=="With_LPSs",]
b1=data[data$clades==clade & data$status=="Without_LPSs",]
b1=b1[!is.na(as.numeric(b1$omega)), ]
a1=a1[!is.na(as.numeric(a1$omega)), ]
yax=max(quantile(a1$omega,probs=0.75,na.rm=T),quantile(b1$omega,probs=0.75,na.rm=T))+0.2
pval=signif(as.numeric(wilcox.test(a1$omega,b1$omega,alternative="less")[3]),digits=3)
segments(xpos,yax,xpos+1,yax,lty=3,col="purple")
text(xpos+0.5,yax+0.6,labels=pval,cex=1.5)
a1mean=round(mean(a1$omega),digits=2)
b1mean=round(mean(b1$omega),digits=2)
yminax=min(quantile(a1$omega,probs=0.25,na.rm=T),quantile(b1$omega,probs=0.25,na.rm=T))-0.2
text(xpos,yminax,bquote(bar("x") == .(a1mean)),cex=1.1)
text(xpos+1,yminax,bquote(bar("x") == .(b1mean)),cex=1.1)
text(xpos,yminax-0.25,bquote("n" == .(dim(a1)[1])),cex=1.1)
text(xpos+1,yminax-0.25,bquote("n" == .(dim(b1)[1])),cex=1.1)
xpos=xpos+2
}
dev.off()
######################Rscript ends here####
