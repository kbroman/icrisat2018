# congenic line

bgcolor <- rgb(0, 0, 80, maxColorValue=255)
color <- broman::brocolors("crayons")[c("Cornflower", "Blush")]

library(qtl)
data(map10)
L <- summary(map10)[1:20,2]
#L <- read.csv("/Users/kbroman/Projects/MouseMaps/lengths.csv")[,4]
L[L<51] <- 51
L[21] <- 10

pdf("../Figs/congenic.pdf", width=9.75, height=6.5,
    pointsize=16, onefile=TRUE)

par(mar=rep(0.1,4),las=1,fg="white",col="white",col.axis="white",col.lab="white",
    bg=bgcolor, bty="n")
plot(0,0,xlim=c(0,864),ylim=c(0,480),xaxt="n",yaxt="n",xlab="",ylab="",type="n")

u <- par("usr")
xcen <- seq(u[1]+3,u[2]-3,len=22)
xcen <- apply(cbind(xcen[-1],xcen[-22]),1,mean)

text(xcen,475,c(1:19,"X","Y"),cex=1.2,adj=c(0.5,1))

mult <- 430/max(L)
for(i in 1:21) {
    rect(xcen[i]-10,430,xcen[i]+10,430-L[i]*mult,col=color[1],border=color[1], lend=1, ljoin=1)
}

rect(xcen[2]-10,430-37*mult, xcen[2]+10,430-(37+20)*mult,
     col=color[2], border=color[2],
     lend=1, ljoin=1)


dev.off()
