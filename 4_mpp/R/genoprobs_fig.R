# illustration of DO genotype reconstruction

blue_background <- TRUE # make TRUE if you want a blue background

# libraries
library(qtl2)     # install with install.packages("qtl2", repos="https://rqtl.org/qtl2cran")
library(broman)
blue <- rgb(0, 0, 80, maxColorValue=255)
gray90 <- ifelse(blue_background, "gray40", "gray90")

cache_dir <- "_cache"
if(!dir.exists(cache_dir)) dir.create(cache_dir)
file <- file.path(cache_dir, "do_geno.RData")
if(file.exists(file)) {
    load(file)
} else {
    doex_file <- "https://raw.githubusercontent.com/rqtl/qtl2data/master/DOex/DOex.zip"
    DOex <- read_cross2(doex_file)

    markers <- marker_names(DOex[,"2"])[79:93]
    DOex <- pull_markers(DOex[c(2, 1, 9),"2"], markers)

    g <- DOex$geno[[1]]
    fg <- DOex$founder_geno[[1]]
    pr <- calc_genoprob(DOex, err=0.02)
    pmap <- DOex$pmap[[1]] - min(DOex$pmap[[1]])

    save(g, fg, pr, pmap, file=file)
}


### version A: just founders + mouse HS-1
pdf("../Figs/genoprobsA.pdf", width=9.75, height=6.5, pointsize=16)

left_gap <- 2
right_gap <- 0.5
point_cex <- 1.7

par(mar=rep(0.1,4), bty="n")
if(blue_background) {
    par(fg="white",col="white",col.axis="white",col.lab="white",bg=blue)
} else {
    par(fg=blue,col=blue,col.axis=blue,col.lab=blue,bg="white")
}

plot(0, 0, type="n", xlim=range(pmap)+c(-left_gap, right_gap), ylim=c(0, 100), xlab="", ylab="", yaxt="n", xaxt="n",
     yaxs="i", xaxs="i")

# founder genotypes
fg_y <- rev(seq(68, 96, length=nrow(fg)))

# labels
text(-left_gap*2/3, mean(fg_y), "Founders", srt=90, adj=0.5, cex=1.2)
text(-left_gap/5, fg_y, LETTERS[1:8], adj=c(1,0.5), col=qtl2::CCcolors)

# dots at genotypes
for(i in 1:nrow(fg)) {
    points(pmap, rep(fg_y[i], length(pmap)), pch=21, bg=c("white", "gray", blue)[fg[i,]],
           cex=point_cex)
}

# DO mouse 1: hom with no crossovers
g_y <- c(55, 49.5, 46.5)
text(-left_gap/5, g_y[1], "HS-1", adj=c(1,0.5))
points(pmap, rep(g_y[1], length(pmap)), pch=21, bg=c("white", "gray", blue)[g[1,]], cex=point_cex)

dev.off()



### version B: add inferred genotypes for mouse HS-1
pdf("../Figs/genoprobsB.pdf", width=9.75, height=6.5, pointsize=16)

left_gap <- 2
right_gap <- 0.5
point_cex <- 1.7

par(mar=rep(0.1,4), bty="n")
if(blue_background) {
    par(fg="white",col="white",col.axis="white",col.lab="white",bg=blue)
} else {
    par(fg=blue,col=blue,col.axis=blue,col.lab=blue,bg="white")
}

plot(0, 0, type="n", xlim=range(pmap)+c(-left_gap, right_gap), ylim=c(0, 100), xlab="", ylab="", yaxt="n", xaxt="n",
     yaxs="i", xaxs="i")

# founder genotypes
fg_y <- rev(seq(68, 96, length=nrow(fg)))

u <- par("usr")
dx <- diff(u[1:2])*0.012
dy <- diff(u[3:4])*0.012/6.5*9.75
rect(min(pmap)-dx, fg_y[7]-dy, max(pmap)+dx, fg_y[7]+dy, col=gray90, border=blue, lend=1, ljoin=1)

# labels
text(-left_gap*2/3, mean(fg_y), "Founders", srt=90, adj=0.5, cex=1.2)
text(-left_gap/5, fg_y, LETTERS[1:8], adj=c(1,0.5), col=qtl2::CCcolors)

# dots at genotypes
for(i in 1:nrow(fg)) {
    points(pmap, rep(fg_y[i], length(pmap)), pch=21, bg=c("white", "gray", blue)[fg[i,]],
           cex=point_cex)
}

# DO mouse 1: hom with no crossovers
g_y <- c(55, 49.5, 46.5)
text(-left_gap/5, g_y[1], "HS-1", adj=c(1,0.5))
points(pmap, rep(g_y[1], length(pmap)), pch=21, bg=c("white", "gray", blue)[g[1,]], cex=point_cex)

for(i in 1:2) {
    points(pmap, rep(g_y[i+1], length(pmap)), pch=21, bg=qtl2::CCcolors[7], cex=point_cex)
}

dev.off()



### version C: mouse HS-2
pdf("../Figs/genoprobsC.pdf", width=9.75, height=6.5, pointsize=16)

left_gap <- 2
right_gap <- 0.5
point_cex <- 1.7

par(mar=rep(0.1,4), bty="n")
if(blue_background) {
    par(fg="white",col="white",col.axis="white",col.lab="white",bg=blue)
} else {
    par(fg=blue,col=blue,col.axis=blue,col.lab=blue,bg="white")
}

plot(0, 0, type="n", xlim=range(pmap)+c(-left_gap, right_gap), ylim=c(0, 100), xlab="", ylab="", yaxt="n", xaxt="n",
     yaxs="i", xaxs="i")

# founder genotypes
fg_y <- rev(seq(68, 96, length=nrow(fg)))

u <- par("usr")
dx <- diff(u[1:2])*0.012
dy <- diff(u[3:4])*0.012/6.5*9.75

# labels
text(-left_gap*2/3, mean(fg_y), "Founders", srt=90, adj=0.5, cex=1.2)
text(-left_gap/5, fg_y, LETTERS[1:8], adj=c(1,0.5), col=qtl2::CCcolors)

# dots at genotypes
for(i in 1:nrow(fg)) {
    points(pmap, rep(fg_y[i], length(pmap)), pch=21, bg=c("white", "gray", blue)[fg[i,]],
           cex=point_cex)
}

# DO mouse 1: hom with no crossovers
g_y <- c(55, 49.5, 46.5)
text(-left_gap/5, g_y[1], "HS-1", adj=c(1,0.5))
points(pmap, rep(g_y[1], length(pmap)), pch=21, bg=c("white", "gray", blue)[g[1,]], cex=point_cex)

for(i in 1:2) {
    points(pmap, rep(g_y[i+1], length(pmap)), pch=21, bg=qtl2::CCcolors[7], cex=point_cex)
}


# DO mouse 2: het with no crossovers
g2_y <- c(34, 28.5, 25.5)
text(-left_gap/5, g2_y[1], "HS-2", adj=c(1,0.5))
points(pmap, rep(g2_y[1], length(pmap)), pch=21, bg=c("white", "gray", blue)[g[2,]], cex=point_cex)

dev.off()



### version D: inferred genotypes for mouse HS-2
pdf("../Figs/genoprobsD.pdf", width=9.75, height=6.5, pointsize=16)

left_gap <- 2
right_gap <- 0.5
point_cex <- 1.7

par(mar=rep(0.1,4), bty="n")
if(blue_background) {
    par(fg="white",col="white",col.axis="white",col.lab="white",bg=blue)
} else {
    par(fg=blue,col=blue,col.axis=blue,col.lab=blue,bg="white")
}

plot(0, 0, type="n", xlim=range(pmap)+c(-left_gap, right_gap), ylim=c(0, 100), xlab="", ylab="", yaxt="n", xaxt="n",
     yaxs="i", xaxs="i")

# founder genotypes
fg_y <- rev(seq(68, 96, length=nrow(fg)))

u <- par("usr")
dx <- diff(u[1:2])*0.012
dy <- diff(u[3:4])*0.012/6.5*9.75
rect(min(pmap)-dx, fg_y[1]-dy, max(pmap)+dx, fg_y[1]+dy, col=gray90, border=blue, lend=1, ljoin=1)
rect(min(pmap)-dx, fg_y[3]-dy, max(pmap)+dx, fg_y[3]+dy, col=gray90, border=blue, lend=1, ljoin=1)

# labels
text(-left_gap*2/3, mean(fg_y), "Founders", srt=90, adj=0.5, cex=1.2)
text(-left_gap/5, fg_y, LETTERS[1:8], adj=c(1,0.5), col=qtl2::CCcolors)

# dots at genotypes
for(i in 1:nrow(fg)) {
    points(pmap, rep(fg_y[i], length(pmap)), pch=21, bg=c("white", "gray", blue)[fg[i,]],
           cex=point_cex)
}

# DO mouse 1: hom with no crossovers
g_y <- c(55, 49.5, 46.5)
text(-left_gap/5, g_y[1], "HS-1", adj=c(1,0.5))
points(pmap, rep(g_y[1], length(pmap)), pch=21, bg=c("white", "gray", blue)[g[1,]], cex=point_cex)

for(i in 1:2) {
    points(pmap, rep(g_y[i+1], length(pmap)), pch=21, bg=qtl2::CCcolors[7], cex=point_cex)
}


# DO mouse 2: het with no crossovers
g2_y <- c(34, 28.5, 25.5)
text(-left_gap/5, g2_y[1], "HS-2", adj=c(1,0.5))
points(pmap, rep(g2_y[1], length(pmap)), pch=21, bg=c("white", "gray", blue)[g[2,]], cex=point_cex)

for(i in 1:2) {
    points(pmap, rep(g2_y[i+1], length(pmap)), pch=21, bg=qtl2::CCcolors[c(1,3)[i]], cex=point_cex)
}


dev.off()



### version E: add mouse HS-3
pdf("../Figs/genoprobsE.pdf", width=9.75, height=6.5, pointsize=16)

left_gap <- 2
right_gap <- 0.5
point_cex <- 1.7

par(mar=rep(0.1,4), bty="n")
if(blue_background) {
    par(fg="white",col="white",col.axis="white",col.lab="white",bg=blue)
} else {
    par(fg=blue,col=blue,col.axis=blue,col.lab=blue,bg="white")
}

plot(0, 0, type="n", xlim=range(pmap)+c(-left_gap, right_gap), ylim=c(0, 100), xlab="", ylab="", yaxt="n", xaxt="n",
     yaxs="i", xaxs="i")

# founder genotypes
fg_y <- rev(seq(68, 96, length=nrow(fg)))

u <- par("usr")
dx <- diff(u[1:2])*0.012
dy <- diff(u[3:4])*0.012/6.5*9.75

# labels
text(-left_gap*2/3, mean(fg_y), "Founders", srt=90, adj=0.5, cex=1.2)
text(-left_gap/5, fg_y, LETTERS[1:8], adj=c(1,0.5), col=qtl2::CCcolors)

# dots at genotypes
for(i in 1:nrow(fg)) {
    points(pmap, rep(fg_y[i], length(pmap)), pch=21, bg=c("white", "gray", blue)[fg[i,]],
           cex=point_cex)
}

# DO mouse 1: hom with no crossovers
g_y <- c(55, 49.5, 46.5)
text(-left_gap/5, g_y[1], "HS-1", adj=c(1,0.5))
points(pmap, rep(g_y[1], length(pmap)), pch=21, bg=c("white", "gray", blue)[g[1,]], cex=point_cex)

for(i in 1:2) {
    points(pmap, rep(g_y[i+1], length(pmap)), pch=21, bg=qtl2::CCcolors[7], cex=point_cex)
}


# DO mouse 2: het with no crossovers
g2_y <- c(34, 28.5, 25.5)
text(-left_gap/5, g2_y[1], "HS-2", adj=c(1,0.5))
points(pmap, rep(g2_y[1], length(pmap)), pch=21, bg=c("white", "gray", blue)[g[2,]], cex=point_cex)

for(i in 1:2) {
    points(pmap, rep(g2_y[i+1], length(pmap)), pch=21, bg=qtl2::CCcolors[c(1,3)[i]], cex=point_cex)
}


# DO mouse 3: het with one exchange
g3_y <- c(13, 7.5, 4.5)
text(-left_gap/5, g3_y[1], "HS-3", adj=c(1,0.5))
points(pmap, rep(g3_y[1], length(pmap)), pch=21, bg=c("white", "gray", blue)[g[3,]], cex=point_cex)

dev.off()



### version F: inferred genotypes for mouse HS-3
pdf("../Figs/genoprobsF.pdf", width=9.75, height=6.5, pointsize=16)

left_gap <- 2
right_gap <- 0.5
point_cex <- 1.7

par(mar=rep(0.1,4), bty="n")
if(blue_background) {
    par(fg="white",col="white",col.axis="white",col.lab="white",bg=blue)
} else {
    par(fg=blue,col=blue,col.axis=blue,col.lab=blue,bg="white")
}

plot(0, 0, type="n", xlim=range(pmap)+c(-left_gap, right_gap), ylim=c(0, 100), xlab="", ylab="", yaxt="n", xaxt="n",
     yaxs="i", xaxs="i")

# founder genotypes
fg_y <- rev(seq(68, 96, length=nrow(fg)))

u <- par("usr")
dx <- diff(u[1:2])*0.012
dy <- diff(u[3:4])*0.012/6.5*9.75
rect(min(pmap)-dx, fg_y[6]-dy, max(pmap)+dx, fg_y[6]+dy, col=gray90, border=blue, lend=1, ljoin=1)
rect(min(pmap)-dx, fg_y[5]-dy, pmap[6]+dx, fg_y[5]+dy, col=gray90, border=blue, lend=1, ljoin=1)
rect(pmap[7]-dx, fg_y[7]-dy, max(pmap)+dx, fg_y[7]+dy, col=gray90, border=blue, lend=1, ljoin=1)

# labels
text(-left_gap*2/3, mean(fg_y), "Founders", srt=90, adj=0.5, cex=1.2)
text(-left_gap/5, fg_y, LETTERS[1:8], adj=c(1,0.5), col=qtl2::CCcolors)

# dots at genotypes
for(i in 1:nrow(fg)) {
    points(pmap, rep(fg_y[i], length(pmap)), pch=21, bg=c("white", "gray", blue)[fg[i,]],
           cex=point_cex)
}

# DO mouse 1: hom with no crossovers
g_y <- c(55, 49.5, 46.5)
text(-left_gap/5, g_y[1], "HS-1", adj=c(1,0.5))
points(pmap, rep(g_y[1], length(pmap)), pch=21, bg=c("white", "gray", blue)[g[1,]], cex=point_cex)

for(i in 1:2) {
    points(pmap, rep(g_y[i+1], length(pmap)), pch=21, bg=qtl2::CCcolors[7], cex=point_cex)
}


# DO mouse 2: het with no crossovers
g2_y <- c(34, 28.5, 25.5)
text(-left_gap/5, g2_y[1], "HS-2", adj=c(1,0.5))
points(pmap, rep(g2_y[1], length(pmap)), pch=21, bg=c("white", "gray", blue)[g[2,]], cex=point_cex)

for(i in 1:2) {
    points(pmap, rep(g2_y[i+1], length(pmap)), pch=21, bg=qtl2::CCcolors[c(1,3)[i]], cex=point_cex)
}


# DO mouse 3: het with one exchange
g3_y <- c(13, 7.5, 4.5)
text(-left_gap/5, g3_y[1], "HS-3", adj=c(1,0.5))
points(pmap, rep(g3_y[1], length(pmap)), pch=21, bg=c("white", "gray", blue)[g[3,]], cex=point_cex)

# first 6 is EF ?
# last  9 is FG ?
points(pmap, rep(g3_y[2], length(pmap)), pch=21, bg=qtl2::CCcolors[6], cex=point_cex)
points(pmap, rep(g3_y[3], length(pmap)), pch=21, bg=qtl2::CCcolors[rep(c(5,7),c(6,9))], cex=point_cex)

dev.off()
