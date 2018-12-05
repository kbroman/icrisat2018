# figure 1: reproduce Fig 5 from Gatti et al (2014)
#
# "We regressed log neutrophil counts on founder allele dosages at
# each marker using a kinship correction with sex and log white blood cell counts as covariates"

# taken from figure for R/qtl2 paper
file <- "cache/maps_n_phe.RData"
load(file)

# genome scan with full model
file <- "cache/out_full.rds"
out_full <- readRDS(file)

# genome scan with additive model
file <- "cache/out_add.rds"
out_add <- readRDS(file)

# GWAS at SNPs
file <- "cache/out_snps.rds"
out_snps <- readRDS(file)


# load permutation results
operm_full <- readRDS("Perms/operm_full.rds")
operm_add <- readRDS("Perms/operm_add.rds")
operm_snps <- readRDS("Perms/operm_snps.rds")

# calculate thresholds
thr_full <- summary(operm_full)
thr_add <- summary(operm_add)
thr_snps <- summary(operm_snps)

# ylim
ymx_full <- thr_full$A/thr_add$A*maxlod(out_add)*1.04
ymx_add <- maxlod(out_add)*1.04
ymx_snps <- thr_snps$A/thr_add$A*maxlod(out_add)*1.04

# make the plots
altcolor <- "green4"
linecolor <- "violetred"
panel_lab_adj <- c(0.12, 0.06)
panel_lab_cex <- 1.3

res <- 256
png("../Figs/do_scan.png", height=7.5*res, width=10*res, pointsize=14, res=res)
par(mfrow=c(3,1))
par(mar=c(2.1, 4.1, 3.6, 1.1))
bgcolor <- rgb(0, 0, 96, maxColorValue=255)
par(bg=bgcolor, fg="white", col="white",col.axis="white",col.lab="white", col.main="lightblue", cex.main=2)

plot(out_full, pmap, xlab="", ylim=c(0, ymx_full), altcol=altcolor, main="Full model")
u <- par("usr")
endA <- xpos_scan1(pmap, thechr=19, thepos=max(pmap[[19]]))+25/2
segments(u[1], thr_full$A, endA, thr_full$A, col=linecolor, lty=2)
segments(endA, thr_full$X, u[2], thr_full$X, col=linecolor, lty=2)
u <- par("usr")
text(u[1]-diff(u[1:2])*panel_lab_adj[1], u[4]+diff(u[3:4])*panel_lab_adj[2], "A", font=2, xpd=TRUE, cex=panel_lab_cex)

plot(out_add, pmap, xlab="", ylim=c(0, ymx_add), altcol=altcolor, main="Additive alleles")
segments(u[1], thr_add$A, endA, thr_add$A, col=linecolor, lty=2)
segments(endA, thr_add$X, u[2], thr_add$X, col=linecolor, lty=2)
u <- par("usr")
text(u[1]-diff(u[1:2])*panel_lab_adj[1], u[4]+diff(u[3:4])*panel_lab_adj[2], "B", font=2, xpd=TRUE, cex=panel_lab_cex)

plot(out_snps$lod, out_snps$snpinfo, altcol=altcolor, xlab="",
     ylim=c(0, ymx_snps), main="SNP associations")
segments(u[1], thr_snps$A, endA, thr_snps$A, col=linecolor, lty=2)
segments(endA, thr_snps$X, u[2], thr_snps$X, col=linecolor, lty=2)
u <- par("usr")
text(u[1]-diff(u[1:2])*panel_lab_adj[1], u[4]+diff(u[3:4])*panel_lab_adj[2], "C", font=2, xpd=TRUE, cex=panel_lab_cex)

dev.off()
