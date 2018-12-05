# Intro to QTL mapping with R/qtl
# Live coding at ICRISAT workshop 2018-12-05
#
# largely following http://rqtl.org/tutorials/rqtltour2.pdf
# latest version of this script at https://bit.ly/introqtl2018

# load R/qtl
library(qtl)

# import some QTL mapping data
# example data file at http://rqtl.org/sug.csv
sug <- read.cross("csv", "http://rqtl.org", "sug.csv",
                  genotypes=c("CC", "CB", "BB"), 
                  alleles=c("C", "B"))

# summary of data
summary(sug)
plot(sug)
plotMissing(sug) # black pixels = missing genotypes
plotMap(sug)

nind(sug)
nmar(sug)
totmar(sug)

# histograms of phenotypes
plotPheno(sug, pheno.col=1)
plotPheno(sug, pheno.col="bw")

# first step in QTL analysis: calculate genotype probabilities
sug <- calc.genoprob(sug, step=1)

# 2nd step in QTL analysis: do the LOD score calculations by interval mapping
out.em <- scanone(sug)
plot(out.em)
phenames(sug)

# permutation test
operm <- scanone(sug, n.perm=1000)
operm
plot(operm)

# 5% and 10% significance thresholds
summary(operm) 
# 20% significance threshold
summary(operm, alpha=0.2)

# significant peaks in QTL results
summary(out.em, perms=operm, alpha=0.5, pvalues=TRUE)

# genome scan with multiple traits
phenames(sug)
out.all <- scanone(sug, pheno.col=1:4)
head(out.all)
dim(out.all)

# plot the LOD curves
plot(out.all, lodcolumn=1:3)
plot(out.all, lodcolumn=4, col="green", add=TRUE)
legend("topleft", lwd=2, 
       col=c("black", "blue", "red", "green"),
       phenames(sug)[1:4])

# permutation test for all 4 traits
operm.all <- scanone(sug, pheno.col = 1:4, n.perm=1000, n.cluster=8)
summary(operm.all)

summary(out.all, perms=operm.all, alpha=0.2,
        format="tabByChr", pvalues=TRUE)

# lod support intervals
lodint(out.all, lodcolumn=1, chr=7)
lodint(out.all, lodcolumn=1, chr=7, drop=2)
lodint(out.all, lodcolumn=1, chr=7, drop=1)
lodint(out.all, lodcolumn=1, chr=15)

# Haley-Knott regression
out.all.hk <- scanone(sug, pheno.col=1:4,
                      method="hk")
# permutation test of that
operm.all.hk <- scanone(sug, pheno.col=1:4,
                        method="hk", 
                        n.perm=1000,
                        n.cluster=8)
plot(out.all, out.all.hk, lodcolumn=1, 
     lty=1:2, col=c("slateblue", "violetred"))

# imputation method
sug <- sim.geno(sug, step=1, n.draws=32)
out.all.imp <- scanone(sug, method="imp",
                       pheno.col=1:4)
plot(out.all, out.all.hk, out.all.imp)


# [lunch]

# install the R/qtlcharts package
install.packages("qtlcharts")

# load the R/qtlcharts package
library(qtlcharts)

# interactive plot of genetic map
iplotMap(sug)

# interactive LOD curve plot
iplotScanone(out.all)
iplotScanone(out.all, sug, chr=c(2,7,11,15))
iplotScanone(out.all, sug, chr=c(2,7,11,15), 
             lodcolumn=2, pheno.col=2)
iplotScanone(out.all, sug, chr=c(2,7,11,15),pxgtype="raw")

# show marker names
plotMap(sug, show.marker.names=TRUE)

# non-parametric genome scan for the bp trait
out.np <- scanone(sug, model="np")
plot(out.np)
plot(out.em, col="Orchid", lty=2, add=TRUE)

# create a binary trait as bp > median
bp <- pull.pheno(sug, pheno.col="bp")
bp_bin <- (bp > median(bp, na.rm=TRUE))*1
out.bin <- scanone(sug, pheno.col=bp_bin,
                   model="binary")
plot(out.bin, col="green", lty=3, add=TRUE)

# qtl scan with covariates
bw <- pull.pheno(sug, pheno.col="bw")
out.hw.bwadd <- scanone(sug, pheno.col="heart_wt",
                        addcovar=bw)
plot(out.hw.bwadd)

out.hw.bwint <- scanone(sug, pheno.col="heart_wt",
                        addcovar=bw, intcovar=bw)
plot(out.hw.bwint)

out.hw.bwi <- out.hw.bwint - out.hw.bwadd
plot(out.hw.bwi)

# two-dimensional scan 
