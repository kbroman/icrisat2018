## Multiple QTL mapping

## 1: load data + libraries

# load R/qtl and R/qtlcharts
library(qtl)
library(qtlcharts)

# load example data set from web
sug <- read.cross("csv", "http://rqtl.org", "sug.csv",
                  genotypes=c("CC", "CB", "BB"), alleles=c("C", "B"))





## 12. Two-dimensional scan
# use a more-coarse step size
sug <- calc.genoprob(sug, step=2.5)

# the 2d scan
out2 <- scantwo(sug, method="hk")

# permutations are really important, but take a very long time
load(url("http://rqtl.org/various.RData"))
summary(operm2)

# summary of scantwo results, with alpha=0.2
summary(out2, perms=operm2, alpha=0.2)

# just chr 7, 12, 15
plot(out2, lower="fv1", upper="int")

# interactive plot
iplotScantwo(out2, sug)

## 13. include marker as covariate
plot(out.hk)
max(out.hk)

# go back to more-fine step size
sug <- calc.genoprob(sug, step=1)

# form covariate matrix for the chr 7 QTL
X <- formMarkerCovar(sug, "7@47.7")

# for an intercross, this is a matrix with 2 columns
head(X)

# scan conditional on the chr 7 QTL
out.c7 <- scanone(sug, method="hk", addcovar=X)
plot(out.hk, out.c7, col=c("slateblue", "orchid"))

## 14. multiple-QTL models
# make a QTL object just with the chr 7 QTL
qtl <- makeqtl(sug, 7, 47.7, what="prob")

# fit the single-QTL model
out.fq <- fitqtl(sug, qtl=qtl, method="hk")
summary(out.fq)

# summary without p-values
summary(out.fq, pvalues=FALSE)

# scan for an additional QTL
out.aq <- addqtl(sug, qtl=qtl, method="hk")
plot(out.aq)

# add to the QTL object
max(out.aq)
qtl <- addtoqtl(sug, qtl=qtl, chr=15, pos=13)

# fit qtl model
out.fq2 <- fitqtl(sug, qtl=qtl, method="hk")
summary(out.fq2, pvalues=FALSE)

# refine the QTL model
rqtl <- refineqtl(sug, qtl=qtl, method="hk")

# plot QTL on the genetic map
plot(rqtl)

# plot LOD profiles
plotLodProfile(rqtl)

# add the single-QTL scan
plot(out.hk, chr=c(7,15), add=TRUE, col="orchid", lty=2)

## 15. stepwise QTL analysis
summary(operm)
out.sq <- stepwiseqtl(sug, method="hk", max.qtl=5,
                      additive.only=TRUE,
                      penalties=c(3.45, Inf, Inf))
out.sq

(pen <- calc.penalties(operm2))
out.sqi <- stepwiseqtl(sug, method="hk", max.qtl=5,
                       penalties=pen)
out.sqi
