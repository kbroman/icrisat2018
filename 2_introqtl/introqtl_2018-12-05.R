# Intro to QTL mapping with R/qtl
# Live coding at ICRISAT workshop 2018-12-05
#
# largely following http://rqtl.org/tutorials/rqtltour2.pdf
# latest version of this script at https://bit.ly/introqtl2018

# load R/qtl
library(qtl)

# import some QTL mapping data
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

# 5% and 10% significance thresholds
summary(operm) 
# 20% significance threshold
summary(operm, alpha=0.2)

# significant peaks in QTL results
summary(out.em, perms=operm, alpha=0.2)
