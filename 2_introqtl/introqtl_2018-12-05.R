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
