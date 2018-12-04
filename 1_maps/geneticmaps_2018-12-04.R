# Genetic Map Construction
# Live coding at ICRISAT workshop 2018-12-04
#
# largely following http://rqtl.org/tutorials/geneticmaps.pdf
# latest version of this script at https://bit.ly/genmaps2018

# load R/qtl library
library(qtl)

# load example data
data(mapthis)

# summary
summary(mapthis)
nind(mapthis)
nmar(mapthis)

## data diagnostics
# plot of missing data pattern
plotMissing(mapthis)

# plot of the genotypes
geno.image(mapthis)

# number of genotypes per individual
