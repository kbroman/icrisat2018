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
nt_ind <- ntyped(mapthis)
plot(nt_ind)

# number of genotypes per marker
nt_mar <- ntyped(mapthis, "mar")
plot(nt_mar)

# omit individuals that have > 50% missing data
mapthis <- subset(mapthis, ind=(nt_ind > 50))
nind(mapthis)

# omit markers that have < 200 genotypes
bad_markers <- names(nt_mar)[nt_mar < 200]
bad_markers
mapthis <- drop.markers(mapthis, bad_markers)
nmar(mapthis)

# look for sample duplicates
cg <- comparegeno(mapthis)
cg[1:5,1:5]
hist(cg, breaks=100)
