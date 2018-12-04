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
rug(cg)
# get the row and column of values > 0.9
which(cg > 0.9 & !is.na(cg), arr.ind=TRUE)

# how similar are they?
g <- pull.geno(mapthis)
table(g[144,], g[292,]) # note 1=AA, 2=AB, 3=BB
cg[144,292] # share 96.7% of genotypes
table(g[214,], g[216,]) # note 1=AA, 2=AB, 3=BB
cg[214,216] # share % of genotypes
table(g[238,], g[288,]) # note 1=AA, 2=AB, 3=BB
cg[238,288] # share % of genotypes

# omit one of each pair of duplicates
mapthis <- subset(mapthis, ind= -c(292, 216, 288))
nind(mapthis)

# dendrogram for the individuals
# distance 
d <- 1 - cg
diag(d) <- 0
d[1:5,1:5]
hclust_out <- hclust( as.dist(d) )
plot(hclust_out)

# look for segregation distortion
gt <- geno.table(mapthis)
head(gt)
gt[gt$P.value < 0.05/nrow(gt), ]
# throw out markers with P.values < 10^-10
bad_markers <- rownames(gt)[gt$P.value < 1e-10]
mapthis <- drop.markers(mapthis, bad_markers)

# genotype frequencies in the individuals
