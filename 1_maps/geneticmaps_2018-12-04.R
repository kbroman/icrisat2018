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
g <- pull.geno(mapthis)
dim(g)
gfreq <- matrix(ncol=3, nrow=nrow(g))
for(i in 1:3) {
   gfreq[,i] <- rowMeans( g == i, na.rm =TRUE)
}
head(gfreq)
plot(gfreq[,1])
plot(gfreq[,2])
plot(gfreq[,3])

pA <- gfreq[,1] + gfreq[,2]/2
plot(pA)
hist(pA, breaks=50)

## pairwise linkages between markers
mapthis <- est.rf(mapthis)
checkAlleles(mapthis)

# pull out all recombination fractions and LOD scores
rf <- pull.rf(mapthis)
lod <- pull.rf(mapthis, "lod")
plot(as.numeric(rf), as.numeric(lod))
rf[1:5, 1:5]

# there are messed up markers
# but we'll proceed and fix them later
lg <- formLinkageGroups(mapthis, max.rf=0.35, min.lod=10)
head(lg)
table(lg[,2])

# play around particularly with min.lod
lg_lod3 <- formLinkageGroups(mapthis, max.rf=0.35, min.lod=3)
table(lg_lod3[,2])

lg_lod20 <- formLinkageGroups(mapthis, max.rf=0.35, min.lod=20)
table(lg_lod20[,2])

lg_lod50 <- formLinkageGroups(mapthis, max.rf=0.35, min.lod=50)
table(lg_lod50[,2])

# reorganize the markers into linkage groups
mapthis2 <- formLinkageGroups(mapthis, max.rf=0.35, min.lod=10, 
                              reorgMarkers=TRUE)
plotRF(mapthis2)

# looks like linkage groups 5, 7, 8, 9, 10, 11, 12 need alleles swapped

# look at a marker on lg 4 and on lg 5
mn4 <- markernames(mapthis2, chr=4)
mn5 <- markernames(mapthis2, chr=5)
mn4
mn5
# marker x marker genotype table
geno.crosstab(mapthis2, mn4[1], mn5[1])

# swap the alleles at markers on lg 5, 7-12
toswitch <- markernames(mapthis2, chr=c(5, 7:12))
mapthis2 <- switchAlleles(mapthis2, toswitch)
