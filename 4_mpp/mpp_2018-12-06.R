# Analysis of Multiparent Populations with R/qtl2
# Live coding at ICRISAT workshop 2018-12-06
#
# largely following https://kbroman.org/qtl2/assets/vignettes/user_guide.html
# latest version of this script at https://bit.ly/rqtl2_2018

# load R/qtl2
library(qtl2)

# load data set from yesterday
library(qtl)
sug <- read.cross("csv", "http://rqtl.org", "sug.csv",
                  genotypes=c("CC", "CB", "BB"),
                  alleles=c("C", "B"))

# convert the data set to R/qtl2 format
sug2 <- convert2cross2(sug)

# summary functions
summary(sug2)
n_ind(sug2)
n_mar(sug2)
tot_mar(sug2)
n_pheno(sug2)

# calculate genotype probabilities
gmap <- insert_pseudomarkers(sug2$gmap, step=1)
pr <- calc_genoprob(sug2, gmap)

# do the genome scan (Haley-Knott regression)
out <- scan1(pr, sug2$pheno[,1:4])  # genome scan on phe 1-4
head(out) # output is just a matrix

# plot the LOD curves for the first phenotype
plot(out, gmap)
plot(out, gmap, lodcolumn=2, col="orchid", add=TRUE)

# find the QTL peaks
find_peaks(out, gmap, threshold=4)
# also include 1.5-LOD support intervals
find_peaks(out, gmap, threshold=4, drop=1.5)

# permutation test (cores=0 means use all available CPUs)
operm <- scan1perm(pr, sug2$pheno, n_perm=400, cores=0)
summary(operm)

# find_peaks with permutation thresholds
thr <- summary(operm)
find_peaks(out, gmap, threshold=thr[1:4], drop=1.5)

# scan with linear model mixed 
## first calculate kinship matrix
k <- calc_kinship(pr)
## kinship matrices by "leave one chromosome out" method (LOCO)
k_loco <- calc_kinship(pr, "loco")

# use these in scan1() function (cores = 0 means use all CPUs)
out_lmm <- scan1(pr, sug2$pheno[,1:4], k, cores=0)
out_loco <- scan1(pr, sug2$pheno[,1:4], k_loco, cores=0)

# plot the results for first phenotype
plot(out, gmap)
plot(out_lmm, gmap, col="orchid", add=TRUE)
plot(out_loco, gmap, col="green3", add=TRUE, lty=2)

# plot the other three traits, in separate panels
par(mfrow=c(3,1), mar=c(3.1,1.1,2.1,0.1))
for(i in 2:4) {
    plot(out, gmap, lodcolumn=i, main=columnes(out)[i])
    plot(out_lmm, gmap, col="orchid", add=TRUE, lodcolumn=i)
    plot(out_loco, gmap, col="green3", add=TRUE, lty=2, lodcolumn=i)
}
