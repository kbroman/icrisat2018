# this will load qtl2geno, qtl2scan, qtl2plot
library(qtl2)

# if you don't have the sug data anymore, reload it
library(qtl)
sug <- read.cross("csv", "http://rqtl.org", "sug.csv",
                  genotypes=c("CC", "CB", "BB"),
                  alleles=c("C", "B"))

# convert data to R/qtl2 format
sug2 <- convert2cross2(sug)

# summary
summary(sug2)

# other summary functions
n_ind(sug2)
n_pheno(sug2)
n_mar(sug2)
tot_mar(sug2)

# prep for calculating genotype probs
gmap <- insert_pseudomarkers(sug2$gmap, step=1)

# calculate genotype probabilities
pr <- calc_genoprob(sug2, gmap)

# genome scan by Haley-Knott regression
out <- scan1(pr, sug2$pheno[,1:4])

# plot the LOD curves for the first phenotype
plot(out, gmap)

# add lod curves for 2nd phenotype
plot(out, gmap, lod=2, col="orchid", add=TRUE)

# find peaks above 4
find_peaks(out, gmap, threshold=4)

# also include 1.5-LOD support intervals
find_peaks(out, gmap, threshold=4, drop=1.5)

# permutations
operm <- scan1perm(pr, sug2$pheno[,1:4], n_perm=400, cores=0)

# 5% significance threshold
summary(operm)

# find peaks again, with threshold=3.5
find_peaks(out, gmap, threshold=3.5, drop=1.5)

## linear mixed models
# calculate kinship matrix
k <- calc_kinship(pr)

# calculate kinship by "loco" method
k_loco <- calc_kinship(pr, "loco")

# genome scan by LMM
out_lmm <- scan1(pr, sug2$pheno[,1:4], k)

# genome scan by LMM with "loco" method
out_loco <- scan1(pr, sug2$pheno[,1:4], k_loco)

# plot them together
(ymax <- max(c(out, out_lmm, out_loco)))
plot(out, gmap, ylim=c(0, 7))
plot(out_lmm, gmap, col="orchid", lty=2, add=TRUE)
plot(out_loco, gmap, col="green3", lty=2, add=TRUE)

## estimating QTL effects
# chr 7, first phenotype
eff7 <- scan1coef(pr[,7], sug2$pheno[,1])
plot(eff7, gmap[7], columns=1:3, col=c("slateblue", "orchid", "green3"))

# also include the LOD curve
plot(eff7, gmap[7], columns=1:3, col=c("slateblue", "orchid", "green3"), scan1_output=out)

# chr 15, first phenotype
eff15 <- scan1coef(pr[,15], sug2$pheno[,1])
plot(eff15, gmap[15], columns=1:3, col=c("slateblue", "orchid", "green3"), scan1_output=out)

# estimating effects with LOCO method
eff7_loco <- scan1coef(pr[,7], sug2$pheno[,1], k_loco[7])
plot(eff7, gmap[7], columns=1:3, col=c("slateblue", "orchid", "green3"))
plot(eff7_loco, gmap[7], columns=1:3,
     col=c("slateblue", "orchid", "green3"), lty=2, add=TRUE)

## how about the additive and dominance effects?
# need to provide a matrix of contrasts
contr <- cbind(mu=c(1,1,1), a=c(-0.5, 0, 0.5), d=c(-0.5, 1, -0.5))
eff7ad <- scan1coef(pr[,7], sug2$pheno[,1], contrasts=contr)
plot(eff7ad, gmap[7], columns=2:3, scan1_output=out)

# chr 15
eff15ad <- scan1coef(pr[,15], sug2$pheno[,1], contrasts=contr)
plot(eff15ad, gmap[15], columns=2:3, scan1_output=out)

## phenotype:genotype relationship at the inferred QTL
# peak location on chr 7
max(out, gmap, chr=7)
# inferred genotypes
g7 <- maxmarg(pr, gmap, chr=7, pos=47.71, return_char=TRUE, minprob=0.75)
# plot genotype x phenotype
plot_pxg(g7, sug2$pheno[,1], ylab="blood pressure")
# don't sort the genotypes
plot_pxg(g7, sug2$pheno[,1], ylab="blood pressure", sort=FALSE)
# add confidence intervals
plot_pxg(g7, sug2$pheno[,1], ylab="blood pressure", sort=FALSE, SEmult=2)
# make it horizontal
plot_pxg(g7, sug2$pheno[,1], xlab="blood pressure", sort=FALSE, SEmult=2, swap_axes=TRUE)
# just show the confidence intervals
plot_pxg(g7, sug2$pheno[,1], xlab="blood pressure", sort=FALSE, SEmult=2, omit_points=TRUE)


# peak location on chr 15
max(out, gmap, chr=15)
# inferred genotypes
g15 <- maxmarg(pr, gmap, chr=15, pos=11.96, return_char=TRUE, minprob=0.75)
# plot genotype x phenotype
plot_pxg(g15, sug2$pheno[,1], ylab="blood pressure", sort=FALSE, SEmult=2)
# just the confidence intervals
plot_pxg(g15, sug2$pheno[,1], xlab="blood pressure", sort=FALSE, SEmult=2, omit_points=TRUE)


######################################################################

## Diversity outbreds
DOex <- read_cross2("http://rqtl.org/DOex.zip")
summary(DOex)

# insert a few pseudomarkers, so that gaps are < 1 cM
do_gmap <- insert_pseudomarkers(DOex$gmap, step=1, stepwidth="max")

# genetic -> physical map
do_pmap <- interp_map(do_gmap, DOex$gmap, DOex$pmap)

# calculate genotype probabilities
pr <- calc_genoprob(DOex, do_gmap, cores=0)

# plot some reconstructed genomes
plot_genoprob(pr, do_pmap, ind=1, chr="2")
# show only genotypes that achieve prob > 0.25
plot_genoprob(pr, do_pmap, ind=1, chr="2", threshold=0.25)
plot_genoprob(pr, do_pmap, ind=2, chr="2", threshold=0.25)

# infer genotypes
m <- maxmarg(pr, minprob=0.5)
v <- guess_phase(DOex, m)

# plot reconstructed genome
plot_onegeno(v, do_pmap, ind=1)
plot_onegeno(v, do_pmap, ind=2)
plot_onegeno(v, do_pmap, ind=53)

# reduce to allele dosages
apr <- genoprob_to_alleleprob(pr)

# calculate kinship matrices using "loco" method
# but note we just have chromosomes 2, 3, X
k <- calc_kinship(apr, "loco")

# let's use sex as an additive covariate
# make sure the IDs are in the names
sex <- (DOex$covar$Sex=="male")*1
names(sex) <- rownames(DOex$covar)

# QTL scan
out <- scan1(apr, DOex$pheno[,1], k, sex)
plot(out, do_gmap)

# there's a strong peak on chromosome 2, let's look at effects
eff2 <- scan1coef(apr[,"2"], DOex$pheno[,1], k["2"], sex)
plot_coefCC(eff2, do_gmap["2"])

# add a legend
plot_coefCC(eff2, do_gmap["2"], legend="bottomleft")

# also include LOD curves
plot_coefCC(eff2, do_gmap["2"], legend="bottomleft", scan1_output=out)

# we can also use "BLUP"
blup2 <- scan1blup(apr[,"2"], DOex$pheno[,1], k["2"], sex)
plot_coefCC(blup2, do_gmap["2"], legend="bottomleft", scan1_output=out)

# phenotype x genotype at a particular position

# phenotype average x genotype at a particular position

## snp associations in region
# first find the position
marker <- rownames(max(out, DOex$gmap, chr="2"))
peak_Mbp <- DOex$pmap[["2"]][marker]

# connect to SNP database
qsnps <- create_variant_query_func(system.file("extdata", "cc_variants_small.sqlite", package="qtl2"))

# scan chr 2 region
out_snps <- scan1snps(pr, do_pmap, DOex$pheno[,1], query_func=qsnps, keep_all_snps=TRUE, chr=2)


# plot all SNPs; highlight the ones within 1.5 LOD of top
plot(out_snps$lod, out_snps$snpinfo, drop_hilit=1.5)

# top SNP
top_snps(out_snps$lod, out_snps$snpinfo, drop=0, show_all_snps=FALSE)

# top SNP + equivalent ones
top_snps(out_snps$lod, out_snps$snpinfo, drop=0)

# indexed SNPs within 1.5 of top
top_snps(out_snps$lod, out_snps$snpinfo, show_all_snps=FALSE)

# also mouse genes in here


#####################################################################

## MAGIC lines

file <- paste0("https://raw.githubusercontent.com/rqtl/",
               "qtl2data/master/ArabMAGIC/arabmagic.zip")
arab <- read_cross2(file)

gmap <- insert_pseudomarkers(arab$gmap, step=1)
pmap <- interp_map(gmap, arab$gmap, arab$pmap)

pr <- calc_genoprob(arab, gmap, error_prob=0.002, cores=0)

# genotype probabilities along a chromosome

# reconstructed genome
color <- sample(colors(), 19)

# genome scan

# estimated effects at a particular locus

# phenotype x genotype


# assembly the snp info from the founder genotypes
install.packages("qtl2convert", repos="http://rqtl.org/qtl2cran")
library(qtl2convert)

snpinfo <- map_list_to_df(arab$pmap, marker="snp")
fg <- do.call("cbind", arab$founder_geno)
sdp <- calc_sdp(t(fg))

snpinfo <- cbind(snpinfo[names(sdp),], sdp)
head(snpinfo)

snpinfo <- index_snps(pmap, snpinfo)

out_snps <- scan1snps(pr, pmap, arab$pheno, snpinfo=snpinfo, keep_all_snps=TRUE)
plot(out_snps$lod, out_snps$snpinfo, gap=5, lod=8)
