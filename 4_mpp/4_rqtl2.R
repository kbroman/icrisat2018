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
map <- insert_pseudomarkers(sug2$gmap, step=1)

# calculate genotype probabilities
pr <- calc_genoprob(sug2, map)

# genome scan by Haley-Knott regression
out <- scan1(pr, sug2$pheno[,1:4])

# plot the LOD curves for the first phenotype
plot(out, map)

# add lod curves for 2nd phenotype
plot(out, map, lod=2, col="orchid", add=TRUE)

# find peaks above 4
find_peaks(out, map, threshold=4)

# permutations
operm <- scan1perm(pr, sug2$pheno[,1:4], n_perm=200, cores=0)

# 5% significance threshold
summary(operm)

# find peaks again, with threshold=3.5
find_peaks(out, map, threshold=3.5)

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
plot(out, map, ylim=c(0, 7))
plot(out_lmm, map, col="orchid", lty=2, add=TRUE)
plot(out_loco, map, col="orange", lty=2, add=TRUE)

## X chromosome: may need special covariates, 
## and in R/qtl2 need to provide them directly
# but here: no X chromosome, and no need for special covariates
# (but no harm in doing so)
(Xcovar <- get_x_covar(sug2))
out <- scan1(pr, sug2$pheno, Xcovar=Xcovar)

## estimating QTL effects
# chr 7, first phenotype
eff7 <- scan1coef(pr[,7], sug2$pheno[,1])
plot(eff7, map[7], columns=1:3, col=c("slateblue", "orchid", "green3"))

# chr 15, first phenotype
eff15 <- scan1coef(pr[,15], sug2$pheno[,1])
plot(eff15, map[15], columns=1:3, col=c("slateblue", "orchid", "green3"))

# estimating effects with LOCO method
eff7_loco <- scan1coef(pr[,7], sug2$pheno[,1], k_loco[7])
plot(eff7, map[7], columns=1:3, col=c("slateblue", "orchid", "green3"))
plot(eff7_loco, map[7], columns=1:3, 
     col=c("slateblue", "orchid", "green3"), lty=2, add=TRUE)

## how about the additive and dominance effects?
# need to provide a matrix of contrasts
contr <- cbind(mu=c(1,1,1), a=c(-0.5, 0, 0.5), d=c(-0.5, 1, -0.5))
eff7ad <- scan1coef(pr[,7], sug2$pheno[,1], contrasts=contr)
plot(eff7ad, map[7], columns=2:3)

# chr 15
eff15ad <- scan1coef(pr[,15], sug2$pheno[,1], contrasts=contr)
plot(eff15ad, map[15], columns=2:3)

## Diversity outbreds
DOex <- read_cross2("http://rqtl.org/DOex.zip")
summary(DOex)

# insert a few pseudomarkers, so that gaps are < 1 cM
do_map <- insert_pseudomarkers(DOex$gmap, step=1, stepwidth="max")

# calculate genotype probabilities
pr <- calc_genoprob(DOex, do_map, cores=0)

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
plot(out, do_map)

# there's a strong peak on chromosome 2, let's look at effects
eff2 <- scan1coef(apr[,"2"], DOex$pheno[,1], k["2"], sex)
plot_coefCC(eff2, do_map["2"])

# add a legend
legend("bottomleft", names(CCcolors), col=CCcolors, lwd=2, ncol=2,
       bg="gray90")

# we can also use "BLUP"
blup2 <- scan1blup(apr[,"2"], DOex$pheno[,1], k["2"], sex)
plot_coefCC(blup2, do_map["2"])

# add a legend
legend("bottomleft", names(CCcolors), col=CCcolors, lwd=2, ncol=2,
       bg="gray90")

## snp associations in region
# first find the position
marker <- rownames(max(out, DOex$gmap, chr="2"))
peak_Mbp <- DOex$pmap[["2"]][marker]

# download SNPs for the region
url <- "http://rqtl.org/c2_snpinfo.rds"
file <- basename(url)
download.file(url, file)
snpinfo <- readRDS(file)

# each row is a SNP; columns 5-12 have founder genotypes as 1/3
head(snpinfo)

# first calculate "strain distribution pattern" encoding
snpinfo$sdp <- calc_sdp(snpinfo[,-(1:4)])
head(snpinfo)

# create index for the SNPs
snpinfo <- index_snps(DOex$pmap, snpinfo)
head(snpinfo)

# convert allele dosages to SNP allele dosages
snp_pr <- genoprob_to_snpprob(apr, snpinfo)

# calculate LOD scores at each distinct SNP
out_snp <- scan1(snp_pr, DOex$pheno[,1], k[["2"]], sex)

# plot the results
plot(out_snp, snpinfo)

# plot just indexed SNPs
plot(out_snp, snpinfo, show_all_snps=FALSE)

# plot all SNPs; highlight the ones within 1.5 LOD of top
plot(out_snp, snpinfo, drop.hilit=1.5)

# top SNP
top_snps(out_snp, snpinfo, drop=0, show_all_snps=FALSE)

# top SNP + equivalent ones
top_snps(out_snp, snpinfo, drop=0)

# indexed SNPs within 1.5 of top
top_snps(out_snp, snpinfo, show_all_snps=FALSE)
