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




