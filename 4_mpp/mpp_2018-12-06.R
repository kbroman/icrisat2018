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
