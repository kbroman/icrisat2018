# illustration of a DO genome

blue_background <- TRUE # make TRUE if you want a blue background

# libraries
library(simcross) # install with devtools::install_github("kbroman/simcross")
library(qtl)      # install with install.packages("qtl")
library(broman)   # install with install.packages("broman")

blue <- rgb(0, 0, 80, maxColorValue=255)

n_gen <- 20

cache_dir <- "_cache"
if(!dir.exists(cache_dir)) dir.create(cache_dir)
file <- file.path(cache_dir, "dosim.RData")
if(file.exists(file)) {
    load(file)
} else {
    set.seed(71465056)

    # generate a DO pedigree
    doped <- sim_do_pedigree(ngen=n_gen, nkids_per=1)

    # chr lengths -> map (from the G2F1 map, Liu et al 2012
    chr_L <- c(98.46, 103.741, 81.891, 88.605, 89.235, 78.76, 88.98,
               76.095, 75.052, 77.349, 87.885, 63.75, 67.21, 66.403, 58.876,
               57.427, 60.689, 59.054, 56.886, 78.781)

    map <- sim.map(chr_L, n.mar=2, anchor=TRUE, eq=TRUE)

    # simulate genotypes
    dosim <- sim_from_pedigree_allchr(doped, map)
    # subset to the last generation
    dosim <- lapply(dosim, function(a) a[doped$do & doped$gen==n_gen])

    # sex
    dosex <- setNames(doped$sex[doped$do & doped$gen==n_gen],
                      doped$id[doped$do & doped$gen==n_gen])

    # collapse alleles
    dosim <- lapply(dosim, collapse_do_alleles)

    # save to file
    save(dosim, doped, dosex, file=file)
}

chr_range <- vapply(dosim, function(a) range(a[[1]]$mat$locations), c(1,1))
chr_L <- apply(chr_range, 2, diff)

ind <- 2

pdf("../Figs/do_genome.pdf", width=9.75, height=6.5, pointsize=16)
par(mar=rep(0.1,4), bty="n")
if(blue_background) {
    par(fg="white",col="white",col.axis="white",col.lab="white",bg=blue)
} else {
    par(fg=blue,col=blue,col.axis=blue,col.lab=blue,bg="white")
}
top_gap <- 7.5
bottom_gap <- 2.5
plot(0,0,type="n", xlim=c(0.5, 20.5), ylim=c(-bottom_gap, 100+top_gap),
     xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i")
for(i in 1:20) {
    this_L <- 100*chr_L[i]/max(chr_L)

    # kludge to deal with X chr in males
    if(i==20 && dosex[ind]==1) {
        chrlength <- c(this_L, this_L/10)
    } else {
        chrlength <- this_L
    }

    plot_ind(dosim[[i]][[ind]], c(i, 100-this_L/2), chrwidth=0.2, gap=0.1,
             chrlength=chrlength, col=simcross::CCcolors())
}
text(1:20, rep(100 + top_gap/2, 20), c(1:19,"X"), adj=0.5)
dev.off()
