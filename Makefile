R_OPTS=--no-save --no-restore --no-init-file --no-site-file # --vanilla, but without --no-environ

all: docs/geneticmaps.pdf \
	 docs/2_introqtl.pdf \
	 docs/3_multiqtl.pdf \
	 docs/4_mpp.pdf \
	 docs/5a_clearcode.pdf \
	 docs/5a_clearcode_withnotes.pdf \
	 docs/5b_rpack.pdf \
	 docs/5b_rpack_withnotes.pdf

docs/geneticmaps.pdf: 1_maps/geneticmaps.dvi 1_maps/geneticmaps.R
	cd 1_maps;dvipdf geneticmaps.dvi geneticmaps.pdf
	mv 1_maps/geneticmaps.pdf $@

1_maps/geneticmaps.dvi: 1_maps/geneticmaps.tex
	cd 1_maps;latex geneticmaps;latex geneticmaps;latex geneticmaps

1_maps/geneticmaps.tex: 1_maps/geneticmaps.Rnw 1_maps/clean_sweave.pl
	echo "library(tools); Sweave(\"geneticmaps.Rnw\", pdf=FALSE)" | cd 1_maps; R --no-save --no-restore --quiet
	cd 1_maps; clean_sweave.pl geneticmaps

1-maps/geneticmaps.R: 1_maps/geneticmaps.Rnw
	echo "library(tools); Stangle(\"geneticmaps.Rnw\")" | cd 1_maps;R --no-save --no-restore --quiet

docs/2_introqtl.pdf: 2_introqtl/2_introqtl.tex
	cd 2_introqtl;pdflatex $(<F)
	mv 2_introqtl/2_introqtl.pdf $@

docs/3_multiqtl.pdf: 3_multiqtl/3_multiqtl.tex
	cd 3_multiqtl;pdflatex $(<F)
	mv 3_multiqtl/3_multiqtl.pdf $@

docs/4_mpp.pdf: 4_mpp/4_mpp.tex
	cd 4_mpp;pdflatex $(<F)
	mv 4_mpp/4_mpp.pdf $@

docs/5a_clearcode.pdf: 5a_clearcode/5a_clearcode.tex latex/header.tex
	cd 5a_clearcode;xelatex 5a_clearcode
	mv 5a_clearcode/5a_clearcode.pdf $@

docs/5a_clearcode_withnotes.pdf: 5a_clearcode/5a_clearcode_withnotes.tex latex/header.tex
	cd 5a_clearcode;xelatex 5a_clearcode_withnotes
	cd 5a_clearcode;pdfnup 5a_clearcode_withnotes.pdf --nup 1x2 --no-landscape --paper letterpaper --frame true --scale 0.9
	mv 5a_clearcode/5a_clearcode_withnotes-nup.pdf $@

5a_clearcode/5a_clearcode_withnotes.tex: 5a_clearcode/5a_clearcode.tex
	ruby/createVersionWithNotes.rb $< $@

docs/5b_rpack.pdf: 5b_packages/5b_rpack.tex latex/header.tex
	cd 5b_packages;xelatex 5b_rpack
	mv 5b_packages/5b_rpack.pdf $@

docs/5b_rpack_withnotes.pdf: 5b_packages/5b_rpack_withnotes.tex latex/header.tex
	cd 5b_packages;xelatex 5b_rpack_withnotes
	cd 5b_packages;pdfnup 5b_rpack_withnotes.pdf --nup 1x2 --no-landscape --paper letterpaper --frame true --scale 0.9
	mv 5b_packages/5b_rpack_withnotes-nup.pdf $@

5b_packages/5b_rpack_withnotes.tex: 5b_packages/5b_rpack.tex
	ruby/createVersionWithNotes.rb $< $@
