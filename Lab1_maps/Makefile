all: geneticmaps.pdf

geneticmaps.pdf: geneticmaps.dvi Makefile geneticmaps.R
	dvipdf geneticmaps.dvi geneticmaps.pdf

geneticmaps.dvi: geneticmaps.tex Makefile
	latex geneticmaps
	latex geneticmaps
	latex geneticmaps

geneticmaps.tex: geneticmaps.Rnw Makefile clean_sweave.pl
	echo "library(tools); Sweave(\"geneticmaps.Rnw\", pdf=FALSE)" | R --no-save --no-restore --quiet
	clean_sweave.pl geneticmaps

geneticmaps.R: geneticmaps.Rnw Makefile
	echo "library(tools); Stangle(\"geneticmaps.Rnw\")" | R --no-save --no-restore --quiet

clean:
	\rm -r *~ *.bak *.dvi *.log *.aux *.tex
