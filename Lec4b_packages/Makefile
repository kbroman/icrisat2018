R_OPTS=--no-save --no-restore --no-init-file --no-site-file # --vanilla, but without --no-environ

pdf: ../docs/rpack.pdf ../docs/rpack_withnotes.pdf

../docs/rpack.pdf: rpack.tex header.tex
	xelatex rpack
	mv rpack.pdf $@

../docs/rpack_withnotes.pdf: rpack_withnotes.tex header.tex
	xelatex rpack_withnotes
	pdfnup rpack_withnotes.pdf --nup 1x2 --no-landscape --paper letterpaper --frame true --scale 0.9
	mv rpack_withnotes-nup.pdf $@

rpack_withnotes.tex: rpack.tex createVersionWithNotes.rb
	createVersionWithNotes.rb rpack.tex rpack_withnotes.tex
