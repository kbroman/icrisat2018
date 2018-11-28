R_OPTS=--no-save --no-restore --no-init-file --no-site-file # --vanilla, but without --no-environ

pdf: ../docs/lec4b_rpack.pdf ../docs/lec4b_rpack_withnotes.pdf

../docs/lec4b_rpack.pdf: lec4b_rpack.tex header.tex
	xelatex lec4b_rpack
	mv lec4b_rpack.pdf $@

../docs/lec4b_rpack_withnotes.pdf: lec4b_rpack_withnotes.tex header.tex
	xelatex lec4b_rpack_withnotes
	pdfnup lec4b_rpack_withnotes.pdf --nup 1x2 --no-landscape --paper letterpaper --frame true --scale 0.9
	mv lec4b_rpack_withnotes-nup.pdf $@

lec4b_rpack_withnotes.tex: lec4b_rpack.tex createVersionWithNotes.rb
	createVersionWithNotes.rb lec4b_rpack.tex lec4b_rpack_withnotes.tex
