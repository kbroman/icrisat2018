R_OPTS=--no-save --no-restore --no-init-file --no-site-file # --vanilla, but without --no-environ

pdf: ../docs/lec4a_clearcode.pdf ../docs/lec4a_clearcode_withnotes.pdf

../docs/lec4a_clearcode.pdf: lec4a_clearcode.tex header.tex
	xelatex lec4a_clearcode
	mv lec4a_clearcode.pdf $@

../docs/lec4a_clearcode_withnotes.pdf: lec4a_clearcode_withnotes.tex header.tex
	xelatex lec4a_clearcode_withnotes
	pdfnup lec4a_clearcode_withnotes.pdf --nup 1x2 --no-landscape --paper letterpaper --frame true --scale 0.9
	mv lec4a_clearcode_withnotes-nup.pdf $@

lec4a_clearcode_withnotes.tex: lec4a_clearcode.tex createVersionWithNotes.rb
	createVersionWithNotes.rb lec4a_clearcode.tex lec4a_clearcode_withnotes.tex
