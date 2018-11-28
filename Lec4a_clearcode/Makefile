R_OPTS=--no-save --no-restore --no-init-file --no-site-file # --vanilla, but without --no-environ

pdf: clearcode.pdf clearcode_withnotes.pdf

clearcode.pdf: clearcode.tex header.tex
	xelatex clearcode

clearcode_withnotes.pdf: clearcode_withnotes.tex header.tex
	xelatex clearcode_withnotes
	pdfnup clearcode_withnotes.pdf --nup 1x2 --no-landscape --paper letterpaper --frame true --scale 0.9
	mv clearcode_withnotes-nup.pdf clearcode_withnotes.pdf

clearcode_withnotes.tex: clearcode.tex createVersionWithNotes.rb
	createVersionWithNotes.rb clearcode.tex clearcode_withnotes.tex
