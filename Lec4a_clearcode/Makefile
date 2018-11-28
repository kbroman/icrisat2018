R_OPTS=--no-save --no-restore --no-init-file --no-site-file # --vanilla, but without --no-environ

pdf: ../docs/clearcode.pdf ../docs/clearcode_withnotes.pdf

../docs/clearcode.pdf: clearcode.tex header.tex
	xelatex clearcode
	mv clearcode.pdf $@

../docs/clearcode_withnotes.pdf: clearcode_withnotes.tex header.tex
	xelatex clearcode_withnotes
	pdfnup clearcode_withnotes.pdf --nup 1x2 --no-landscape --paper letterpaper --frame true --scale 0.9
	mv clearcode_withnotes-nup.pdf $@

clearcode_withnotes.tex: clearcode.tex createVersionWithNotes.rb
	createVersionWithNotes.rb clearcode.tex clearcode_withnotes.tex
