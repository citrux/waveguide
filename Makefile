all: scheme text

text:
	latexmk -pdf text.tex
	cp text.pdf ~/Dropbox/Public/works/electrodynamics_mw.pdf

scheme:
	asy scheme.asy

clean:
	rm -f scheme.pdf
	latexmk -C

