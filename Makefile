all: scheme text

text:
	latexmk -pdf -cd report/text.tex
	cp report/text.pdf ~/Dropbox/Public/works/electrodynamics_mw.pdf

scheme:
	cd report && asy scheme.asy

clean:
	rm -f report/scheme.pdf
	cd report && latexmk -C

