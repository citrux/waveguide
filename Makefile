all: scheme plots text

text:
	latexmk -pdf text.tex
	cp text.pdf ~/Dropbox/Public/works/electrodynamics_mw.pdf

dispersion_relation:
	g++ -std=c++11 -lmgl dispersion_relation.cpp -o dispersion_relation

plots: dispersion_relation
	./dispersion_relation

scheme:
	asy scheme.asy

clean:
	rm -f dispersion_relation scheme.pdf
	latexmk -C

