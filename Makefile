text:
	latexmk -pdf text.tex
	cp text.pdf ~/Dropbox/Public/works/electrodynamics_mw.pdf

dispersion_relation:
	g++ -std=c++11 -lmgl dispersion_relation.cpp -o dispersion_relation

clean:
	rm $(NAME)
	latexmk -C

