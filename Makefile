%.pdf: %.tex
	pdflatex $^

rings.pdf: rings.md
	pandoc --template=./compact.latex -o $@ $<

all: $(patsubst %.tex, %.pdf, $(wildcard *.tex)) rings


