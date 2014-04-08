%.pdf: %.tex
	pdflatex $^

all: $(patsubst %.tex, %.pdf, $(wildcard *.tex))
