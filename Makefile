LATEX = pdflatex
BIBTEX = bibtex

SRC = $(wildcard src-*.tex)

.PHONY: all bib clean

all: paper.pdf

paper.pdf: $(SRC) macros.tex paper.tex
	$(LATEX) paper.tex

bib:
	$(LATEX) paper.tex
	$(BIBTEX) paper.aux
	$(LATEX) paper.tex
	$(LATEX) paper.tex

clean:
	@find . -name "*.aux" -type f -delete
	@find . -name "*.bbl" -type f -delete
	@find . -name "*.blg" -type f -delete
	@find . -name "*.log" -type f -delete
	@find . -name "*.out" -type f -delete
	@find . -name "*.toc" -type f -delete
	@find . -name "*.pdf" -type f -delete
	@echo "Deleted .aux .bbl .blg .log .out .toc .pdf files."
