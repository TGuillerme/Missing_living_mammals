all: Missing_mammals-draft.pdf

Missing_mammals-draft.pdf: Missing_mammals-draft.tex
	pdflatex Missing_mammals-draft
	bibtex Missing_mammals-draft
	pdflatex Missing_mammals-draft
	pdflatex Missing_mammals-draft

tidy:
	rm -f Missing_mammals-draft.{aux,bbl,blg,log}

clean: tidy
	rm -f Missing_mammals-draft.pdf

.PHONY: all tidy clean
