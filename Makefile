# Project Makefile.
#
# The project is now a Quarto book (see _quarto.yml) that renders ms.qmd and
# supplement.qmd together with shared cross-referencing. Main objects use flat
# built-in numbering (Figure 1, Table 1); supplement figures/tables use custom
# Quarto floats (`suppfig-*`, `supptbl-*`) for Figure S1 / Table S1 labels.
#
# Targets:
#   all        — alias for `book` (renders the book in all formats listed in _quarto.yml)
#   book       — render the book (HTML + DOCX) into _book/
#   ms         — render ONLY ms.qmd (book context, single chapter)
#   supplement — render ONLY supplement.qmd (book context, single chapter)
#   clean      — remove rendered outputs

.PHONY: all book ms supplement clean

all: book

book:
	quarto render

ms:
	quarto render ms/ms.qmd

supplement:
	quarto render ms/supplement.qmd

clean:
	rm -rf _book ms/ms.html ms/ms.docx ms/supplement.html ms/supplement.docx
