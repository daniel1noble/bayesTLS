# Project Makefile.
# `make` (default target = `all`) renders both ms.qmd and supplement.qmd
# in every format listed in their YAML headers (currently html + docx).
#
# Targets:
#   all        — render ms.qmd and supplement.qmd in all formats
#   ms         — render only ms.qmd
#   supplement — render only supplement.qmd
#   clean      — remove rendered HTML/DOCX outputs (caches preserved)

.PHONY: all ms supplement clean

all: ms supplement

ms:
	quarto render ms/ms.qmd

supplement:
	quarto render ms/supplement.qmd

clean:
	rm -f ms/ms.html ms/ms.docx ms/supplement.html ms/supplement.docx
