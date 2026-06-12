# Project Makefile.
#
# Two documents, each rendered to HTML, DOCX, and PDF:
#
#   ms    ms/ms.qmd          → _output/ms.{html,docx,pdf}    (manuscript)
#   supp  ms/supplement.qmd  → _output/supp.{html,docx,pdf}  (supplement)
#
# The two documents reference each other via plain text (e.g. "see
# Equation 7 of the manuscript"); there is no merged combined document.
#
# There is intentionally NO project-level `_quarto.yml`. Each .qmd carries
# its own complete YAML, so renders are fully independent.
#
# ---------------------------------------------------------------------------
# Dropbox: this project lives inside a Dropbox-synced directory. Quarto's
# render produces several intermediate files (`*.tex`, `*.log`, `*_files/`,
# `*_cache/`, the .pdf itself) alongside the source. When Dropbox syncs
# those mid-render, it produces "conflicted copy" duplicates.
#
# To prevent that entirely, the render is done OUTSIDE Dropbox in
# $(BUILDROOT) (default: ~/Library/Caches/tls-render). Sources are rsync'd
# in, `bib/` and `output/` are symlinked, and only the final output is
# copied back into the Dropbox-synced project's `_output/`. Build caches
# persist in $(BUILDROOT) across renders — `make clean` doesn't touch
# them. Run `make build-clean` to wipe build caches.
# ---------------------------------------------------------------------------

OUTDIR    := _output
BUILDROOT := $(HOME)/Library/Caches/tls-render
BUILD_MS  := $(BUILDROOT)/ms

.PHONY: all clean build-clean \
        ms supp \
        ms-html ms-docx ms-pdf \
        supp-html supp-docx supp-pdf \
        sync-sources

all: ms supp

# Sync sources to BUILDROOT. `bib/`, `output/`, `data/`, and `pics/` are
# symlinked from the project so that:
#   - qmd files using `../bib/...` YAML paths resolve to BUILDROOT/bib
#   - `here::here("output", "models")` from chunks evaluated in BUILD_MS
#     resolves to the project's `output/models/` (brms-cached fits)
#   - `here::here("data", ...)` resolves to the project's `data/` (shrimp xlsx)
#   - `here::here("pics", ...)` resolves to the project's `pics/` (Figure 4
#     species illustrations: danio life stages, shrimp, euca, dros, etc.)
#   - `here::here("inst", "extdata", ...)` resolves to the project's `inst/`
#     (Figure 5 reads the Rezende 2020 Chilean hourly temperature CSV)
#
# R/ symlinks are intentionally NOT created here: the analytical functions
# now live in the bayesTLS R package, which the supplement loads via
# `library(bayesTLS)`. Install the package once before rendering:
#   remotes::install_github("daniel1noble/bayesTLS")
#
# Symlinks at BUILDROOT serve the YAML `../bib/...` style; symlinks at
# BUILD_MS serve `here::here(...)` (which detects no project root in the
# build dir and falls back to cwd = BUILD_MS).
sync-sources:
	@mkdir -p $(BUILD_MS) $(OUTDIR)
	rsync -a --delete \
	  --exclude '*.pdf' --exclude '*.html' --exclude '*.docx' \
	  --exclude '*.tex' --exclude '*.log' \
	  --exclude '*_files/' --exclude '*_cache/' --exclude '.quarto/' \
	  ms/ $(BUILD_MS)/
	@ln -sfn $(abspath bib)    $(BUILDROOT)/bib
	@ln -sfn $(abspath output) $(BUILDROOT)/output
	@ln -sfn $(abspath data)   $(BUILDROOT)/data
	@ln -sfn $(abspath pics)   $(BUILDROOT)/pics
	@ln -sfn $(abspath inst)   $(BUILDROOT)/inst
	@ln -sfn $(abspath bib)    $(BUILD_MS)/bib
	@ln -sfn $(abspath output) $(BUILD_MS)/output
	@ln -sfn $(abspath data)   $(BUILD_MS)/data
	@ln -sfn $(abspath pics)   $(BUILD_MS)/pics
	@ln -sfn $(abspath inst)   $(BUILD_MS)/inst

# Render <source.qmd> to <format> with output filename <basename.ext>.
# All work happens in $(BUILD_MS); only the final artefact is copied to
# the project's $(OUTDIR).
# $(4) is an optional extra-arg slot for --metadata-file overrides.
define RENDER
	$(MAKE) sync-sources
	cd $(BUILD_MS) && quarto render $(notdir $(1)) --to $(2) --output $(3) $(4)
	cp $(BUILD_MS)/$(3) $(OUTDIR)/$(3)
endef

# Standalone supp targets need "S" labels for figures/tables/equations
# (Figure S1, Table S1, …). HTML/DOCX use crossref prefix overrides in
# `_supp-overrides.yml`; PDF instead uses the supp-labels Lua filter
# (`_supp-pdf-overrides.yml`) which renames LaTeX counters so numbers
# render as "S1, S2, …". Keeping these out of supplement.qmd's own YAML
# keeps them from polluting any other render that includes its content.
SUPP_META_HTMLDOCX := --metadata-file _supp-overrides.yml
SUPP_META_PDF      := --metadata-file _supp-pdf-overrides.yml

# DOCX-specific (ms only): Quarto auto-renders the YAML `abstract:` right
# after the title, putting it BEFORE the body content. To get the order
# title → authors+affiliations → abstract → keywords (matching the PDF),
# we suppress the YAML abstract for DOCX with `--metadata abstract=""`
# and provide the abstract inside the content-visible-docx block.
DOCX_META := --metadata abstract="" --metadata abstract-title=""

# ---- ms.qmd (manuscript) --------------------------------------------------
ms: ms-html ms-docx ms-pdf
ms-html: ; $(call RENDER,ms/ms.qmd,html,ms.html)
ms-docx: ; $(call RENDER,ms/ms.qmd,docx,ms.docx,$(DOCX_META))
ms-pdf:  ; $(call RENDER,ms/ms.qmd,pdf,ms.pdf)

# ---- supplement.qmd (supplement) ------------------------------------------
supp: supp-html supp-docx supp-pdf

supp-html:
	$(call RENDER,ms/supplement.qmd,html,supp.html,$(SUPP_META_HTMLDOCX))
	$(call STRIP_SUPP_NBSP_HTML,$(OUTDIR)/supp.html)
	# Publish the supplement as the GitHub Pages landing page. supp.html is
	# self-contained (embed-resources), so the root copy needs no asset folder.
	cp $(OUTDIR)/supp.html index.html

supp-docx:
	$(call RENDER,ms/supplement.qmd,docx,supp.docx,$(SUPP_META_HTMLDOCX))
	$(call STRIP_SUPP_NBSP_DOCX,$(OUTDIR)/supp.docx)

supp-pdf:
	$(call RENDER,ms/supplement.qmd,pdf,supp.pdf,$(SUPP_META_PDF))

# Strip the non-breaking space (U+00A0, UTF-8 0xc2 0xa0) Quarto inserts
# between the "Figure S" / "Table S" / "Equation S" prefix and the
# number that follows. Pandoc's HTML writer hardcodes this nbsp; there's
# no Quarto option to suppress it for the standard fig/tbl/eq crossref
# types. We strip it post-render so the supplement reads "Figure S1"
# (the standard journal supplement format), not "Figure S 1".
define STRIP_SUPP_NBSP_HTML
	@perl -i -pe 'BEGIN{binmode STDIN,":raw";binmode STDOUT,":raw"} \
	  s/(Figure S)\xc2\xa0/$$1/g; \
	  s/(Table S)\xc2\xa0/$$1/g; \
	  s/(Equation S)\xc2\xa0/$$1/g' $(1)
endef

# DOCX is a zip; the relevant text lives in `word/document.xml`. The
# helper script extracts, strips nbsp between "Figure S"/"Table S"/
# "Equation S" and the digit that follows, then repackages the archive
# in place.
define STRIP_SUPP_NBSP_DOCX
	@python3 scripts/strip-docx-nbsp.py $(1)
endef

# ---- maintenance ----------------------------------------------------------

# Clean project-side outputs only. Build caches in $(BUILDROOT) are kept
# so the next render can reuse them.
clean:
	rm -rf $(OUTDIR)
	@find ms -maxdepth 1 \( -name '*.pdf' -o -name '*.html' -o -name '*.docx' \
	  -o -name '*.tex' -o -name '*.log' \) -delete 2>/dev/null || true
	@find . -maxdepth 2 -name "*conflicted copy*" -delete 2>/dev/null || true

# Wipe the entire out-of-tree build directory (caches and intermediates).
build-clean:
	rm -rf $(BUILDROOT)
