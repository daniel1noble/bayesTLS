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

# Sync sources to BUILDROOT. `bib/` and `output/` are symlinked from the
# project so that:
#   - qmd files using `../bib/...` paths resolve correctly
#   - `here::here("output", "models")` resolves to the project's
#     `output/models/` directory, where brms-cached .rds fits live —
#     otherwise brms re-fits the models every render because it can't
#     find them in the build dir.
sync-sources:
	@mkdir -p $(BUILD_MS) $(OUTDIR)
	rsync -a --delete \
	  --exclude '*.pdf' --exclude '*.html' --exclude '*.docx' \
	  --exclude '*.tex' --exclude '*.log' \
	  --exclude '*_files/' --exclude '*_cache/' --exclude '.quarto/' \
	  ms/ $(BUILD_MS)/
	@ln -sfn $(abspath bib) $(BUILDROOT)/bib
	@ln -sfn $(abspath output) $(BUILDROOT)/output

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
# (Figure S1, Table S1, …). The crossref overrides live in
# `ms/_supp-overrides.yml` and are loaded via --metadata-file. Keeping
# them out of supplement.qmd's own YAML keeps them from polluting the
# manuscript render.
SUPP_META := --metadata-file _supp-overrides.yml

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
supp-html: ; $(call RENDER,ms/supplement.qmd,html,supp.html,$(SUPP_META))
supp-docx: ; $(call RENDER,ms/supplement.qmd,docx,supp.docx,$(SUPP_META))
supp-pdf:  ; $(call RENDER,ms/supplement.qmd,pdf,supp.pdf,$(SUPP_META))

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
