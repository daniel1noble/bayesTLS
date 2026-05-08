# Project Makefile.
#
# Three documents, each rendered to HTML, DOCX, and PDF:
#
#   ms       ms/ms.qmd          → _output/ms.{html,docx,pdf}      (manuscript only)
#   supp     ms/supplement.qmd  → _output/supp.{html,docx,pdf}    (supplement only)
#   ms-supp  ms/ms_supp.qmd     → _output/ms_supp.{html,docx,pdf} (combined ms + supp)
#
# `ms_supp.qmd` is a thin wrapper that uses Quarto `{{< include >}}` to combine
# the manuscript and the supplement. The supplement section is wrapped in
# `:::{.supp-section}` so the local Lua filter `ms/_extensions/supp-labels/`
# applies "S" prefixes to figure/table/equation numbers in that section only.
#
# There is intentionally NO project-level `_quarto.yml`. Each .qmd carries its
# own complete YAML, so renders are fully independent.
#
# ---------------------------------------------------------------------------
# Dropbox: this project lives inside a Dropbox-synced directory. Quarto's
# render produces several intermediate files (`*.tex`, `*.log`, `*_files/`,
# `*_cache/`, the .pdf itself) alongside the source. When Dropbox syncs
# those files mid-render, it produces "conflicted copy" duplicates.
#
# To prevent that entirely, the render is done OUTSIDE Dropbox in
# $(BUILDROOT) (default: ~/Library/Caches/tls-render). Sources are rsync'd
# in, `bib/` is symlinked, and only the final output is copied back into
# the Dropbox-synced project's `_output/`. Build caches persist in
# $(BUILDROOT) across renders — a clean of the project does not invalidate
# them. Run `make build-clean` to wipe build caches.
# ---------------------------------------------------------------------------

OUTDIR    := _output
BUILDROOT := $(HOME)/Library/Caches/tls-render
BUILD_MS  := $(BUILDROOT)/ms

.PHONY: all clean build-clean \
        ms supp ms-supp \
        ms-html ms-docx ms-pdf \
        supp-html supp-docx supp-pdf \
        ms-supp-html ms-supp-docx ms-supp-pdf \
        sync-sources

all: ms supp ms-supp

# Sync sources to BUILDROOT (excluding outputs and caches that are
# specific to the build dir). The `bib/` directory is symlinked so that
# qmd files using relative paths like `../bib/tdt_problems.bib` resolve
# correctly. Output and intermediate files (.pdf, .tex, .log, *_files,
# *_cache) inside ms/ are NOT copied — those are produced anew (or read
# from the build dir's own caches).
sync-sources:
	@mkdir -p $(BUILD_MS) $(OUTDIR)
	rsync -a --delete \
	  --exclude '*.pdf' --exclude '*.html' --exclude '*.docx' \
	  --exclude '*.tex' --exclude '*.log' \
	  --exclude '*_files/' --exclude '*_cache/' --exclude '.quarto/' \
	  ms/ $(BUILD_MS)/
	@ln -sfn $(abspath bib) $(BUILDROOT)/bib

# Render <source.qmd> to <format> with output filename <basename.ext>.
# All work happens in $(BUILD_MS); only the final artefact is copied to
# the project's $(OUTDIR).
# Explicit --output is needed for ms_supp.qmd because supplement.qmd's
# `output-file: supp` YAML leaks through the {{< include >}} shortcode and
# would otherwise name the combined doc "supp".
# $(4) is an optional --metadata-file argument for supp targets only.
define RENDER
	$(MAKE) sync-sources
	cd $(BUILD_MS) && quarto render $(notdir $(1)) --to $(2) --output $(3) $(4)
	cp $(BUILD_MS)/$(3) $(OUTDIR)/$(3)
endef

# Standalone supp targets need crossref prefix overrides for HTML/DOCX/PDF
# native S labels. Keeping them in a separate YAML (instead of supplement.qmd's
# own YAML) prevents them from leaking into the combined ms_supp render.
SUPP_META := --metadata-file _supp-overrides.yml

# ---- ms.qmd (manuscript only) ---------------------------------------------
ms: ms-html ms-docx ms-pdf
ms-html: ; $(call RENDER,ms/ms.qmd,html,ms.html)
ms-docx: ; $(call RENDER,ms/ms.qmd,docx,ms.docx)
ms-pdf:  ; $(call RENDER,ms/ms.qmd,pdf,ms.pdf)

# ---- supplement.qmd (supplement only) -------------------------------------
supp: supp-html supp-docx supp-pdf
supp-html: ; $(call RENDER,ms/supplement.qmd,html,supp.html,$(SUPP_META))
supp-docx: ; $(call RENDER,ms/supplement.qmd,docx,supp.docx,$(SUPP_META))
supp-pdf:  ; $(call RENDER,ms/supplement.qmd,pdf,supp.pdf,$(SUPP_META))

# ---- ms_supp.qmd (combined) -----------------------------------------------
ms-supp: ms-supp-html ms-supp-docx ms-supp-pdf
ms-supp-html: ; $(call RENDER,ms/ms_supp.qmd,html,ms_supp.html)
ms-supp-docx: ; $(call RENDER,ms/ms_supp.qmd,docx,ms_supp.docx)
ms-supp-pdf:  ; $(call RENDER,ms/ms_supp.qmd,pdf,ms_supp.pdf)

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
