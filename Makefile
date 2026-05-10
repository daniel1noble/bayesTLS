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
# specific to the build dir). `bib/` and `output/` are symlinked from the
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
# Explicit --output is needed for ms_supp.qmd because supplement.qmd's
# `output-file: supp` YAML leaks through the {{< include >}} shortcode and
# would otherwise name the combined doc "supp".
# $(4) is an optional --metadata-file argument for supp targets only.
define RENDER
	$(MAKE) sync-sources
	cd $(BUILD_MS) && quarto render $(notdir $(1)) --to $(2) --output $(3) $(4)
	cp $(BUILD_MS)/$(3) $(OUTDIR)/$(3)
endef

# Standalone supp targets need "S" labels for figures/tables/equations.
# All three formats (HTML/DOCX/PDF) get them via crossref prefix overrides
# in `_supp-overrides.yml`. The Lua filter `supp-labels` is reserved for
# the combined ms_supp render only (where the supp section needs S labels
# but the manuscript section must not). The override is kept out of
# supplement.qmd's own YAML to prevent it leaking into the combined doc
# via the `{{< include >}}` shortcode.
SUPP_META := --metadata-file _supp-overrides.yml

# DOCX-specific: Quarto auto-renders the YAML `abstract:` right after the
# title, putting it BEFORE the body content (where our authors+affiliations
# block lives). To get the order title → authors+affil → abstract → keywords
# (matching the PDF), we suppress the YAML abstract for DOCX with
# `--metadata abstract=""` and provide the abstract inside our content-visible
# block in the body. PDF and HTML still use the YAML abstract.
DOCX_META := --metadata abstract="" --metadata abstract-title=""

# For the combined ms_supp doc: supplement.qmd's `title:` and `subtitle:`
# leak into the wrapper's metadata via `{{< include >}}`. Override at the
# command line so the combined doc shows the manuscript title.
COMBINED_TITLE := A flexible modelling framework for estimating thermal sensitivity across life
COMBINED_META  := --metadata title="$(COMBINED_TITLE)" --metadata subtitle=""

# ---- ms.qmd (manuscript only) ---------------------------------------------
ms: ms-html ms-docx ms-pdf
ms-html: ; $(call RENDER,ms/ms.qmd,html,ms.html)
ms-docx: ; $(call RENDER,ms/ms.qmd,docx,ms.docx,$(DOCX_META))
ms-pdf:  ; $(call RENDER,ms/ms.qmd,pdf,ms.pdf)

# ---- supplement.qmd (supplement only) -------------------------------------
supp: supp-html supp-docx supp-pdf
supp-html: ; $(call RENDER,ms/supplement.qmd,html,supp.html,$(SUPP_META))
supp-docx: ; $(call RENDER,ms/supplement.qmd,docx,supp.docx,$(SUPP_META))
supp-pdf:  ; $(call RENDER,ms/supplement.qmd,pdf,supp.pdf,$(SUPP_META))

# ---- ms_supp.qmd (combined) -----------------------------------------------
ms-supp: ms-supp-html ms-supp-docx ms-supp-pdf
ms-supp-html: ; $(call RENDER,ms/ms_supp.qmd,html,ms_supp.html,$(COMBINED_META))
ms-supp-docx: ; $(call RENDER,ms/ms_supp.qmd,docx,ms_supp.docx,$(DOCX_META) $(COMBINED_META))
ms-supp-pdf:  ; $(call RENDER,ms/ms_supp.qmd,pdf,ms_supp.pdf,$(COMBINED_META))

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
