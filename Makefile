# Project Makefile.
#
# Three documents, each rendered to HTML, DOCX, and PDF:
#
#   ms       ms/ms.qmd          → _output/ms.{html,docx,pdf}      (manuscript only)
#   supp     ms/supplement.qmd  → _output/supp.{html,docx,pdf}    (supplement only)
#   ms-supp  ms/ms_supp.qmd     → _output/ms_supp.{html,docx,pdf} (combined ms + supp)
#
# `ms_supp.qmd` is a thin wrapper that uses Quarto `{{< include >}}` to pull in
# the manuscript and the supplement. The supplement section is wrapped in
# `:::{.supp-section}` so the local Lua filter `ms/_extensions/supp-labels/`
# applies "S" prefixes to figure/table/equation numbering in that section only.
#
# There is intentionally NO project-level `_quarto.yml`. Each .qmd carries its
# own complete YAML, so renders are fully independent.
#
# Dropbox: this project lives inside Dropbox. To avoid sync conflicts during
# rendering, intermediate output goes to `$(TMPDIR)`, which is marked
# `com.dropbox.ignored=1`. Final outputs are then moved to `$(OUTDIR)/`,
# which IS synced. Run `make dropbox-ignore` once after cloning to apply
# the xattr to all build dirs.

OUTDIR := _output
TMPDIR := _tmp

.PHONY: all clean dropbox-ignore \
        ms supp ms-supp \
        ms-html ms-docx ms-pdf \
        supp-html supp-docx supp-pdf \
        ms-supp-html ms-supp-docx ms-supp-pdf

all: ms supp ms-supp

# Render <source.qmd> to <format> with output filename <basename.ext>,
# directing Quarto's intermediate write to $(TMPDIR) (Dropbox-ignored), then
# move to $(OUTDIR)/.
# Explicit --output is needed for ms_supp.qmd because supplement.qmd's
# `output-file: supp` YAML leaks through the {{< include >}} shortcode and
# would otherwise name the combined doc "supp".
# $(4) is an optional --metadata-file argument for supp targets only.
define RENDER
	mkdir -p $(OUTDIR) $(TMPDIR)
	@xattr -w com.dropbox.ignored 1 $(TMPDIR) 2>/dev/null || true
	quarto render $(1) --to $(2) --output $(3) --output-dir $(abspath $(TMPDIR)) $(4)
	mv $(TMPDIR)/$(3) $(OUTDIR)/$(3)
	@if [ -d $(TMPDIR)/$(basename $(notdir $(3)))_files ]; then \
	  rm -rf $(TMPDIR)/$(basename $(notdir $(3)))_files; \
	fi
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
clean:
	rm -rf $(OUTDIR) $(TMPDIR) ms/*_files ms/*_cache ms/.quarto ms/*.tex ms/*.log

# Tell Dropbox to skip syncing build artefacts. Run once after cloning.
# `_output/` is intentionally NOT included — final outputs should sync.
dropbox-ignore:
	@mkdir -p $(TMPDIR)
	@for d in $(TMPDIR) ms/ms_cache ms/supplement_cache ms/ms_supp_cache \
	         ms/ms_files ms/supplement_files ms/ms_supp_files \
	         ms/.quarto; do \
	  if [ -e "$$d" ]; then \
	    xattr -w com.dropbox.ignored 1 "$$d" && echo "ignored: $$d"; \
	  fi; \
	done
