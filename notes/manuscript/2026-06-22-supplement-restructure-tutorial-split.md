# Plan: split the tutorial out of the supplement; should the SI become a Book?

**Status:** planning / for discussion (before the CTmax/z reconfiguration work).
**Date:** 2026-06-22.

## The problem

`ms/supplement.qmd` is **6,901 lines** and does *triple duty*:

1. **Journal SI** for Ecology Letters (must submit as one PDF/DOCX).
2. **Package tutorial** (how to use `bayesTLS`: function tour, simulation
   walkthrough, worked case studies, heat-injury demo).
3. **The "online vignette"** — the GitHub Pages site at
   `daniel1noble.github.io/bayesTLS` is currently just `supp.html` copied to
   `index.html` by the Makefile. There is **no pkgdown site and no `vignettes/`**
   in the repo today, despite the README advertising an "online vignette".

It is hard to navigate and conflates teaching material with results that support
the paper's claims.

## Two pieces of history/context that constrain the options

- **Book mode was already tried and abandoned.** `CLAUDE.md` §8a:
  *"There is intentionally NO project-level `_quarto.yml`… This avoids the TOC
  pollution and format conflicts that book mode caused."* So "turn the supplement
  into a Book" partially **reverts a deliberate decision**. Before going back to a
  book we need to know exactly which conflicts bit last time (TOC handling, the
  per-doc YAML, the "S"-label crossref overrides, the authblk title page, the
  tinytable PDF preamble, the out-of-tree Dropbox build) and confirm a book
  re-solves rather than reopens them.
- **The current build is tuned for two standalone docs** (`ms.qmd`,
  `supplement.qmd`): out-of-tree render to `~/Library/Caches/tls-render`,
  `--metadata-file _supp-overrides.yml` for "Figure S1" labels, hardcoded
  `authblk` title pages, tinytable dependency registration, brms `file=` caches in
  `output/models`. Any restructure has to carry these along.

## Content map (what is tutorial vs SI)

Rough split of the current supplement:

- **Tutorial (package usage), ≈ lines 165–1585:** Introduction; *What each
  function does* (`#sec-function-tour`); *A Tutorial with Simulations* (simulate
  data → two-stage pipeline `#sec-twostage` → joint 4PL `#sec-joint` → derive
  z/CTmax/T_crit `#sec-joint-derive` → compare pipelines `#sec-compare` → worked
  temperature-on-shape example `#sec-joint-extend`); heat-injury accumulation +
  validation (`#sec-hi-validation`).
- **SI proper (supports the paper), ≈ 1586 →:** *Extended Simulation Results*
  (`#sec-extendsims`: equivalence baseline, likelihood misspecification,
  asymptote drift, coverage sweeps, experimental-design impacts); the empirical
  **case studies** (shrimp, zebrafish, snow gum, *Drosophila*) and heat-injury
  forecasts the manuscript cites; derivations.
- **Overlap (needs a decision):** the case studies are *both* results (SI) and
  great worked examples (tutorial). Also the new
  `notes/2026-06-16-direct-ctmax-z-parameterisation.qmd` tutorial could become a
  tutorial article.

## Options

### Option 1 (recommended) — Tutorial → pkgdown articles; SI stays standalone
- Stand up **pkgdown** for the package and move the tutorial into a few
  **articles** under `vignettes/articles/` (pkgdown builds these; `R CMD check` /
  CRAN do **not**, so the heavy cached brms fits are fine). pkgdown gives native
  multi-page navigation (navbar + per-page sidebar) **and** integrates the
  function reference — exactly the "easier to navigate" goal, and it ships with
  the package.
- The SI (`supplement.qmd`) keeps its **standalone single-file** architecture
  (journal PDF/DOCX), now much smaller. No `_quarto.yml`, no book mode → preserves
  the decision that fixed the earlier problems.
- pkgdown becomes the real owner of `daniel1noble.github.io/bayesTLS`, replacing
  the `supp.html`→`index.html` hack and making the README's "online vignette"
  claim true.
- **Pros:** idiomatic R-package docs; navigation solved; SI stays journal-ready;
  avoids reintroducing book mode; fixes the site/vignette gap.
- **Cons:** set up pkgdown + a deploy workflow; decide where the case studies go;
  cross-references that used Quarto `@` within the supplement become within-article
  refs.

### Option 2 — SI itself becomes a Quarto Book
- `_quarto.yml` with `project: book`; chapters = separate `.qmd`s; multi-page HTML
  nav + a single combined PDF/DOCX.
- **Pros:** clean chapter separation; nice HTML nav for the SI.
- **Cons:** reintroduces the `_quarto.yml`/book mode that caused TOC pollution +
  format conflicts (must re-solve); the S-labels, authblk title page, tinytable
  PDF deps and out-of-tree build all need re-engineering for book mode; a journal
  SI is submitted as one PDF anyway, so multi-page HTML nav buys little for the
  *SI*. Higher risk/effort.

### Option 3 — Hybrid
- Tutorial → pkgdown articles (Option 1), **and** if the remaining SI is still
  unwieldy, split it into a couple of standalone `.qmd`s (e.g.
  `supplement-simulations.qmd`, `supplement-casestudies.qmd`) that each render
  independently — modular files without adopting book mode.

## Recommendation

**Option 1 (escalate to Option 3 if the trimmed SI is still too big).** pkgdown is
the low-risk, idiomatic way to give the tutorial real navigation and ship it with
the package, and it directly fixes the missing-vignette situation. Reserve a
Quarto **book** only if we specifically want the *SI itself* as a navigable
multi-page website — which a journal SI doesn't really need, and which reopens the
book-mode problems we deliberately left behind.

## Decisions needed from Daniel (to discuss)

1. **Tutorial home:** pkgdown articles (recommended) vs a standalone Quarto book
   vs CRAN vignettes (CRAN = heavy-fit build problem).
2. **Case studies:** tutorial only, SI only, or both (SI = concise results,
   tutorial = full worked workflow)?
3. **Confirm** the SI must remain a single submittable file for Ecology Letters
   (implies *not* book mode for the SI).
4. Should the **direct-CTmax/z note** and the **heat-injury demo** become tutorial
   articles too?
5. **Site ownership:** stand up pkgdown to own `daniel1noble.github.io/bayesTLS`
   (so the SI HTML is just a render artifact, not the Pages landing page)?

## Effort / risk notes

- pkgdown setup: low–moderate (`usethis::use_pkgdown()`, `_pkgdown.yml`, a deploy
  GitHub Action). Content move = mostly relocating existing chunks + fixing
  cross-refs.
- Heavy fits: articles reuse the existing `output/models` `file=` caches (or
  precompute); `vignettes/articles/` keeps them out of CRAN's build.
- Cross-document refs: stay plain-text across ms ↔ SI ↔ tutorial (existing
  convention); more of them to manage after the split.
- Whatever we choose, the S-labels / authblk / tinytable / out-of-tree-build
  machinery must be preserved for the SI.
