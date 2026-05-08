--[[
  supp-prefix.lua

  Reset and rename figure/table/equation counters in LaTeX/PDF output for
  content inside a `:::{.supp-section}` div, so the supplement section
  prints "Figure S1, S2, …", "Table S1, …", "Equation (S1), …".

  Used by:
    - supplement.qmd (standalone supp.pdf) — body wrapped in one outer
      `.supp-section` div so the entire doc gets S-numbered counters.
    - ms_supp.qmd (combined ms+supp) — only the supplement include is
      wrapped, so the manuscript section keeps regular labels and the
      supplement section gets S labels with numbering reset to S1.

  HTML and DOCX get "S" labels via Quarto's native `crossref:` prefix
  overrides, applied for those formats through `--metadata-file
  ms/_supp-overrides.yml` for standalone supp renders, and accepted as
  continuous numbering in the combined ms_supp HTML/DOCX (HTML is
  primarily a navigation surface; the formal "S" deliverable is the PDF).

  The filter only injects counters once on the first .supp-section div
  encountered. Nested wrappers (e.g. supp.qmd's own .supp-section being
  included inside ms_supp.qmd's .supp-section) are harmless — only the
  outermost transition resets counters.
--]]

local is_latex = FORMAT and FORMAT:match('latex')
local injected = false

local LATEX_RESET = [[
\setcounter{figure}{0}
\renewcommand{\thefigure}{S\arabic{figure}}
\setcounter{table}{0}
\renewcommand{\thetable}{S\arabic{table}}
\setcounter{equation}{0}
\renewcommand{\theequation}{S\arabic{equation}}
]]

function Div(el)
  if not is_latex then return nil end
  if injected then return nil end
  if not (el.classes and el.classes:includes('supp-section')) then return nil end
  injected = true
  local new_content = pandoc.List({ pandoc.RawBlock('latex', LATEX_RESET) })
  for _, c in ipairs(el.content) do
    new_content:insert(c)
  end
  el.content = new_content
  return el
end
