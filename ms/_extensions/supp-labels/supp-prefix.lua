--[[
  supp-prefix.lua

  PDF/LaTeX-only filter. Without it, the supplement's PDF renders figure
  captions as "Figure S~\arabic{figure}" (where ~ is a non-breaking
  space) — the prefix override "Figure S" plus the hardcoded nbsp Pandoc
  uses between prefix and number gives a visible "Figure S 1".

  This filter injects raw LaTeX at the start of the document body to:
    - Reset the figure/table/equation counters to 0
    - Redefine `\thefigure` (etc.) to print "S<arabic>" so the rendered
      number itself becomes "S1, S2, …"

  Combined with `\figurename{Figure}` (the LaTeX default — NOT the
  prefix override), captions read "Figure S1" with a single non-breaking
  space (visually a normal space) — the standard journal supplement
  format.

  The filter does nothing for HTML or DOCX. For HTML, the Makefile
  post-processes the rendered file with sed to strip the nbsp. For DOCX,
  the nbsp renders as a normal-looking space and is left as-is.
--]]

if not (FORMAT and FORMAT:match('latex')) then
  return {}
end

local LATEX_RESET = [[
\setcounter{figure}{0}
\renewcommand{\thefigure}{S\arabic{figure}}
\setcounter{table}{0}
\renewcommand{\thetable}{S\arabic{table}}
\setcounter{equation}{0}
\renewcommand{\theequation}{S\arabic{equation}}
]]

function Pandoc(doc)
  doc.blocks:insert(1, pandoc.RawBlock('latex', LATEX_RESET))
  return doc
end
