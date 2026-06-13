.onAttach <- function(libname, pkgname) {
  v <- utils::packageVersion("bayesTLS")
  packageStartupMessage(
    "bayesTLS ", v, "\n",
    "Please cite: Noble DWA, Arnold PA & Pottier P (2026) A flexible modelling\n",
    "  framework for estimating thermal tolerance and sensitivity. Manuscript in\n",
    "  preparation. Run  citation(\"bayesTLS\")  for the full entry.\n\n",
    "Tutorial & online vignette: https://daniel1noble.github.io/bayesTLS/"
  )
}
