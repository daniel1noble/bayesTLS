.onAttach <- function(libname, pkgname) {
  v <- utils::packageVersion("bayesTLS")
  packageStartupMessage(
    "bayesTLS ", v, "\n",
    "Please cite: Noble DWA, Arnold PA, Nakagawa S & Pottier P (2026) A flexible modelling\n",
    "  framework for estimating thermal tolerance and sensitivity. Manuscript in\n",
    "  preparation.\n",
    "bayesTLS fits with brms and Stan -- please also cite brms (Burkner 2017) and the Stan\n",
    "  backend you use (RStan or cmdstanr). Run  citation(\"bayesTLS\")  for all entries.\n\n",
    "Tutorial & online vignette: https://daniel1noble.github.io/bayesTLS/"
  )
}
