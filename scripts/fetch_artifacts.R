#!/usr/bin/env Rscript
# Fetch cached artifacts from the project's OSF deposit (node c6dxy) and unpack
# them under output/, so a fresh clone can render the paper / reproduce results
# without re-fitting models or re-running the simulation.
#
#   Rscript scripts/fetch_artifacts.R           # default tiers: results + models (~205 MB)
#   Rscript scripts/fetch_artifacts.R full      # only the 2.6 GB raw/draws re-derivation tier
#                                               # (make data-full fetches all three tiers)
#   Rscript scripts/fetch_artifacts.R results   # one named tier
#   Rscript scripts/fetch_artifacts.R --force   # re-download even if already present
#
# Check-then-download: a tier whose files are already present locally is skipped
# (override with --force). No token needed -- c6dxy is public. Each download is
# verified against the md5 recorded in data-raw/osf_manifest.json when present.
# The zips store output-relative paths (sim_twostage/..., data/..., models/...),
# so they unpack straight into output/.

suppressPackageStartupMessages({
  library(curl); library(httr); library(jsonlite); library(tools)
})

NODE     <- "c6dxy"
OUT      <- here::here("output")
MANIFEST <- here::here("data-raw", "osf_manifest.json")

# Tier -> zip name + marker globs (output-relative). A tier counts as already
# present only if EVERY marker matches a local file. The models markers are the
# fits ms/ms.qmd and ms/supplement.qmd read that only models.zip provides: git
# tracks a few small .rds files in output/models/, so the old "models/*.rds"
# marker matched on a fresh clone and the fits were never downloaded. Add a fit
# here when a document starts reading it. Keep tiers in sync with osf_publish.R.
MODEL_FITS <- c(
  "aphid_species", "aphid_species_age",
  "fit_dsuzukii_coma_4pl", "fit_dsuzukii_coma_linear",
  "fit_dsuzukii_mortality_4pl_F", "fit_dsuzukii_mortality_4pl_M",
  "fit_dsuzukii_mortality_sex_direct", "fit_dsuzukii_mortality_sex_joint",
  "fit_dsuzukii_prod_hurdle_gamma", "fit_dsuzukii_prod_hurdle_lognormal",
  "fit_leaf_recovery_4pl_cens", "fit_leaf_recovery_4pl_cens_run_array",
  "sim_4pl_direct", "zf_oxygen_normhyper",
  # alternative dose-response curves (supplement curve-comparison section)
  paste0(rep(c("zf_oxygen_normhyper", "aphid_species",
               "fit_dsuzukii_mortality_4pl_F", "fit_dsuzukii_mortality_4pl_M",
               "fit_leaf_recovery_4pl_cens_run_array"), each = 3),
         c("_lognormal", "_weibull1", "_weibull2")),
  # prior-sensitivity refits (supplement prior-sensitivity section)
  paste0(rep(c("zf_oxygen_normhyper", "aphid_species"), each = 2), c("_prior_weak", "_prior_strong")),
  # 2PL refits with both asymptotes fixed at 0 and 1 (supplement 2PL comparison)
  paste0(c("zf_oxygen_normhyper", "aphid_species", "fit_dsuzukii_mortality_4pl_F",
           "fit_dsuzukii_mortality_4pl_M"), "_fixed01"),
  # sublethal fly 4PLs behind manuscript Figure 6 (supplement Section 10)
  "fit_dsuzukii_coma_4pl_temp", "fit_dsuzukii_coma_4pl_tempsex",
  "fit_dsuzukii_fertility_4pl", "fit_dsuzukii_fertility_4pl_temp", "fit_dsuzukii_fertility_4pl_tempsex"
)
TIERS <- list(
  results = list(zip = "results.zip",  marker = "sim_twostage/per_sim_*.rds"),
  models  = list(zip = "models.zip",
                 marker = file.path("models", paste0(MODEL_FITS, ".rds"))),
  full    = list(zip = "sim_full.zip", marker = "sim_twostage/draws_*.rds")
)
DEFAULT <- c("results", "models")

# ---- args ------------------------------------------------------------------
argv  <- commandArgs(trailingOnly = TRUE)
force <- "--force" %in% argv
tiers <- setdiff(argv, "--force")
if (!length(tiers)) tiers <- DEFAULT
bad <- setdiff(tiers, names(TIERS))
if (length(bad)) stop("unknown tier(s): ", paste(bad, collapse = ", "),
                      " (have: ", paste(names(TIERS), collapse = ", "), ")",
                      call. = FALSE)

# ---- helpers ---------------------------------------------------------------
have_it <- function(marker)
  all(vapply(marker, function(m) length(Sys.glob(file.path(OUT, m))) > 0, logical(1)))

read_manifest <- function() if (file.exists(MANIFEST))
  fromJSON(MANIFEST, simplifyVector = FALSE) else list(tiers = list())

osf_find_file <- function(name) {
  r <- GET(sprintf("https://api.osf.io/v2/nodes/%s/files/osfstorage/", NODE))
  stop_for_status(r)
  d <- fromJSON(content(r, as = "text", encoding = "UTF-8"),
                simplifyVector = FALSE)$data
  for (f in d) if (f$attributes$name == name) return(f)
  NULL
}

# ---- main ------------------------------------------------------------------
man <- read_manifest()
for (t in tiers) {
  spec <- TIERS[[t]]
  if (!force && have_it(spec$marker)) {
    message(sprintf("[%s] already present - skip (use --force to re-download)", t))
    next
  }
  f <- osf_find_file(spec$zip)
  if (is.null(f)) {
    message(sprintf("[%s] %s is not on OSF yet - skip", t, spec$zip))
    next
  }
  exp_md5 <- man$tiers[[t]]$zip_md5
  tmp <- tempfile(fileext = ".zip")
  message(sprintf("[%s] downloading %s (%.1f MB)...", t, spec$zip,
                  f$attributes$size / 1e6))
  curl_download(f$links$download, tmp)                 # public: no auth header
  if (!is.null(exp_md5)) {
    got <- unname(md5sum(tmp))
    if (!identical(got, exp_md5))
      stop(sprintf("[%s] md5 mismatch: got %s, manifest expects %s. Aborting.",
                   t, got, exp_md5), call. = FALSE)
    message(sprintf("[%s] md5 verified", t))
  } else {
    message(sprintf("[%s] no manifest md5 recorded - skipping integrity check", t))
  }
  dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
  utils::unzip(tmp, exdir = OUT)
  unlink(tmp)
  message(sprintf("[%s] unpacked into %s", t, OUT))
}
message("done.")
