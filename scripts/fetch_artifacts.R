#!/usr/bin/env Rscript
# Fetch cached artifacts from the project's OSF deposit (node c6dxy) and unpack
# them under output/, so a fresh clone can render the paper / reproduce results
# without re-fitting models or re-running the simulation.
#
#   Rscript scripts/fetch_artifacts.R           # default tiers: results + models (~85 MB)
#   Rscript scripts/fetch_artifacts.R full      # ALSO the 2.6 GB raw/draws re-derivation tier
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

# Tier -> zip name + a marker glob (output-relative). A tier counts as already
# present if the marker matches any local file. Keep in sync with osf_publish.R.
TIERS <- list(
  results = list(zip = "results.zip",  marker = "sim_twostage/per_sim_*.rds"),
  models  = list(zip = "models.zip",   marker = "models/*.rds"),
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
have_it <- function(marker) length(Sys.glob(file.path(OUT, marker))) > 0

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
