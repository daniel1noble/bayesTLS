#!/usr/bin/env Rscript
# Publish cached artifacts to the project's OSF deposit (node c6dxy) so OSF
# mirrors the latest local results / models / simulation outputs.
#
#   OSF_PAT=<token> Rscript scripts/osf_publish.R                # changed tiers
#   OSF_PAT=<token> Rscript scripts/osf_publish.R results models # named tiers
#   OSF_PAT=<token> Rscript scripts/osf_publish.R --force        # all, ignore cache
#   Rscript scripts/osf_publish.R --status                       # dry run, no token
#
# Idempotent: a tier is re-zipped + uploaded only when its SOURCE files changed
# since the last publish. "Changed" is a cheap signature over the sorted set of
# (relative path, size, mtime) for the tier's globs -- so re-running after a
# render that touched nothing uploads nothing. The signature, the uploaded zip's
# md5, and the publish time are recorded per tier in data-raw/osf_manifest.json
# (committed; the fetch script reads the same file for integrity checks).
#
# The token must be an OSF personal access token (scope osf.full_write) in the
# OSF_PAT environment variable (e.g. a line `OSF_PAT=...` in ~/.Renviron). It is
# NEVER printed, logged, or written to disk. OSF tokens are account-wide -- treat
# OSF_PAT like a password and keep it out of the repo and any CI log.

suppressPackageStartupMessages({
  library(httr); library(jsonlite); library(tools)
})

NODE     <- "c6dxy"
OUT      <- here::here("output")
STAGE    <- file.path(tempdir(), "osf_stage")
MANIFEST <- here::here("data-raw", "osf_manifest.json")
dir.create(STAGE, showWarnings = FALSE, recursive = TRUE)

# Tier -> zip name + the source globs it bundles (relative to output/). All zips
# are built from inside output/, so their internal paths are output-relative
# (sim_twostage/..., data/..., models/...) and the fetch script unzips into
# output/. Keep these globs in sync with the fetch manifest and the design note
# (notes/2026-06-23-osf-artifact-fetch.qmd).
TIERS <- list(
  results = list(zip = "results.zip", globs = c(
    "sim_twostage/per_sim_*.rds", "sim_twostage/summary_*.rds",
    "sim_twostage/diffs_*.rds",   "sim_twostage/diag_*.rds", "data/*.rds")),
  models  = list(zip = "models.zip", globs = "models/*.rds"),
  full    = list(zip = "sim_full.zip", globs = c(
    "sim_twostage/raw", "sim_twostage/draws_*.rds", "sim_twostage/meta_*.rds"))
)

# ---- args ------------------------------------------------------------------
argv   <- commandArgs(trailingOnly = TRUE)
force  <- "--force"  %in% argv
status <- "--status" %in% argv
tiers  <- setdiff(argv, c("--force", "--status"))
if (!length(tiers)) tiers <- names(TIERS)
bad <- setdiff(tiers, names(TIERS))
if (length(bad)) stop("unknown tier(s): ", paste(bad, collapse = ", "),
                      " (have: ", paste(names(TIERS), collapse = ", "), ")",
                      call. = FALSE)

# ---- helpers ---------------------------------------------------------------
# Expand a tier's globs to actual files (recursing into any directory, e.g. raw/).
tier_files <- function(globs) {
  hits <- unlist(lapply(globs, function(g) Sys.glob(file.path(OUT, g))))
  files <- unlist(lapply(hits, function(p)
    if (dir.exists(p)) list.files(p, recursive = TRUE, full.names = TRUE) else p))
  sort(unique(files))
}

# Cheap content signature: hash of sorted (output-relative path, size, mtime).
# Avoids reading 1.4 GB of raw/ just to decide whether anything changed.
tier_signature <- function(files) {
  if (!length(files)) return(NA_character_)
  info <- file.info(files)
  rel  <- sub(paste0("^", OUT, "/"), "", files)
  digest_str <- paste(rel, info$size, as.numeric(info$mtime), collapse = "\n")
  unname(md5sum(textConnection_to_tmp(digest_str)))
}
# md5sum() needs a file path; hash the signature string via a temp file.
textConnection_to_tmp <- function(s) {
  tf <- tempfile(); writeLines(s, tf); tf
}

read_manifest <- function() if (file.exists(MANIFEST))
  fromJSON(MANIFEST, simplifyVector = FALSE) else list(node = NODE, tiers = list())
write_manifest <- function(m) writeLines(toJSON(m, auto_unbox = TRUE,
                                                pretty = TRUE), MANIFEST)

# OSF API: find a file in the node's osfstorage root by name (public read needs
# no token; pass one when available for private nodes / rate limits).
osf_find_file <- function(name, token = NULL) {
  h <- if (!is.null(token)) add_headers(Authorization = paste("Bearer", token))
       else add_headers()
  r <- GET(sprintf("https://api.osf.io/v2/nodes/%s/files/osfstorage/", NODE), h)
  stop_for_status(r)
  d <- fromJSON(content(r, as = "text", encoding = "UTF-8"),
                simplifyVector = FALSE)$data
  for (f in d) if (f$attributes$name == name) return(f)
  NULL
}

# Create (new file) or update (new version of existing file) on osfstorage.
osf_upload <- function(path, token) {
  name <- basename(path)
  existing <- osf_find_file(name, token)
  url <- if (is.null(existing))
    sprintf("https://files.osf.io/v1/resources/%s/providers/osfstorage/?kind=file&name=%s",
            NODE, utils::URLencode(name, reserved = TRUE))
  else existing$links$upload                       # PUT here -> new version
  r <- PUT(url, add_headers(Authorization = paste("Bearer", token)),
           body = upload_file(path), timeout(60 * 60))
  stop_for_status(r)
  body <- fromJSON(content(r, as = "text", encoding = "UTF-8"),
                   simplifyVector = FALSE)
  list(action = if (is.null(existing)) "created" else "updated",
       md5 = body$data$attributes$extra$hashes$md5)
}

# ---- main ------------------------------------------------------------------
man   <- read_manifest()
token <- Sys.getenv("OSF_PAT", "")
if (!status && !nzchar(token))
  stop("OSF_PAT is not set. Create an osf.full_write token at ",
       "osf.io/settings/tokens and put `OSF_PAT=...` in ~/.Renviron.",
       call. = FALSE)

for (t in tiers) {
  spec  <- TIERS[[t]]
  files <- tier_files(spec$globs)
  if (!length(files)) { message(sprintf("[%s] no source files found - skip", t)); next }
  sig   <- tier_signature(files)
  prev  <- man$tiers[[t]]
  total <- sum(file.info(files)$size)
  unchanged <- !is.null(prev) && identical(prev$signature, sig)

  if (status) {
    message(sprintf("[%s] %d files, %.1f MB, %s", t, length(files),
                    total / 1e6, if (unchanged) "up to date" else "CHANGED"))
    next
  }
  if (unchanged && !force) {
    message(sprintf("[%s] up to date (%d files) - skip", t, length(files)))
    next
  }

  zip_path <- file.path(STAGE, spec$zip)
  if (file.exists(zip_path)) unlink(zip_path)
  message(sprintf("[%s] zipping %d files (%.1f MB)...", t, length(files), total / 1e6))
  rel <- sub(paste0("^", OUT, "/"), "", files)
  old_wd <- setwd(OUT); on.exit(setwd(old_wd), add = TRUE)         # internal paths
  utils::zip(zip_path, rel, flags = "-rqX")                        # -X: no extra attrs
  setwd(old_wd); on.exit()

  message(sprintf("[%s] uploading %s to OSF %s...", t, spec$zip, NODE))
  res <- osf_upload(zip_path, token)
  message(sprintf("[%s] %s (osf md5 %s)", t, res$action, res$md5))

  man$node <- NODE
  man$tiers[[t]] <- list(zip = spec$zip, globs = as.list(spec$globs),
                         signature = sig, zip_md5 = res$md5,
                         n_files = length(files),
                         size_bytes = total)
  write_manifest(man)
}
message("done.")
