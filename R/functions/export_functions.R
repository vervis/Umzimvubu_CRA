# ============================================================================
# R/functions/export_functions.R
#
# Stage-80 export helpers. Only function definitions: tar_source("R/functions")
# loads them.
#
# export_pages(): publishes a saveWidget() map (non-selfcontained) to the
#   GitHub Pages folder (default docs/):
#   1. Copies <name>.html to docs/index.html, rewriting references to the
#      dependency folder <name>_files/ as index_files/.
#   2. Replaces docs/index_files/ with a fresh copy of <name>_files/.
#   3. Ensures docs/.nojekyll exists (so Pages serves the _files folder).
#   Returns all written paths, for a format = "file" target.
# ============================================================================

export_pages <- function(html, dest = "docs", page = "index") {
  src_lib <- sub("\\.html$", "_files", html)
  if (!file.exists(html) || !dir.exists(src_lib))
    stop("Missing map or its dependency folder: ", html, " / ", src_lib)

  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  out <- file.path(dest, paste0(page, ".html"))
  lib <- file.path(dest, paste0(page, "_files"))

  txt <- readChar(html, file.size(html), useBytes = TRUE)
  txt <- gsub(paste0(basename(src_lib), "/"), paste0(basename(lib), "/"), txt,
              fixed = TRUE, useBytes = TRUE)
  writeChar(txt, out, eos = NULL, useBytes = TRUE)

  unlink(lib, recursive = TRUE)
  dir.create(lib)
  ok <- file.copy(list.files(src_lib, full.names = TRUE), lib, recursive = TRUE)
  if (!all(ok)) stop("Copying ", src_lib, " to ", lib, " failed")

  nj <- file.path(dest, ".nojekyll")
  if (!file.exists(nj)) file.create(nj)

  c(out, nj, list.files(lib, recursive = TRUE, full.names = TRUE))
}
