# Shared setup for the article pages: project-local R library, and the HumanDLPFC h5ad
# export the Python tabs read (created once; get-started.qmd uses the same file).
local({
  root <- Sys.getenv("QUARTO_PROJECT_DIR", unset = normalizePath(".."))
  if (dir.exists(file.path(root, "_rlib"))) .libPaths(c(file.path(root, "_rlib"), .libPaths()))
  reticulate::py_run_file(file.path(root, "scripts", "py-logging-shim.py"))
  h5ad <- file.path(root, "_cache", "HumanDLPFC.h5ad")
  if (!file.exists(h5ad)) {
    old <- setwd(root)
    on.exit(setwd(old))
    status <- system2(file.path(R.home("bin"), "Rscript"),
                      c("scripts/export-h5ad.R", "_cache/HumanDLPFC.h5ad"))
    stopifnot(status == 0)
  }
})
