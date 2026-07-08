#!/usr/bin/env Rscript
# Mechanical check phases for the deploy-bioc skill. Run from the repo root of
# any Bioconductor-track R package.
#
#   Rscript .claude/skills/deploy-bioc/deploy.R <phase>
#
# Phases:
#   document   devtools::document() -- regenerate man/ and NAMESPACE
#   check      devtools::check()    -- R CMD check, log + summary
#   bioccheck  BiocCheck::BiocCheck -- Bioconductor linter, log + summary
#   all        document -> check -> bioccheck (stops if check reports an ERROR)
#
# The package name, vignette, etc. are read from the repo -- nothing here is
# hard-coded to a specific package, so the script works in any Bioc package.
#
# Logs land in $DEPLOY_OUT (default: a tempdir, path printed on start).
# Exit code is non-zero when a phase surfaces blocking problems -- check or
# BiocCheck errors/warnings, or an aborted build -- so the script doubles as a
# gate. Environmental BiocCheck failures (support-site network call, leftover
# output folder) are reported but NOT counted as blocking.

args  <- commandArgs(trailingOnly = TRUE)
phase <- if (length(args)) args[1] else "all"
valid <- c("document", "check", "bioccheck", "all")
if (!phase %in% valid)
  stop("Unknown phase '", phase, "'. Use one of: ", paste(valid, collapse = ", "))

if (!file.exists("DESCRIPTION"))
  stop("No DESCRIPTION in the working directory -- run deploy.R from the package root.")
pkg <- tryCatch(unname(read.dcf("DESCRIPTION")[, "Package"]),
                error = function(e) NA_character_)

out <- Sys.getenv("DEPLOY_OUT", unset = tempfile("deploy-"))
dir.create(out, showWarnings = FALSE, recursive = TRUE)
message("== deploy.R pkg=", pkg, " phase=", phase, "  logs -> ", out)

fail <- 0L  # accumulates blocking problems across phases

run_document <- function() {
  message("\n== document ==")
  suppressPackageStartupMessages(devtools::document())
  message("document: OK (man/ and NAMESPACE regenerated -- review the diff)")
}

run_check <- function() {
  message("\n== check ==")
  log <- file.path(out, "check.log")
  # error_on = "never": for check *findings*, return the result so we can
  # summarise + fix. But a broken vignette/build fails inside pkgbuild::build()
  # *before* check runs, and that throws regardless of error_on -- catch it so
  # we report a clean blocking error instead of a raw backtrace.
  # quiet = FALSE so the full R CMD check output -- crucially the vignette/build
  # error when the build aborts -- streams to stdout and lands in the redirect.
  res <- tryCatch(
    devtools::check(error_on = "never", quiet = FALSE),
    error = function(e) e)
  if (inherits(res, "condition")) {
    msg <- conditionMessage(res)
    writeLines(c("BUILD/CHECK ABORTED:", msg), log)
    message("check: build aborted (usually a vignette or install failure)")
    message("  ", gsub("\n", "\n  ", msg))
    message("  full log: ", log)
    fail <<- fail + 1L
    return(invisible(FALSE))
  }
  writeLines(c(res$errors, res$warnings, res$notes), log)
  ne <- length(res$errors); nw <- length(res$warnings); nn <- length(res$notes)
  message(sprintf("check: %d error(s) | %d warning(s) | %d note(s)  (%s)",
                  ne, nw, nn, log))
  for (blk in c(res$errors, res$warnings, res$notes)) {
    message("----"); message(blk)
  }
  # BiocCheck requires a package that passes check cleanly: block on warnings too.
  if (ne + nw > 0L) fail <<- fail + ne + nw
  invisible(ne + nw == 0L)
}

run_bioccheck <- function() {
  message("\n== bioccheck ==")
  # BiocCheck writes a "<pkg>.BiocCheck" output folder into the working dir and
  # then flags any leftover one on the *next* run. Clear it before and after so
  # the checkBiocCheckOutputFolder false-positive never fires. The glob matches
  # whatever the package is called.
  clean_bc <- function() unlink(list.files(".", pattern = "\\.BiocCheck$",
                                           full.names = TRUE), recursive = TRUE)
  clean_bc()
  # no-check-version-num = TRUE: the deploy skill bumps the version by hand
  # (step 6), so let BiocCheck skip the version-number check.
  res <- BiocCheck::BiocCheck(".", `no-check-version-num` = TRUE,
                              `quit-with-status` = FALSE)
  # error/warning/note are named lists keyed by check id; getNum() is the
  # authoritative count. Each value is a list of offending locations.
  n   <- res$getNum()
  fmt <- function(cat) {
    x <- res[[cat]]
    vapply(names(x), function(nm)
      paste0(nm, "\n    ", paste(unlist(x[[nm]]), collapse = "\n    ")),
      character(1))
  }
  errs <- fmt("error"); warns <- fmt("warning"); notes <- fmt("note")
  log <- file.path(out, "bioccheck.log")
  writeLines(c("== ERRORS ==", errs, "", "== WARNINGS ==", warns, "",
               "== NOTES ==", notes), log)
  message(sprintf("bioccheck: %d error(s) | %d warning(s) | %d note(s)  (%s)",
                  n[["error"]], n[["warning"]], n[["note"]], log))

  # These two checks fail for environmental reasons, not code problems:
  #   checkSupportReg  - network call to the Bioc support site (flaky HTTP)
  #   checkBiocCheckOutputFolder - leftover output folder (we clean it anyway)
  env_only <- c("checkSupportReg", "checkBiocCheckOutputFolder")
  real_err <- setdiff(names(res$error), env_only)
  env_err  <- intersect(names(res$error), env_only)
  if (length(real_err)) { message("-- ERRORS --");   for (m in fmt("error")[real_err]) message("  * ", m) }
  if (length(env_err))  { message("-- ERRORS (environmental, non-blocking) --"); for (m in env_err) message("  * ", m) }
  if (n[["warning"]])   { message("-- WARNINGS --"); for (m in warns) message("  * ", m) }
  if (n[["note"]])      { message("-- NOTES --");    for (m in notes) message("  * ", m) }
  clean_bc()
  if (length(real_err) > 0L || n[["warning"]] > 0L) fail <<- fail + length(real_err) + n[["warning"]]
  invisible(length(real_err) == 0L && n[["warning"]] == 0L)
}

switch(phase,
  document  = run_document(),
  check     = run_check(),
  bioccheck = run_bioccheck(),
  all = {
    run_document()
    if (!run_check()) {
      message("\ncheck reported ERROR(s) -- fix these before bioccheck.")
    } else {
      run_bioccheck()
    }
  }
)

if (fail > 0L) {
  message("\nDEPLOY GATE: ", fail, " blocking problem(s). Fix and re-run.")
  quit(status = 1L)
}
message("\ndeploy.R phase '", phase, "' completed with no blocking problems.")
