# Convert the package's Rd files to Quarto pages under reference/r/ with altdoc's Rd converter
# (altdoc is installed in the git-ignored _rlib), plus an index grouped by topic.
# Examples run at build time and show their output, except pages that reference-examples.yml
# excludes (shown, not run, with the reason on the page).
.libPaths(c(normalizePath("_rlib"), .libPaths()))
suppressPackageStartupMessages(library(altdoc))
`%||%` <- function(a, b) if (is.null(a)) b else a

man_dir <- "../man"
out_dir <- "reference/r"
unlink(out_dir, recursive = TRUE)
dir.create(out_dir, recursive = TRUE)

# altdoc's converter looks for its settings next to the package: give it a stub.
stub <- tempfile("altdoc-stub")
dir.create(file.path(stub, "altdoc"), recursive = TRUE)
file.create(file.path(stub, "altdoc", "quarto_website.yml"))
file.copy("../DESCRIPTION", file.path(stub, "DESCRIPTION"))

rd_files <- list.files(man_dir, pattern = "\\.Rd$", full.names = TRUE)

# altdoc:::.rd2qmd is private API (altdoc is pinned in build_reference.sh). If it is ever
# missing, fall back to base R's Rd2HTML, embedded as raw HTML in a qmd page.
have_rd2qmd <- exists(".rd2qmd", envir = asNamespace("altdoc"), inherits = FALSE)
if (!have_rd2qmd) message("altdoc:::.rd2qmd not found; using tools::Rd2HTML fallback")
for (f in rd_files) {
  if (have_rd2qmd) {
    altdoc:::.rd2qmd(f, target_dir = out_dir, path = stub)
  } else {
    html <- tempfile(fileext = ".html")
    tools::Rd2HTML(f, out = html)
    h <- readLines(html)
    h <- h[seq(grep("<body", h)[1] + 1, tail(grep("</body>", h), 1) - 1)]
    writeLines(c("```{=html}", h, "```"),
               file.path(out_dir, sub("\\.Rd$", ".qmd", basename(f))))
  }
}

# Which examples run is decided by reference-examples.yml (see its header). Every Rd file must be
# listed there. Pages that run load the shared fitted object (built once per site build by
# build-shared-fit-r.R) in a hidden setup chunk when `shared: true`.
spec <- yaml::read_yaml("reference-examples.yml")$r
all_pages <- sub("\\.Rd$", "", basename(rd_files))
unlisted <- setdiff(all_pages, names(spec))
if (length(unlisted)) stop("not listed in reference-examples.yml (r): ", paste(unlisted, collapse = ", "), call. = FALSE)
stale <- setdiff(names(spec), all_pages)
if (length(stale)) stop("listed in reference-examples.yml but no such Rd file: ", paste(stale, collapse = ", "), call. = FALSE)

setup_chunk <- function(shared) {
  c("```{r}", "#| include: false",
    'if (dir.exists("../../_rlib")) .libPaths(c(normalizePath("../../_rlib"), .libPaths()))',
    "suppressPackageStartupMessages({library(SpaNorm); library(SpatialExperiment); library(ggplot2)})",
    if (shared) c(
      '# the shared fitted object, built once per site build; the examples call data(HumanDLPFC)',
      '# to load the raw data, so data() is masked for this page',
      'HumanDLPFC <- readRDS("../../_cache/HumanDLPFC-fit.rds")',
      "data <- function(...) invisible(NULL)"),
    "set.seed(1)", "```", "")
}

shared_note <- c("::: {.callout-note}", "The code below runs when the site is built, on a small pre-fitted subset of `HumanDLPFC` (301 genes x 4015 spots; `SpaNorm(sample.p = 0.05, df.tps = 2)`, `SpaNormSVG()` and `SpaNormPCA()` already applied) that is fitted once per build. The outputs are those of the code shown.", ":::", "")

for (f in rd_files) {
  nm <- sub("\\.Rd$", "", basename(f))
  qmd <- file.path(out_dir, paste0(nm, ".qmd"))
  x <- readLines(qmd)
  sp <- spec[[nm]]
  has_ex <- any(grepl("^### Examples", x))
  if (isTRUE(sp$run) && !has_ex) stop(nm, ": run: true but the Rd has no example", call. = FALSE)
  if (isTRUE(sp$run)) {
    if (!is.null(sp$code)) {   # replace the Rd example (a pointer such as example(SpaNorm))
      i <- grep("^### Examples", x)[1]
      x <- c(x[seq_len(i)], "```{r}", sub("\n$", "", sp$code), "```", "")
    }
    x <- gsub("eval\\s*=\\s*(TRUE|FALSE)", "eval=TRUE", x)
    x <- gsub("message\\s*=\\s*(TRUE|FALSE)", "message=TRUE", x)
    if (isTRUE(sp$shared)) x <- append(x, shared_note, after = grep("^### Examples", x)[1])
    hdr <- c("---", sprintf("title: \"%s\"", nm), "toc: false", "execute:", "  eval: true",
             "  error: false", "  fig-width: 7", "  fig-height: 5", "  fig-dpi: 144", "---", "",
             setup_chunk(isTRUE(sp$shared)))
  } else {
    if (!nzchar(sp$reason %||% "")) stop(nm, ": run: false needs a reason in reference-examples.yml", call. = FALSE)
    x <- gsub("eval\\s*=\\s*TRUE", "eval=FALSE", x)
    if (has_ex) x <- append(x, c("", "::: {.callout-note}", paste0("This example is shown but not run when the site is built: ", sp$reason, "."), ":::", ""), after = grep("^### Examples", x)[1])
    hdr <- c("---", sprintf("title: \"%s\"", nm), "toc: false", "execute:", "  eval: false", "---", "")
  }
  writeLines(c(hdr, x), qmd)
}

# no chunk may be left executable on a page that is not listed as run: true
for (f in list.files(out_dir, pattern = "\\.qmd$", full.names = TRUE)) {
  nm <- sub("\\.qmd$", "", basename(f))
  if (nm != "index" && !isTRUE(spec[[nm]]$run) && any(grepl("eval\\s*[=:]\\s*(TRUE|true)", readLines(f)))) {
    stop("executable chunk left in ", f, call. = FALSE)
  }
}

# index grouped by topic; Rd title taken from the Rd files
groups <- list(
  "Normalisation" = c("SpaNorm", "SpaNormFit", "HumanDLPFC"),
  "SVGs and PCA" = c("SpaNormSVG", "topSVGs", "getSVGResults", "SpaNormPCA"),
  "Polishing" = c("polishSpaNorm", "isPolished", "polishNB", "nbProfilePsi",
                  "nbNewtonSolver", "nbMuFloor", "nbBlockDesign"),
  "Model fitting and helpers" = c("fitNB", "calculateMu", "filterGenes", "fastSizeFactors",
                                   "tpsBasis", "qlDispersion", "qlMomentTable",
                                   "nbUnitDeviance", "nbDevianceMoments"),
  "Plotting" = c("plotSpatial", "plotCovariate")
)
groups[["GPU and backend"]] <- setdiff(all_pages, unlist(groups))
stopifnot(all(unlist(groups) %in% all_pages))

rd_title <- function(nm) {
  rd <- tools::parse_Rd(file.path(man_dir, paste0(nm, ".Rd")))
  tag <- vapply(rd, function(e) attr(e, "Rd_tag"), "")
  trimws(gsub("\\s+", " ", paste(unlist(rd[tag == "\\title"]), collapse = "")))
}

idx <- c("---", "title: \"R reference\"", "toc: false", "---", "")
for (g in names(groups)) {
  idx <- c(idx, paste("##", g), "", "| | |", "| --- | --- |")
  for (nm in sort(groups[[g]])) {
    idx <- c(idx, sprintf("| [`%s`](%s.qmd) | %s |", nm, nm, rd_title(nm)))
  }
  idx <- c(idx, "")
}
writeLines(idx, file.path(out_dir, "index.qmd"))
message("wrote ", length(rd_files), " R reference pages")
