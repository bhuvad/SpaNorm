# Convert the package's Rd files to Quarto pages under reference/r/ with altdoc's Rd converter
# (altdoc is installed in the git-ignored _rlib), plus an index grouped by topic.
# Examples are shown but not run on the site.
.libPaths(c(normalizePath("_rlib"), .libPaths()))
suppressPackageStartupMessages(library(altdoc))

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

# front matter (page title = Rd name), examples shown but not executed
for (f in rd_files) {
  nm <- sub("\\.Rd$", "", basename(f))
  qmd <- file.path(out_dir, paste0(nm, ".qmd"))
  x <- readLines(qmd)
  x <- gsub("eval\\s*=\\s*TRUE", "eval=FALSE", x)
  writeLines(c("---", sprintf("title: \"%s\"", nm), "toc: false",
               "execute:", "  eval: false", "---", "", x), qmd)
}

# examples are never run on the site: no chunk may be left executable
for (f in list.files(out_dir, pattern = "\\.qmd$", full.names = TRUE)) {
  if (any(grepl("eval\\s*[=:]\\s*(TRUE|true)", readLines(f)))) {
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
all_pages <- sub("\\.Rd$", "", basename(rd_files))
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
