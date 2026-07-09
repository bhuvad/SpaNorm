if (interactive() && Sys.getenv("RSTUDIO") == "") {
  Sys.setenv(TERM_PROGRAM = "vscode")
  source(file.path(Sys.getenv(
    if (.Platform$OS.type == "windows") "USERPROFILE" else "HOME"
  ), ".vscode-R", "init.R"))

  # R >= 4.6 workaround: in some environments the VS Code hook is assigned
  # but not executed, which leaves .vsc.attach unavailable.
  if (
    !exists(".vsc.attach", mode = "function") &&
      exists(".First.sys", envir = globalenv(), inherits = FALSE)
  ) {
    try(get(".First.sys", envir = globalenv())(), silent = TRUE)
  }
}

options(vsc.rstudioapi = TRUE)

options(languageserver.formatting_style = function(options) {
  style <- styler::tidyverse_style(indent_by = options$tabSize)
  style$token$force_assignment_op <- NULL
  style
})
