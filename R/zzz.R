# ggplot2::aes() references these as data-frame columns via non-standard
# evaluation, but R CMD check's static analysis flags bare NSE symbols as
# undefined globals (see plotSpatial.R). Declaring them here is a no-op at
# runtime; it only silences the check NOTE.
utils::globalVariables(c("x", "y", "colour", "r"))
