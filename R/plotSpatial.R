#' Plot spatial transcriptomic annotations per spot
#'
#' @param spe a SpatialExperiment object.
#' @param what a character, specifying what aspect should be plot, "annotation",
#'   "expression", or reduced dimension ("reduceddim").
#' @param ... additional aesthetic mappings or fixed parameters (e.g., shape = ".").
#' @param assay a character or numeric, specifying the assay to plot (default is the first assay).
#' @param dimred a character or numeric, specifying the reduced dimension to plot (default is the first reduced dimension).
#' @param img a logical, indicating whether the tissue image (if present) should be plot (default = FALSE).
#' @param crop a logical, indicating whether the image should be cropped to the spatial coordinates (default = FALSE).
#' @param imgAlpha a numeric, specifying the alpha value for the image (default = 1).
#' @param rl a numeric, specifying the relative size of the text (default = 1).
#' @param circles a logical, indicating whether the spots should be plotted as circles (default = FALSE). This can be slower for large datasets.
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' library(SpatialExperiment)
#' library(ggplot2)
#' 
#' # load data
#' data(HumanDLPFC)
#' 
#' # plot spatial region annotations
#' p1 <- plotSpatial(HumanDLPFC, colour = AnnotatedCluster)
#' p1
#' 
#' # change colour scale
#' p1 + scale_colour_brewer(palette = "Paired")
#' 
#' # plot spatial expression
#' plotSpatial(HumanDLPFC, what = "expression", colour = ENSG00000075624) +
#'  scale_colour_viridis_c(option = "F")
#' 
#' # plot logcounts
#' logcounts(HumanDLPFC) <- log2(counts(HumanDLPFC) + 1)
#' plotSpatial(HumanDLPFC, what = "expression", colour = ENSG00000075624, assay = "logcounts") +
#'  scale_colour_viridis_c(option = "F")
#' 
#' # change point shape
#' plotSpatial(HumanDLPFC, what = "expression", colour = ENSG00000075624, assay = "logcounts", shape = 18) +
#'  scale_colour_viridis_c(option = "F")
#' 
plotSpatial <-
  function(spe,
           what = c("annotation", "expression", "reduceddim"),
           ...,
           assay = SummarizedExperiment::assayNames(spe),
           dimred = SingleCellExperiment::reducedDimNames(spe),
           img = FALSE,
           crop = FALSE,
           imgAlpha = 1,
           rl = 1,
           circles = FALSE) {
    stopifnot(is(spe, "SpatialExperiment"))
    stopifnot(imgAlpha >= 0 & imgAlpha <= 1)
    stopifnot(rl > 0)
    what = match.arg(what)
    assay = match.arg(assay)
    dimred = match.arg(dimred)

    #----extract aes----
    aesmap = rlang::enquos(...)
    #compute plot
    aesmap = aesmap[!names(aesmap) %in% c('x0', 'y0', 'x', 'y', 'sf')]#remove x,y mappings if present
    # split aes params into those that are not aes i.e. static parametrisation
    if (length(aesmap) > 0) {
      is_aes = sapply(aesmap, rlang::quo_is_symbolic)
      defaultmap = lapply(aesmap[!is_aes], rlang::eval_tidy)
      aesmap = aesmap[is_aes]
      pltcols = as.character(sapply(aesmap, rlang::quo_get_expr))
    } else {
      defaultmap = list()
    }

    #----prepare spot data----
    plotdf = switch (what,
      annotation = extractAnnotation(spe, pltcols),
      expression = extractExpression(spe, assay, pltcols),
      reduceddim = extractReducedDim(spe, dimred, pltcols)
    )

    #----image data----
    imgdf = NULL
    sf = rep(1, ncol(spe))
    if (img) {
      imgdf = extractImage(spe)
      sf = SpatialExperiment::scaleFactors(spe)
      names(sf) = SpatialExperiment::imgData(spe)$sample_id
      sf = sf[spe$sample_id]
    }

    #----add spatial coordinates and scales----
    spatdf = SpatialExperiment::spatialCoords(spe)
    colnames(spatdf) = c('x', 'y')
    plotdf = cbind(plotdf, spatdf)
    plotdf = as.data.frame(plotdf)
    plotdf$sf = sf

    # crop
    if (crop & img) {
      #compute lims
      xlim = range((plotdf$x * plotdf$sf))
      ylim = range((plotdf$y * plotdf$sf))

      #filter image data - remove out-of-bounds pixels
      imgdf = imgdf[imgdf$x >= xlim[1] &
                  imgdf$x <= xlim[2] &
                  imgdf$y >= ylim[1] &
                  imgdf$y <= ylim[2], ]
    }

    #----plot----
    #initialise plot
    p1 = ggplot2::ggplot()

    #image
    if (img) {
      requirePkg('ggnewscale')
      p1 = p1 +
        ggplot2::geom_raster(ggplot2::aes(x, y, fill = colour),
                             alpha = imgAlpha,
                             data = imgdf) +
        ggplot2::scale_fill_identity() +
        ggnewscale::new_scale_fill()
    }

    #plot spots/circles
    if (circles) {
      requirePkg('ggforce')
      #add radius if not present
      if (!'r' %in% names(aesmap)) {
        aesmap = c(aesmap, r = rlang::quo_set_env(rlang::quo(r), rlang::quo_get_env(aesmap[[1]])))
        plotdf$r = ifelse(is.null(defaultmap$r), 100, defaultmap$r)
        defaultmap = defaultmap[setdiff(names(defaultmap), 'r')]
      }

      #add scale factor
      scquo = paste0(rlang::quo_text(aesmap$r), ' * sf')
      aesmap$r = rlang::set_expr(aesmap$r, rlang::parse_expr(scquo))

      plt_params = list(
        mapping = ggplot2::aes(, , x0 = x * sf, y0 = y * sf,!!!aesmap),
        data = plotdf
      ) |> c(defaultmap)
      p1 = p1 + do.call(ggforce::geom_circle, plt_params)
    } else {
      plt_params = list(
        mapping = ggplot2::aes(x = x * sf, y = y * sf, !!!aesmap),
        data = plotdf
      ) |> c(defaultmap)
      p1 = p1 + do.call(ggplot2::geom_point, plt_params)
    }

    #----theme----
    p1 = p1 +
      custom_theme(rl) +
      ggplot2::labs(x = SpatialExperiment::spatialCoordsNames(spe)[1],
                    y = SpatialExperiment::spatialCoordsNames(spe)[2])

    return(p1)
  }

requirePkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("'%s' needed for this function to work.", pkg), call. = FALSE)
  }
}

custom_theme <- function(rl = 1.1) {
  stopifnot(rl > 0)
  rl = ggplot2::rel(rl)
  ggplot2::theme_minimal() + ggplot2::theme(
    panel.border = ggplot2::element_rect(colour = "black", fill = NA),
    panel.grid = ggplot2::element_blank(),
    axis.title = ggplot2::element_blank(),
    axis.text = ggplot2::element_text(size = rl),
    plot.title = ggplot2::element_text(size = rl * 1.2),
    strip.background = ggplot2::element_rect(fill = NA, colour = "black"),
    strip.text = ggplot2::element_text(size = rl),
    legend.text = ggplot2::element_text(size = rl),
    legend.title = ggplot2::element_text(size = rl, face = "italic"),
    legend.position = "bottom"
  )
}

extractAnnotation <- function(spe, pltcols) {
  df = SummarizedExperiment::colData(spe)
  pltcols = intersect(pltcols, colnames(df))
  df[, pltcols, drop = FALSE]
}

extractExpression <- function(spe, assay, pltcols) {
  emat = SummarizedExperiment::assay(spe, assay)
  pltcols = intersect(pltcols, rownames(emat))
  emat[pltcols, , drop = FALSE] |>
    as.matrix() |>
    t() |>
    as.data.frame()
}

extractReducedDim <- function(spe, dimred, pltcols) {
  dmat = SingleCellExperiment::reducedDim(spe, type = dimred)
  pltcols = intersect(pltcols, colnames(dmat))
  dmat[, pltcols, drop = FALSE] |>
    as.data.frame()
}

extractImage <- function(spe) {
  if (nrow(SpatialExperiment::imgData(spe)) == 0) {
    stop("no image data found in 'imgData'")
  }

  if (any(duplicated(SpatialExperiment::imgData(spe)$sample_id))) {
    stop("multiple records found for some samples in 'imgData'")
  }

  # get img data
  imgdf = SpatialExperiment::imgData(spe)$data
  names(imgdf) = SpatialExperiment::imgData(spe)$sample_id
  imgdf = mapply(function(x, sample_id) {
    x = as.matrix(SpatialExperiment::imgRaster(x))
    nc = ncol(x)
    nr = nrow(x)
    data.frame(
      x = rep(seq(nc), each = nr),
      y = rep(seq(nr), nc),
      colour = as.character(x),
      sample_id = sample_id
    )
  }, imgdf, names(imgdf), SIMPLIFY = FALSE) |>
    do.call(what = rbind)

  return(imgdf)
}