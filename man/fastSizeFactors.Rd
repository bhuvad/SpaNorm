% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fastSizeFactors.R
\name{fastSizeFactors}
\alias{fastSizeFactors}
\alias{fastSizeFactors,SpatialExperiment-method}
\title{Filter genes based on expression}
\usage{
fastSizeFactors(spe)

\S4method{fastSizeFactors}{SpatialExperiment}(spe)
}
\arguments{
\item{spe}{a SpatialExperiment, Seurat, or SpatialFeatureExperiment object containing count data.}
}
\value{
a SpatialExperiment, Seurat, or SpatialFeatureExperiment, containing size factors in the 'sizeFactor' column of the column annotation.
}
\description{
This function computes the size factors using a fast but inaccurate approach. Size factors are computed using the direct estimate of library sizes (sum of all counts). Though fast, this approach does not cater for compositional biases in the data and therefore is less accurate than scran-based estimates.
}
\examples{
data(HumanDLPFC)
HumanDLPFC <- fastSizeFactors(HumanDLPFC)
head(HumanDLPFC$sizeFactor)

}
