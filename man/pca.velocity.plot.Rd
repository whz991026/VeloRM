% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pca.velocity.plot.R
\name{pca.velocity.plot}
\alias{pca.velocity.plot}
\title{PCA-based visualization of the velocities}
\usage{
pca.velocity.plot(
  current,
  projected,
  deltaE,
  nPcs = 4,
  cell.colors = NULL,
  scale = "log",
  plot.cols = min(3, nPcs - 1),
  norm.nPcs = NA,
  pc.multipliers = NULL,
  show.grid.flow = FALSE,
  grid.n = 20,
  grid.sd = NULL,
  arrow.scale = 1,
  min.grid.cell.mass = 1,
  min.arrow.size = NULL,
  pcount = 1,
  arrow.lwd = 1,
  size.norm = FALSE,
  return.details = FALSE,
  plot.grid.points = FALSE,
  fixed.arrow.length = FALSE,
  max.grid.arrow.length = NULL,
  n.cores = 1,
  point.size = 3,
  cell.border.alpha = 0.5,
  arrow_size = 0.3,
  ...
)
}
\arguments{
\item{current}{current matrix}

\item{projected}{projected matrix}

\item{deltaE}{deltaE matrix}

\item{nPcs}{number of successive PCs to visualize}

\item{cell.colors}{a named vector of cell colors for visualization}

\item{scale}{scale to use for expression state transform (default: 'log', other possible values are 'sqrt','linear')}

\item{plot.cols}{number of columns into which to arrange the plots}

\item{norm.nPcs}{optional total number of PCs to use for velocity magnitude normalization}

\item{pc.multipliers}{an optional vector multipliers for the cell PC scores (useful for reorienting the PCs)}

\item{show.grid.flow}{whether a grid flow should be shown}

\item{grid.n}{number of grid points (on each axis)}

\item{grid.sd}{standard deviation of the grid}

\item{arrow.scale}{scale multiplier for the velocity estimates}

\item{min.grid.cell.mass}{minimum cellular mass}

\item{min.arrow.size}{minimum size of an arrow to show}

\item{pcount}{pseudocount}

\item{arrow.lwd}{thickness of arrows to plot}

\item{size.norm}{whether to rescale current and projected states by cell size (default=FALSE)}

\item{return.details}{whether to return detailed output}

\item{plot.grid.points}{whether to show dots at every grid point}

\item{fixed.arrow.length}{whether to use fixed-size arrow}

\item{max.grid.arrow.length}{limit to the size of the arrows that could be shown (when fixed.arrow.length=FALSE)}

\item{n.cores}{number of cores to use in the calculations}

\item{point.size}{size of the point}

\item{cell.border.alpha}{transparency for the cell border}

\item{arrow_size}{size of arrow}

\item{...}{extra parameters are passed to plot() function}
}
\value{
If return.details=F, returns invisible list containing PCA info (epc) and projection of velocities onto the PCs (delta.pcs). If return.details=T, returns an extended list that can be passed into p1 app for velocity visualization.
}
\description{
PCA-based visualization of the velocities
}
