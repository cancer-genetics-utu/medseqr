\name{ggbiplot.n}
\alias{ggbiplot.n}
\title{Make a biplot of PCA output data using ggplot2.}
\usage{
  ggbiplot.n(pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE,
             obs.scale = 1 - scale, var.scale = scale, groups = NULL,
             grp.nam = NULL, ellipse = FALSE, ellipse.prob = 0.68, labels = NULL,
             labels.size = 3, alpha = 1, var.axes = TRUE, circle = FALSE,
             circle.prob = 0.69, varname.size = 3, varname.adjust = 1.5,
             varname.abbrev = FALSE, point.size = 1, axes.title.size = 1,
             legend.text.size = 1, legend.title.size = 1, ellipse.lwd = 1,
             main = NULL, fix.aspect = TRUE, tweak = TRUE, ...)
}
\arguments{
  \item{pcobj}{an object returned by prcomp() or princomp()}

  \item{choices}{\code{length-two numeric}. which PCs to plot; default: 1:2}

  \item{scale}{\code{length-one numeric}. covariance biplot (scale = 1) [default], form biplot (scale = 0). When scale = 1, the inner product between the variables approximates the covariance and the distance between the points approximates the Mahalanobis distance.}

  \item{pc.biplot}{\code{logical}. for compatibility with biplot.princomp()}

  \item{obs.scale}{\code{length-one numeric}. scale factor to apply to observations; default: 1-scale}

  \item{var.scale}{\code{length-one numeric}. scale factor to apply to variables; default: scale}

  \item{groups}{\code{factor}. optional factor variable indicating the groups that the observations belong to. If provided the points will be colored according to groups}

  \item{grp.nam}{\code{character}. optional name of the grouping variable to be used as legend title}

  \item{ellipse}{\code{logical}. draw a normal data ellipse for each group?}

  \item{ellipse.prob}{\code{double}. size of the ellipse in Normal probability}

  \item{labels}{\code{character}. optional vector of labels for the observations}

  \item{labels.size}{\code{double}. size of the text used for the labels}

  \item{alpha}{\code{double}. alpha transparency value for the points (0 = transparent, 1 = opaque)}

  \item{var.axes}{\code{logical}. draw arrows for the variables?}

  \item{circle}{\code{logical}. draw a correlation circle? (only applies when prcomp was called with scale = TRUE and when var.scale = 1)}

  \item{circle.prob}{\code{double}. size of the circle in Normal probability}

  \item{varname.size}{\code{double}. size of the text for variable names}

  \item{varname.adjust}{\code{double}. adjustment factor the placement of the variable names, >= 1 means farther from the arrow}

  \item{varname.abbrev}{\code{logical}. whether or not to abbreviate the variable names}

  \item{point.size}{\code{double}. expansion factor for point size; uses \code{rel()} internally}

  \item{axes.title.size}{\code{double}. expansion factor for axes title sizes; uses \code{rel()} internally}

  \item{legend.text.size}{\code{double}. expansion factor for legend text size; uses \code{rel()} internally}

  \item{legend.title.size}{\code{double}. expansion factor for legend title size; uses \code{rel()} internally}

  \item{ellipse.lwd}{\code{double}. expansion factor for ellipse line width; uses \code{rel()} internally}

  \item{main}{\code{character}. Plot title. If \code{NULL} (default) no title will be added.}

  \item{fix.aspect}{\code{logical}. Should the aspect ratio of the x- and y-axes be kept constant for different plot sizes?}

  \item{tweak}{\code{logical}. Should the plot theme be tweaked? Will apply values set in axes.title.size, legend.text.size and legend.title.size. Defaults to TRUE.}

  \item{...}{currently not in use}
}
\value{
  The final plot object returned by ggplot.
}
\description{
  Make a biplot of PCA output data using ggplot2.
}
