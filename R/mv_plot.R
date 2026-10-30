#' @include multivariate.R mvgaussian_distrib.R ggplot.R
NULL


#' @title Panels of a Multivariate Density
#' @name plot.multivariate_distrib
#'
#' @description
#' `autoplot()` returns a graph of a multivariate distribution as a matrix of
#' panels, and `plot()` prints that graph: the marginal density of each
#' coordinate on the diagonal, contours of each bivariate marginal below it,
#' and the implied correlation printed above it. The density panels are
#' built from marginals, so the graph exists where [mv_marginal()] returns a
#' marginal of several coordinates, which holds for the multivariate gaussian
#' and Student \eqn{t}. For the Dirichlet and the multinomial, whose marginals
#' are returned one coordinate at a time, the graph signals the error of
#' [mv_marginal()].
#'
#' @details
#' # At most three coordinates
#'
#' A panel matrix of \eqn{k} coordinates has \eqn{k^2} panels, and past three
#' the panels are too small to read. A larger distribution signals an error that
#' names the number of coordinates and suggests `which`.
#'
#' # The contour levels
#'
#' The contours are levels of the density and not of the probability content of
#' the enclosed region, which would require the mass of each level set of the
#' density. The levels are the values that `pretty()` selects, about six,
#' over the range of the density on the grid.
#'
#' # One setting only
#'
#' A univariate `plot()` method reads a `theta` component of length \eqn{k} as
#' \eqn{k} settings drawn over one another. Here the picture is already a matrix
#' of panels and has no axis on which to overlay them, so a component longer
#' than one signals an error that names it.
#'
#' # The plotting range
#'
#' Each coordinate is drawn over three standard deviations either side of the
#' location, and the spread is taken from the matrix that the parametrization
#' carries (see [mv_sigma()]). A multivariate Student \eqn{t} with \eqn{\nu \le
#' 2} has an infinite variance, so [variance()] cannot be the source. For an
#' elliptical family the two matrices agree up to a factor, which is enough for
#' a plotting range.
#'
#' The panels are drawn with \pkg{ggplot2} and joined with \pkg{patchwork};
#' both packages are suggested, and the methods signal an error that names
#' them when they are not installed.
#'
#' @param x,object An object inheriting from [multivariate_distrib()].
#' @param theta A named list or vector of parameters in which each component is
#'   a single number. When missing, parameters are drawn by
#'   [generate_random_theta()] and reported in a message. A component of
#'   length greater than one is an error.
#' @param which An integer vector of the coordinates to show, at most three.
#'   Defaults to `NULL`, which means all of them and is itself an error above
#'   three coordinates.
#' @param n_grid Points per axis for the contours. The marginal density curves
#'   on the diagonal use `4 * n_grid` points. Defaults to `80`. The bivariate
#'   panels evaluate the density on an `n_grid * n_grid` grid, so the cost is
#'   quadratic in it.
#' @param legend Logical of length 1: draw the legend below the panel matrix.
#'   Defaults to `TRUE`.
#' @param ... Unused; an argument passed here signals an error.
#'
#' @return `autoplot()` returns a patchwork object; `plot()` prints the same
#'   object and returns it invisibly.
#'
#' @seealso [mv_marginal()], which supplies the density panels,
#'   [plot.distrib_fit()], which draws the same matrix with the data on it,
#'   and [plot.continuous_distrib()] for the one-dimensional case, which does
#'   overlay several settings.
#'
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("patchwork", quietly = TRUE)
#' d <- mvgaussian1_distrib(3)
#' theta <- as.list(stats::setNames(
#'   c(0, 1, -1, 0, 0, 0, 0.6, -0.3, 0.2), d@params
#' ))
#' plot(d, theta)
#'
#' # Two coordinates of a heavy-tailed family are drawn, and its contours are
#' # wider than a gaussian's at the same matrix.
#' t2 <- mvstudent_t1_distrib(3)
#' th2 <- as.list(stats::setNames(c(unlist(theta), 3), t2@params))
#' plot(t2, th2, which = c(1, 2))
#'
#' # Four coordinates would give sixteen panels, so the call signals an error.
#' d4 <- mvgaussian1_distrib(4)
#' th4 <- as.list(stats::setNames(rep(0, d4@n_params), d4@params))
#' try(plot(d4, th4))
#'
#' @keywords internal
#' @export
plot.multivariate_distrib <- function(x, theta, which = NULL,
                                                   n_grid = 80, legend = TRUE,
                                                   ...) {
  p <- autoplot_multivariate_distrib(x, theta, which = which, n_grid = n_grid,
                                     legend = legend, ...)
  print(p)
  invisible(p)
}
S7::method(plot, multivariate_distrib) <- plot.multivariate_distrib

#' @rdname plot.multivariate_distrib
autoplot_multivariate_distrib <- function(object, theta, which = NULL,
                                          n_grid = 80, legend = TRUE, ...) {
  check_plot_args(legend, ...)
  x <- object
  if (missing(theta)) {
    theta <- generate_random_theta(x)
    message(sprintf(
      "Argument 'theta' is missing. The plot uses random parameters: %s",
      paste(names(theta), round(unlist(theta), 3), sep = " = ", collapse = ", ")
    ))
  }
  theta <- align_theta(x, theta)
  # The univariate methods read a vector component as several settings drawn
  # over one another. Here the picture is already a matrix of panels, with no
  # axis left to overlay them on, so the request is rejected rather than
  # answered with the first setting.
  long <- names(theta)[lengths(theta) > 1L]
  if (length(long)) {
    stop(sprintf(paste0("A multivariate distribution is drawn at one setting: ",
                        "%s carr%s more than one value. Plot the settings ",
                        "separately."),
                 paste(sprintf("'%s'", long), collapse = ", "),
                 if (length(long) == 1L) "ies" else "y"), call. = FALSE)
  }
  mv_pairs_panels(x, theta, which, n_grid, data = NULL, legend = legend)
}

S7::method(autoplot, multivariate_distrib) <- autoplot_multivariate_distrib


#' @title Draw the Panel Matrix of a Multivariate Density
#'
#' @description
#' This function is the shared engine of [plot.multivariate_distrib()] and the
#' multivariate branch of [plot.distrib_fit()]. It lays out one panel per
#' ordered pair of the chosen coordinates: the marginal density on the
#' diagonal, contours of the bivariate marginal below it, and the implied
#' correlation printed above it. The density panels come from [mv_marginal()],
#' so the function does not integrate anything; the printed correlation comes
#' from [mv_sigma()]. One legend below the matrix, built by
#' [mv_pairs_legend()], names the elements of every panel.
#'
#' @details
#' When `data` is supplied the diagonal also carries a kernel density estimate
#' of the observed coordinate and the panels below the diagonal carry the
#' observations, so the fitted shape and the sample are read against each other
#' in one frame. The panels above the diagonal show only the correlation. The
#' kernel estimate (`stats::density()`) does not depend on the fitted family and
#' serves as a model-free comparison.
#'
#' Each coordinate is drawn over three standard deviations either side of the
#' location, widened when `data` is supplied to include the observations. The
#' spread is taken from [mv_sigma()], so a Student \eqn{t} with \eqn{\nu \le 2}
#' still has a finite range, where [variance()] returns an infinite matrix.
#'
#' @param d A [multivariate_distrib()] object.
#' @param theta A named list of parameters, already aligned by the caller.
#' @param which An integer vector of distinct coordinates, or `NULL` for all
#'   of them. An empty, non-integer, out-of-range or duplicated `which`, or
#'   more than three coordinates, is an error.
#' @param n_grid Points per axis, a single positive whole number.
#' @param data An \eqn{n \times p} matrix of observations with one column
#'   per coordinate of the distribution, whichever coordinates `which`
#'   selects, or `NULL` for a distribution with no data behind it.
#' @param which_arg The name of the caller's argument that carries `which`,
#'   used in the error message for more than three coordinates. Defaults to
#'   `"which"`; [plot.distrib_fit()] passes `"mv_which"`.
#' @param legend Whether the legend is drawn below the matrix.
#'
#' @return A patchwork object.
#'
#' @seealso [plot.multivariate_distrib()] and [plot.distrib_fit()], its two
#'   callers, and [mv_marginal()] for the panels.
#'
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("patchwork", quietly = TRUE)
#' d <- mvgaussian1_distrib(2)
#' theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0,
#'               sigma_log_L2 = 0, sigma_L2.1 = 0.6)
#'
#' # Without data, the graph shows the density alone.
#' distributions7:::mv_pairs_panels(d, theta, NULL, 40, data = NULL)
#'
#' # With data: the same frame carries a kernel estimate and the observations.
#' set.seed(1)
#' distributions7:::mv_pairs_panels(d, theta, NULL, 40,
#'                                  data = distrib_rng(d, 200, theta))
#'
#' @keywords internal
mv_pairs_panels <- function(d, theta, which, n_grid, data,
                            which_arg = "which", legend = TRUE) {
  p <- d@n_dim
  if (is.null(which)) which <- seq_len(p)
  if (!is.numeric(which) || !length(which) || !all(is.finite(which)) ||
      any(which != round(which)) || any(which < 1) || any(which > p) ||
      anyDuplicated(which)) {
    stop(sprintf("'%s' must be distinct coordinates in 1:%d.", which_arg, p),
         call. = FALSE)
  }
  which <- as.integer(which)
  k <- length(which)

  if (k > 3L) {
    stop(sprintf(paste0(
      "A panel matrix of %d coordinates has %d panels and is not readable.\n",
      "  Choose at most three with '%s', for instance %s = c(1, 2, 3)."
    ), k, k * k, which_arg, which_arg), call. = FALSE)
  }
  need_ggplot(patchwork = TRUE)

  # The range each coordinate is drawn over, from its own marginal: three
  # standard deviations either side of the mean cover the shape without
  # flattening it into the axis.
  # The spread comes from the matrix that the matrix parameter carries rather than from
  # variance(): a Student t with nu <= 2 has no variance, and that is precisely
  # the shape worth drawing. The two agree up to a factor for an elliptical
  # family, which is all a plotting range needs.
  center <- as.numeric(mv_location(d, theta))
  spread <- sqrt(diag(mv_sigma(d, theta)))
  rng <- lapply(which, function(j) {
    lo <- center[j] - 3 * spread[j]
    hi <- center[j] + 3 * spread[j]
    if (!is.null(data)) {
      lo <- min(lo, min(data[, j]))
      hi <- max(hi, max(data[, j]))
    }
    c(lo, hi)
  })
  names(rng) <- as.character(which)

  lab <- paste0("v", which)
  key <- mv_pairs_keys(data)
  .data <- ggplot2::.data
  # every panel carries the same colour scale, so that its legend can be
  # removed and the matrix shows the one legend below it
  scale <- ggplot2::scale_colour_manual(values = key, guide = "none")
  panels <- vector("list", k * k)

  for (a in seq_len(k)) {
    for (b in seq_len(k)) {
      ja <- which[a]
      jb <- which[b]
      if (a == b) {
        m <- mv_marginal(d, theta, ja)
        xs <- seq(rng[[as.character(ja)]][1], rng[[as.character(ja)]][2],
          length.out = n_grid * 4L
        )
        dens <- distrib_pdf(m$distrib, matrix(xs, ncol = 1L), m$theta)
        g <- ggplot2::ggplot() +
          ggplot2::geom_line(data = data.frame(x = xs, y = dens),
                             ggplot2::aes(x = .data$x, y = .data$y,
                                          colour = names(key)[1L]),
                             linewidth = 0.8)
        if (!is.null(data)) {
          kd <- stats::density(data[, ja])
          g <- g + ggplot2::geom_line(data = data.frame(x = kd$x, y = kd$y),
                                      ggplot2::aes(x = .data$x, y = .data$y,
                                                   colour = names(key)[2L]),
                                      linewidth = 0.8)
        }
        g <- g + ggplot2::labs(x = lab[a], y = "density")
      } else if (a > b) {
        m <- mv_marginal(d, theta, c(jb, ja))
        gx <- seq(rng[[as.character(jb)]][1], rng[[as.character(jb)]][2],
          length.out = n_grid
        )
        gy <- seq(rng[[as.character(ja)]][1], rng[[as.character(ja)]][2],
          length.out = n_grid
        )
        grid <- expand.grid(x = gx, y = gy)
        grid$z <- distrib_pdf(m$distrib, as.matrix(grid[, c("x", "y")]),
                              m$theta)
        g <- ggplot2::ggplot()
        if (!is.null(data)) {
          g <- g + ggplot2::geom_point(
            data = data.frame(x = data[, jb], y = data[, ja]),
            ggplot2::aes(x = .data$x, y = .data$y, colour = names(key)[3L]),
            size = 0.7, alpha = 0.35)
        }
        g <- g + ggplot2::geom_contour(
            data = grid,
            ggplot2::aes(x = .data$x, y = .data$y, z = .data$z,
                         colour = names(key)[1L]),
            breaks = pretty(range(grid$z), 6L), linewidth = 0.6) +
          ggplot2::coord_cartesian(xlim = range(gx), ylim = range(gy),
                                   expand = FALSE) +
          ggplot2::labs(x = lab[b], y = lab[a])
      } else {
        # The upper triangle repeats the lower one with the axes swapped, so it
        # is left empty and used for the correlation instead, which is the one
        # number a reader wants from a pair and cannot read off a contour.
        s <- mv_sigma(d, theta)
        r <- s[ja, jb] / sqrt(s[ja, ja] * s[jb, jb])
        g <- ggplot2::ggplot() +
          ggplot2::annotate("text", x = 0.5, y = 0.5, size = 4,
            colour = "grey20",
            label = sprintf("%s of %s and %s\n%.3f",
              if (is.null(data)) "correlation" else "fitted correlation",
              lab[a], lab[b], r)) +
          ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
          ggplot2::labs(x = NULL, y = NULL)
      }
      g <- g + scale + ggplot2::theme_bw()
      if (a < b) {
        g <- g + ggplot2::theme(axis.text = ggplot2::element_blank(),
                                axis.ticks = ggplot2::element_blank(),
                                panel.grid = ggplot2::element_blank())
      }
      panels[[(a - 1L) * k + b]] <- g
    }
  }
  out <- patchwork::wrap_plots(panels, ncol = k, byrow = TRUE)
  if (!legend) return(out)
  leg <- mv_pairs_legend(data)
  patchwork::wrap_plots(out, patchwork::wrap_elements(full = leg), ncol = 1L,
                        heights = grid::unit.c(grid::unit(1, "null"),
                                               sum(leg$heights)))
}


#' The Elements of a Panel Matrix and Their Colours
#'
#' @description
#' Returns the names of the elements that [mv_pairs_panels()] draws, as the
#' names of a vector of colours: the density (its marginal curves and its
#' contours) and, when there are data, the kernel density estimate of each
#' marginal and the observations.
#'
#' @param data The data matrix, or `NULL`.
#'
#' @return A named character vector of colours, of length 1 without data and
#'   3 with data.
#'
#' @keywords internal
mv_pairs_keys <- function(data) {
  if (is.null(data)) {
    return(c("density: marginal curves and contours" = "#1f5fa8"))
  }
  c("fitted density: marginal curves and contours" = "#1f5fa8",
    "kernel density estimate of a marginal of the observations" = "grey30",
    "observations" = "grey30")
}


#' The Legend of a Panel Matrix
#'
#' @description
#' Builds the legend that [mv_pairs_panels()] places below the panel matrix,
#' as a grob: a graph that draws each element of the matrix once is built
#' with its legend below the panel, and the legend is taken from it. The
#' entries are stacked one per row.
#'
#' @param data The data matrix, or `NULL`.
#'
#' @return A gtable holding the legend.
#'
#' @keywords internal
mv_pairs_legend <- function(data) {
  key <- mv_pairs_keys(data)
  .data <- ggplot2::.data
  src <- data.frame(x = c(0, 1), y = c(0, 1))
  g <- ggplot2::ggplot(src, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_line(ggplot2::aes(colour = names(key)[1L]), linewidth = 0.8)
  if (!is.null(data)) {
    g <- g +
      ggplot2::geom_line(ggplot2::aes(colour = names(key)[2L]),
                         linewidth = 0.8) +
      ggplot2::geom_point(ggplot2::aes(colour = names(key)[3L]), alpha = 0.6)
  }
  # one entry per row, as the entries are long and the panels share the width
  g <- g + ggplot2::scale_colour_manual(
      name = NULL, values = key, breaks = names(key),
      guide = ggplot2::guide_legend(ncol = 1L)) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom",
                   legend.key.width = ggplot2::unit(2, "lines"),
                   legend.key.spacing.y = ggplot2::unit(0, "pt"))
  gt <- ggplot2::ggplotGrob(g)
  gt$grobs[[which(gt$layout$name == "guide-box-bottom")]]
}


#' @title Marginal of a Multivariate Gaussian
#' @name mv_marginal.MvGaussianDistrib
#'
#' @description
#' Returns the marginal law of a subset of coordinates, which for a gaussian is
#' again a gaussian: the mean is the subvector and the covariance is the
#' corresponding block, without any integration. The result is a fresh
#' [mvgaussian1_distrib()] of the reduced dimension (an [mvgaussian2_distrib()]
#' when the parent is written on the precision), whose parameters are its own
#' and need not be a subset of the ones passed in.
#'
#' @details
#' The marginal is returned on an unstructured log-Cholesky matrix whatever the
#' parent carried, because a block of a structured matrix need not have the
#' structure of the parent: the leading block of an AR(1) is AR(1), but a block
#' taken at scattered indices need not be. A parent written on the precision
#' gives a marginal written on the precision. The covariance block is computed
#' first and then inverted, since a block of a precision matrix is not the
#' inverse of the corresponding covariance block. The cost is `k(k+1)/2` free
#' values where the parent may have spent fewer.
#'
#' @param distrib An [MvGaussianDistrib] object, from [mvgaussian1_distrib()] or
#'   [mvgaussian2_distrib()].
#' @param theta A named list of parameters in which each component is a single
#'   number.
#' @param which An integer vector of distinct coordinates to keep, between 1 and
#'   \eqn{p}. The generic [mv_marginal()] checks it before dispatch.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list with `distrib`, an `MvGaussianDistrib` of dimension
#'   `length(which)` on a log-Cholesky matrix (a covariance, or a precision when
#'   the parent is written on the precision), and `theta`, its parameters as a
#'   named list.
#'
#' @seealso [mv_marginal.MvStudentTDistrib()], where the degrees of freedom are
#'   carried across unchanged, [plot.multivariate_distrib()], whose density
#'   panels are these marginals, and [mv_marginal()] for the generic.
#'
#' @examples
#' d <- mvgaussian1_distrib(3)
#' theta <- as.list(stats::setNames(
#'   c(1, -2, 0.5, 0.1, -0.2, 0.3, 0.4, -0.1, 0.2), d@params))
#'
#' m <- mv_marginal(d, theta, c(1, 3))
#' m$distrib@n_dim
#'
#' # The covariance block is carried across exactly.
#' all.equal(mv_sigma(m$distrib, m$theta),
#'           mv_sigma(d, theta)[c(1, 3), c(1, 3)], check.attributes = FALSE)
#'
#' # A single coordinate is the univariate normal, against stats::dnorm.
#' m1 <- mv_marginal(d, theta, 2)
#' s <- sqrt(mv_sigma(d, theta)[2, 2])
#' c(ours = distrib_pdf(m1$distrib, -1.4, m1$theta, log = TRUE),
#'   dnorm = dnorm(-1.4, mean = -2, sd = s, log = TRUE))
#'
#' # A precision parametrization marginalizes to a precision, computed from
#' # the covariance block, since the block of the precision is not the
#' # precision of the block.
#' o <- mvgaussian2_distrib(3, parameters7::log_cholesky(3))
#' th_o <- as.list(stats::setNames(unlist(theta), o@params))
#' mo <- mv_marginal(o, th_o, c(1, 3))
#' mo$distrib@inverted
#' all.equal(mv_sigma(mo$distrib, mo$theta),
#'           mv_sigma(o, th_o)[c(1, 3), c(1, 3)], check.attributes = FALSE)
#'
#' @keywords internal
S7::method(mv_marginal, MvGaussianDistrib) <- function(distrib, theta, which, ...) {
  mu <- as.numeric(mv_location(distrib, theta))[which]
  sg <- mv_sigma(distrib, theta)[which, which, drop = FALSE]
  # the marginal keeps the parametrization it came from: a family written in
  # the precision reports a precision, its free values read off the inverse of
  # the covariance submatrix
  md <- if (isTRUE(distrib@inverted)) {
    mvgaussian2_distrib(length(which))
  } else {
    mvgaussian1_distrib(length(which))
  }
  eta <- parameters7::param_free(
    md@param, if (isTRUE(distrib@inverted)) solve(unname(sg)) else unname(sg)
  )
  list(
    distrib = md,
    theta = as.list(stats::setNames(c(mu, unname(eta)), md@params))
  )
}
