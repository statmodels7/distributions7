#' @include multivariate.R mvgaussian_distrib.R
NULL


#' Panels of a Multivariate Density
#'
#' @name plot.multivariate_distrib
#'
#' @description
#' Draws a multivariate distribution as a matrix of panels: the marginal
#' density of each coordinate on the diagonal, and contours of each bivariate
#' marginal below it.
#'
#' @details
#' The picture is built from marginals, so it exists exactly when
#' \code{\link{mv_marginal}} does. Above about three coordinates the panel
#' matrix stops being readable, so a larger distribution is refused with the
#' suggestion of choosing coordinates rather than being drawn illegibly.
#'
#' The contours are drawn at levels of the density itself rather than at
#' probability levels: the equal-density contour is what the density's shape
#' means, and computing a probability level would need the integral this
#' package refuses to approximate in several dimensions.
#'
#' @param x An object inheriting from class \code{\link{multivariate_distrib}}.
#' @param theta A named list or vector of parameters. Generated at random when
#'   missing, as for a univariate distribution.
#' @param which The coordinates to show. Defaults to all of them.
#' @param n_grid Points per axis for the marginal densities and the contours.
#' @param col_fit Colour of the density and of the contours.
#' @param ... Passed to the underlying plotting calls.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{mv_marginal}}, \code{\link{plot.distrib_fit}}
#'
#' @examples
#' d <- mvgaussian_distrib(3)
#' theta <- as.list(stats::setNames(
#'   c(0, 1, -1, 0, 0, 0, 0.6, -0.3, 0.2), d@params
#' ))
#' plot(d, theta)
#'
#' @keywords internal
S7::method(plot, multivariate_distrib) <- function(x, theta, which = NULL,
                                                   n_grid = 80,
                                                   col_fit = "#B22222", ...) {
  if (missing(theta)) {
    theta <- generate_random_theta(x)
    message(sprintf(
      "Argument 'theta' is missing. Using random parameters: %s",
      paste(names(theta), round(unlist(theta), 3), sep = " = ", collapse = ", ")
    ))
  }
  theta <- align_theta(x, theta)
  mv_pairs_panels(x, theta, which, n_grid, col_fit, data = NULL, ...)
  invisible(x)
}


#' Draw the Panel Matrix of a Multivariate Density
#'
#' @description
#' The common engine of \code{\link{plot.multivariate_distrib}} and the
#' multivariate branch of \code{\link{plot.distrib_fit}}: one panel per pair of
#' coordinates, the marginal densities on the diagonal.
#'
#' @details
#' When \code{data} is supplied the diagonal also carries a kernel density
#' estimate of the observed coordinate and the off-diagonal panels carry the
#' observations, so that the fitted shape and the sample are read against each
#' other in the same frame. The kernel estimate is the comparison that does not
#' assume the model: it is what the data say without the family's help.
#'
#' @param d A \code{\link{multivariate_distrib}} object.
#' @param theta A named list of parameters, already aligned.
#' @param which The coordinates to show, or \code{NULL} for all.
#' @param n_grid Points per axis.
#' @param col_fit Colour of the fitted density.
#' @param data An \eqn{n \times p} matrix of observations, or \code{NULL}.
#' @param col_data Colour of the observed summary.
#' @param ... Unused.
#'
#' @return Invisibly \code{NULL}.
#'
#' @keywords internal
mv_pairs_panels <- function(d, theta, which, n_grid, col_fit, data,
                            col_data = "#4682B4", ...) {
  p <- d@n_dim
  if (is.null(which)) which <- seq_len(p)
  which <- as.integer(which)
  k <- length(which)

  if (k > 3L) {
    stop(sprintf(paste0(
      "A panel matrix of %d coordinates has %d panels and is not readable.\n",
      "  Choose at most three with 'which', for instance which = c(1, 2, 3)."
    ), k, k * k), call. = FALSE)
  }

  # The range each coordinate is drawn over, from its own marginal: two and a
  # half standard deviations either side of the mean covers the shape without
  # flattening it into the axis.
  # The spread comes from the matrix the matrix parameter carries rather than from
  # variance(): a Student t with nu <= 2 has no variance, and that is precisely
  # the shape worth drawing. The two agree up to a factor for an elliptical
  # family, which is all a plotting range needs.
  centre <- as.numeric(mv_location(d, theta))
  spread <- sqrt(diag(mv_sigma(d, theta)))
  rng <- lapply(which, function(j) {
    lo <- centre[j] - 3 * spread[j]
    hi <- centre[j] + 3 * spread[j]
    if (!is.null(data)) {
      col <- match(j, which)
      lo <- min(lo, min(data[, col]))
      hi <- max(hi, max(data[, col]))
    }
    c(lo, hi)
  })
  names(rng) <- as.character(which)

  op <- graphics::par(mfrow = c(k, k), mar = c(3.2, 3.2, 1.6, 0.8),
                      mgp = c(2, 0.6, 0))
  on.exit(graphics::par(op), add = TRUE)

  lab <- paste0("v", which)

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
        ymax <- max(dens)
        kd <- NULL
        if (!is.null(data)) {
          kd <- stats::density(data[, a])
          ymax <- max(ymax, max(kd$y))
        }
        graphics::plot(xs, dens,
          type = "l", col = col_fit, lwd = 2,
          xlab = lab[a], ylab = "density", ylim = c(0, ymax * 1.05),
          main = ""
        )
        if (!is.null(kd)) {
          graphics::lines(kd, col = col_data, lwd = 2, lty = 2)
          if (a == 1L) {
            graphics::legend("topright",
              legend = c("fitted", "kernel"), col = c(col_fit, col_data),
              lty = c(1, 2), lwd = 2, bty = "n", cex = 0.8
            )
          }
        }
      } else if (a > b) {
        m <- mv_marginal(d, theta, c(jb, ja))
        gx <- seq(rng[[as.character(jb)]][1], rng[[as.character(jb)]][2],
          length.out = n_grid
        )
        gy <- seq(rng[[as.character(ja)]][1], rng[[as.character(ja)]][2],
          length.out = n_grid
        )
        grid <- as.matrix(expand.grid(gx, gy))
        z <- matrix(distrib_pdf(m$distrib, grid, m$theta), n_grid, n_grid)
        graphics::plot(NA,
          xlim = range(gx), ylim = range(gy),
          xlab = lab[b], ylab = lab[a], main = ""
        )
        if (!is.null(data)) {
          graphics::points(data[, b], data[, a],
            pch = 16, cex = 0.35,
            col = grDevices::adjustcolor(col_data, alpha.f = 0.35)
          )
        }
        graphics::contour(gx, gy, z,
          add = TRUE, col = col_fit, lwd = 1.6,
          drawlabels = FALSE, nlevels = 6
        )
      } else {
        # The upper triangle repeats the lower one with the axes swapped, so it
        # is left empty and used for the correlation instead, which is the one
        # number a reader wants from a pair and cannot read off a contour.
        s <- mv_sigma(d, theta)
        r <- s[ja, jb] / sqrt(s[ja, ja] * s[jb, jb])
        graphics::plot(NA,
          xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
          xlab = "", ylab = "", main = ""
        )
        graphics::text(0.5, 0.5, sprintf("corr\n%.3f", r),
          cex = 1.3, col = col_fit
        )
      }
    }
  }
  invisible(NULL)
}


#' @title Marginal of a Multivariate Gaussian
#' @name mv_marginal.MvGaussianDistrib
#' @description
#' A marginal of a gaussian is a gaussian: the mean is the subvector and the
#' covariance the corresponding block, with no integration to perform.
#' @details
#' The marginal is returned on an unstructured covariance whatever the parent
#' carried, because a block of a structured matrix need not have the parent's
#' structure -- the leading block of an AR(1) is AR(1), but a block of a
#' compound-symmetry matrix taken at scattered indices need not be, and a
#' precision block is not the inverse of the covariance block at all.
#' @param distrib A \code{\link{MvGaussianDistrib}} object.
#' @param theta A named list of parameters.
#' @param which An integer vector of coordinates.
#' @param ... Unused.
#' @return A list with \code{distrib} and \code{theta}.
#' @keywords internal
S7::method(mv_marginal, MvGaussianDistrib) <- function(distrib, theta, which, ...) {
  mu <- as.numeric(mv_location(distrib, theta))[which]
  sg <- mv_sigma(distrib, theta)[which, which, drop = FALSE]
  md <- mvgaussian_distrib(length(which))
  eta <- parameters7::param_free(md@param, unname(sg))
  list(
    distrib = md,
    theta = as.list(stats::setNames(c(mu, unname(eta)), md@params))
  )
}
