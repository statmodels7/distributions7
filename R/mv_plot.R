#' @include multivariate.R mvgaussian_distrib.R
NULL


#' @title Panels of a Multivariate Density
#' @name plot.multivariate_distrib
#'
#' @description
#' Draws a multivariate distribution as a matrix of panels: the marginal
#' density of each coordinate on the diagonal, contours of each bivariate
#' marginal below it, and the implied correlation printed above it. The picture
#' is built entirely from marginals, so it exists exactly where
#' [mv_marginal()] does, which for the two families that ship is everywhere.
#'
#' @details
#' # Why at most three coordinates
#'
#' A panel matrix of \eqn{k} coordinates has \eqn{k^2} panels, and past three
#' the panels are too small to read. A larger distribution is rejected with an
#' error naming the count and suggesting `which`; drawing it illegibly would be
#' the worse answer.
#'
#' # What the contours are levels of
#'
#' The contours are levels of the DENSITY, not of the probability. An
#' equal-density contour is what the shape of the density means, and a
#' probability level would need the orthant integral this package does not
#' approximate in several dimensions. The levels are chosen from the density's
#' own range on the grid.
#'
#' # One setting only
#'
#' A univariate `plot()` method reads a `theta` component of length \eqn{k} as
#' \eqn{k} settings drawn over one another. Here the picture is already a matrix
#' of panels, with no axis left to overlay them on, so a component longer than
#' one is rejected by name.
#'
#' # The range each panel is drawn over
#'
#' Two and a half standard deviations either side of the location, with the
#' spread taken from the matrix the parametrization carries. A multivariate
#' Student t at \eqn{\nu \le 2} has no variance at all, and that is precisely
#' the shape worth drawing, so [variance()] cannot be the source. For an
#' elliptical family the two agree up to a factor, which is all a plotting
#' range needs.
#'
#' @param x An object inheriting from [multivariate_distrib()].
#' @param theta A named list or vector of parameters, each component a single
#'   number. When missing, parameters are drawn by
#'   [generate_random_theta()] and reported in a message. A component of
#'   length greater than one is an error.
#' @param which An integer vector of the coordinates to show, at most three.
#'   Defaults to `NULL`, which means all of them and is itself an error above
#'   three coordinates.
#' @param n_grid Points per axis for the marginal densities and the contours.
#'   Defaults to `80`. The bivariate panels evaluate the density on an
#'   `n_grid * n_grid` grid, so the cost is quadratic in it.
#' @param col_fit Color of the density curves, the contours and the printed
#'   correlations. Defaults to `"#B22222"`.
#' @param ... Passed to the underlying plotting calls.
#'
#' @return `x`, invisibly. Called for the plot it draws.
#'
#' @seealso [mv_marginal()], which supplies every panel,
#'   [plot.distrib_fit()], which draws the same matrix with the data on it,
#'   and [plot.continuous_distrib()] for the one-dimensional case, which does
#'   overlay several settings.
#'
#' @examples
#' d <- mvgaussian_distrib(3)
#' theta <- as.list(stats::setNames(
#'   c(0, 1, -1, 0, 0, 0, 0.6, -0.3, 0.2), d@params
#' ))
#'
#' op <- graphics::par(no.readonly = TRUE)
#' plot(d, theta)
#' graphics::par(op)
#'
#' # Two coordinates of a heavy-tailed family, whose contours are wider than a
#' # gaussian's at the same matrix.
#' t2 <- mvstudent_t_distrib(3)
#' th2 <- as.list(stats::setNames(c(unlist(theta), 3), t2@params))
#' op <- graphics::par(no.readonly = TRUE)
#' plot(t2, th2, which = c(1, 2))
#' graphics::par(op)
#'
#' # Four coordinates would be sixteen panels, so it is refused by name.
#' d4 <- mvgaussian_distrib(4)
#' th4 <- as.list(stats::setNames(rep(0, d4@n_params), d4@params))
#' try(plot(d4, th4))
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
  mv_pairs_panels(x, theta, which, n_grid, col_fit, data = NULL, ...)
  invisible(x)
}


#' @title Draw the Panel Matrix of a Multivariate Density
#'
#' @description
#' The shared engine of [plot.multivariate_distrib()] and the multivariate
#' branch of [plot.distrib_fit()]. It lays out one panel per ordered pair of
#' the chosen coordinates: the marginal density on the diagonal, contours of
#' the bivariate marginal below it, and the implied correlation printed above
#' it. Every panel comes from [mv_marginal()], so nothing is integrated.
#'
#' @details
#' When `data` is supplied the diagonal also carries a kernel density estimate
#' of the observed coordinate and the off-diagonal panels carry the
#' observations, so the fitted shape and the sample are read against each other
#' in one frame. The kernel estimate is the comparison that assumes no model:
#' it is what the data say without the family's help.
#'
#' The range each coordinate is drawn over is two and a half standard
#' deviations either side of the location, with the spread taken from
#' [mv_sigma()], so a Student t at \eqn{\nu \le 2} still has a range where
#' [variance()] would give it none.
#'
#' @param d A [multivariate_distrib()] object.
#' @param theta A named list of parameters, already aligned by the caller.
#' @param which An integer vector of coordinates, or `NULL` for all of them.
#'   More than three is an error.
#' @param n_grid Points per axis, a single positive whole number.
#' @param col_fit Color of the fitted density, the contours and the printed
#'   correlations.
#' @param data An \eqn{n \times p} matrix of observations, or `NULL` for a
#'   distribution with no data behind it.
#' @param col_data Color of the observed summary, used only when `data` is
#'   given. Defaults to `"#4682B4"`.
#' @param ... Unused.
#'
#' @return `NULL`, invisibly. Called for the plot it draws, and it changes
#'   `graphics::par("mfrow")` and the margins, so a caller who needs the
#'   previous settings back saves them first.
#'
#' @seealso [plot.multivariate_distrib()] and [plot.distrib_fit()], its two
#'   callers, and [mv_marginal()] for the panels.
#'
#' @examples
#' d <- mvgaussian_distrib(2)
#' theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0,
#'               sigma_log_L2 = 0, sigma_L2.1 = 0.6)
#'
#' # Without data: the density alone.
#' op <- graphics::par(no.readonly = TRUE)
#' distributions7:::mv_pairs_panels(d, theta, NULL, 40, "#B22222", data = NULL)
#' graphics::par(op)
#'
#' # With data: the same frame carries a kernel estimate and the observations.
#' set.seed(1)
#' op <- graphics::par(no.readonly = TRUE)
#' distributions7:::mv_pairs_panels(d, theta, NULL, 40, "#B22222",
#'                                  data = distrib_rng(d, 200, theta))
#' graphics::par(op)
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
  center <- as.numeric(mv_location(d, theta))
  spread <- sqrt(diag(mv_sigma(d, theta)))
  rng <- lapply(which, function(j) {
    lo <- center[j] - 3 * spread[j]
    hi <- center[j] + 3 * spread[j]
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
#'
#' @description
#' Returns the marginal law of a subset of coordinates, which for a gaussian is
#' again a gaussian: the mean is the subvector and the covariance the
#' corresponding block, with nothing to integrate. The result is a fresh
#' [mvgaussian_distrib()] of the reduced dimension, whose parameters are its
#' own: a caller cannot expect them to be a subset of the ones passed in.
#'
#' @details
#' The marginal is returned on an UNSTRUCTURED covariance whatever the parent
#' carried, because a block of a structured matrix need not have the parent's
#' structure. The leading block of an AR(1) is AR(1); a block of a
#' compound-symmetry matrix taken at scattered indices need not be; and a block
#' of a precision matrix is not the inverse of the corresponding covariance
#' block at all. Returning the unstructured form is correct in every case, at
#' the cost of `k(k+1)/2` free values where the parent may have spent fewer.
#'
#' @param distrib An [MvGaussianDistrib] object, from [mvgaussian_distrib()].
#' @param theta A named list of parameters, each component a single number.
#' @param which An integer vector of coordinates to keep, between 1 and
#'   \eqn{p}. Duplicates and out-of-range values are not checked and reach the
#'   matrix subsetting.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list with `distrib`, an `MvGaussianDistrib` of dimension
#'   `length(which)` on a log-Cholesky covariance, and `theta`, its parameters
#'   as a named list.
#'
#' @seealso [mv_marginal.MvStudentTDistrib()], where the degrees of freedom are
#'   carried across unchanged, [plot.multivariate_distrib()], whose panels are
#'   these marginals, and [mv_marginal()] for the generic.
#'
#' @examples
#' d <- mvgaussian_distrib(3)
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
#' # And a single coordinate is the univariate normal, against stats::dnorm.
#' m1 <- mv_marginal(d, theta, 2)
#' s <- sqrt(mv_sigma(d, theta)[2, 2])
#' c(ours = distrib_pdf(m1$distrib, -1.4, m1$theta, log = TRUE),
#'   dnorm = dnorm(-1.4, mean = -2, sd = s, log = TRUE))
#'
#' # A precision parametrization marginalizes to a covariance, the block of
#' # the precision not being the precision of the block.
#' o <- mvgaussian_distrib(3, omega = parameters7::log_cholesky(3))
#' th_o <- as.list(stats::setNames(unlist(theta), o@params))
#' all.equal(mv_sigma(mv_marginal(o, th_o, c(1, 3))$distrib,
#'                    mv_marginal(o, th_o, c(1, 3))$theta),
#'           mv_sigma(o, th_o)[c(1, 3), c(1, 3)], check.attributes = FALSE)
#'
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
