#' @include distrib.R generics.R
NULL

#' @title Print Method for `distrib` Objects
#' @name print.distrib
#'
#' @description
#' Writes what a distribution object is: its name and dimension, the support,
#' and one row per parameter giving the parameter's name, what it means, the
#' interval it lives in and the link that carries it to the unconstrained
#' scale. A wrapper prints its parent's line first, so
#' `truncated(gamma2_distrib(), upper = 5)` shows both the truncation and the
#' family under it.
#'
#' What it does **not** show is a value: a `distrib` carries a parametrization,
#' not an estimate. For estimates see [print.distrib_fit()].
#'
#' @param x An object inheriting from `distrib`.
#' @param ... Unused, accepted for compatibility with [base::print()].
#'
#' @return `x`, invisibly. Called for the output it writes.
#'
#' @examples
#' gaussian1_distrib()
#' poisson_distrib()
#'
#' # The link is part of the object, so changing it changes what prints.
#' poisson_distrib(link_mu = linkfunctions7::sqrt_link())
#'
#' @seealso [print.distrib_fit()] for a fitted object;
#'   [check_distrib()] to validate one;
#'   [plot.continuous_distrib()] to draw one.
#' @keywords internal
S7::method(print, distrib) <- function(x, ...) {
  # Format the distribution name (Capitalize first letters)
  d_name <- gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", x@distrib_name, perl = TRUE)
  
  # Determine type based on S7 inheritance
  d_type <- if (S7::S7_inherits(x, continuous_distrib)) {
    "Continuous"
  } else if (S7::S7_inherits(x, discrete_distrib)) {
    "Discrete"
  } else if (S7::S7_inherits(x, multivariate_distrib)) {
    sprintf("Continuous, %d-dimensional", x@n_dim)
  } else {
    "Unknown"
  }
  
  cat(sprintf("Distribution: %s\n", d_name))
  cat(sprintf("Type:         %s\n", d_type))
  cat(sprintf("Dimensions:   %s\n", x@dimension))
  
  cat("\nParameters:\n")
  
  max_param_len <- max(nchar(x@params), na.rm = TRUE)
  smooth <- param_smoothness(x)

  for (i in seq_len(x@n_params)) {
    p_name <- x@params[i]

    interpretation <- if (!is.na(x@params_interpretation[p_name])) {
      paste0("(", x@params_interpretation[[p_name]], ")")
    } else {
      ""
    }

    link_obj <- x@link_params[[p_name]]
    link_name <- if (is.null(link_obj)) "none" else link_obj@link_name

    bounds <- x@params_bounds[[p_name]]
    domain_str <- if (is.null(bounds)) "(unspecified)" else paste0("(", bounds[1], ", ", bounds[2], ")")

    smooth_str <- if (isFALSE(smooth[[p_name]])) "  [non-smooth log-likelihood]" else ""

    cat(sprintf(
      "  %-*s %-20s | Link: %-10s | Domain: %s%s\n",
      max_param_len, p_name,
      interpretation,
      link_name,
      domain_str,
      smooth_str
    ))
  }

  invisible(x)
}

#' @title Random Parameters Inside a Family's Own Domain
#' @name generate_random_theta.distrib
#'
#' @description
#' Draws one parameter value per entry of `distrib@params`, each inside that
#' parameter's own interval in `params_bounds` and strictly away from its
#' ends. A bounded interval is sampled across its width, a half-line is
#' sampled from a positive draw offset from the finite end, and an unbounded
#' parameter is sampled around zero.
#'
#' The draws ignore the response entirely, so they are a fallback and never an
#' estimate: being of order one whatever the data, they suit a shape or a
#' probability and are the wrong order for a location or a scale on a response
#' of any size. [distrib_start()] uses them for the starting values
#' after the first, and [check_distrib()] uses them to probe a family at a
#' parameter nobody chose.
#'
#' @param distrib An object inheriting from `distrib`, read for `params` and
#'   `params_bounds`.
#' @param ... Unused.
#'
#' @return A named list with one numeric value per parameter, named and
#'   ordered as `distrib@params`, every value strictly inside its own bounds.
#'
#' @examples
#' set.seed(1)
#' unlist(generate_random_theta(gaussian1_distrib()))
#'
#' # Every draw is admissible by construction: a scale is positive, a
#' # probability is inside (0, 1).
#' set.seed(2)
#' d <- beta1_distrib()
#' th <- generate_random_theta(d)
#' rbind(draw = unlist(th), lower = vapply(d@params_bounds, `[`, numeric(1), 1),
#'       upper = vapply(d@params_bounds, `[`, numeric(1), 2))
#'
#' @seealso [distrib_start()], which draws these for the starts after the
#'   first; [start_from_moments()] for the data-based one;
#'   [check_distrib()], which probes a family at a random parameter.
#' @keywords internal
S7::method(generate_random_theta, distrib) <- function(distrib, ...) {
  theta <- list()
  for (p in distrib@params) {
    b <- distrib@params_bounds[[p]]
    lower <- b[1]
    upper <- b[2]
    
    if (is.finite(lower) && is.finite(upper)) {
      span <- upper - lower
      theta[[p]] <- stats::runif(1, min = lower + 0.1 * span, max = upper - 0.1 * span)
    } else if (is.finite(lower) && is.infinite(upper)) {
      theta[[p]] <- lower + stats::runif(1, min = 0.1, max = 6.0)
    } else if (is.infinite(lower) && is.finite(upper)) {
      theta[[p]] <- upper - stats::runif(1, min = 0.1, max = 6.0)
    } else {
      theta[[p]] <- stats::runif(1, min = -3.0, max = 3.0)
    }
  }
  theta
}

#' @title One Parameter Setting Per Curve
#'
#' @description
#' Turns the `theta` a plot method was given into the list of scalar settings
#' it draws, one per curve. A component of length \eqn{k > 1} asks for \eqn{k}
#' curves, so `list(mu = 0, sigma = c(1, 2, 4))` becomes three settings sharing
#' a mean and differing in scale.
#'
#' @details
#' Every component must have length 1 or the same \eqn{k}. A length that merely
#' **divides** \eqn{k} signals an error: R would recycle it without a word, and
#' a partial setting is far likelier to be a mistake than a request.
#'
#' The meaning is available only because a plot has no data to recycle a
#' parameter against. [distrib_pdf()] and every derivative generic read a
#' vector component as one value per observation, which is why each setting is
#' handed to them as scalars.
#'
#' @param x A distribution object, read for `params`.
#' @param theta A named list, or a named numeric vector, holding one entry per
#'   parameter. A missing parameter signals an error naming what was expected;
#'   so does a component of length zero.
#'
#' @return A list with three components: `settings`, a list of \eqn{k} named
#'   parameter lists each holding scalars; `k`, the number of curves; and
#'   `varying`, a character vector naming the parameters that differ between
#'   settings, which may be empty.
#'
#' @seealso [plot.continuous_distrib()] and [plot.discrete_distrib()], the
#'   callers; [plot_labels()], which turns `varying` into a legend;
#'   [plot_keys()] for the colors and symbols.
#' @keywords internal
plot_settings <- function(x, theta) {
  if (is.numeric(theta) && !is.list(theta)) theta <- as.list(theta)
  if (!all(x@params %in% names(theta))) {
    stop(sprintf("Missing parameters in 'theta'. Expected: %s",
                 paste(x@params, collapse = ", ")), call. = FALSE)
  }
  theta <- theta[x@params]
  lens <- vapply(theta, length, integer(1))
  if (any(lens == 0L)) {
    stop(sprintf("Empty parameters in 'theta': %s",
                 paste(x@params[lens == 0L], collapse = ", ")), call. = FALSE)
  }
  k <- max(lens)
  bad <- x@params[lens != 1L & lens != k]
  if (length(bad)) {
    stop(sprintf(paste0("Every parameter must have length 1 or %d, the number ",
                        "of curves asked for. %s"),
                 k, paste(sprintf("'%s' has length %d", bad, lens[bad]),
                          collapse = ", ")), call. = FALSE)
  }
  settings <- lapply(seq_len(k), function(i) {
    lapply(theta, function(v) if (length(v) == 1L) v[[1L]] else v[[i]])
  })
  list(settings = settings, k = k,
       varying = if (k > 1L) x@params[lens == k] else character(0))
}

#' @title Colors, Line Types and Symbols for a Set of Curves
#'
#' @description
#' Returns the three visual keys a plot method distinguishes its settings by,
#' each recycled to length \eqn{k}. A single curve is black, solid and a filled
#' circle. Several take a qualitative palette, the six base line types and six
#' distinguishable symbols, so that a printed copy with no color is still
#' readable.
#'
#' A continuous family draws with the color and the line type; a discrete one
#' draws with the color and the **symbol**, because dashing a stem is what a
#' line type does there and a dashed stem reads as a broken one.
#'
#' @param k How many settings are drawn, a single positive integer.
#' @param dots The caller's `...`, as a list. A `col`, `lty` or `pch` given
#'   there wins and is recycled over the settings; anything else is ignored.
#'
#' @return A list with components `col`, `lty` and `pch`, each of length `k`.
#'
#' @seealso [plot.continuous_distrib()] and [plot.discrete_distrib()], which
#'   read different pairs of these; [plot_settings()] for `k`.
#' @keywords internal
plot_keys <- function(k, dots = list()) {
  col <- if (!is.null(dots$col)) rep_len(dots$col, k)
         else if (k == 1L) "black"
         else grDevices::hcl.colors(k, palette = "Dark 3")
  lty <- if (!is.null(dots$lty)) rep_len(dots$lty, k)
         else if (k == 1L) 1L
         else rep_len(seq_len(6L), k)
  pch <- if (!is.null(dots$pch)) rep_len(dots$pch, k)
         else if (k == 1L) 16L
         else rep_len(c(16L, 17L, 15L, 18L, 1L, 2L), k)
  list(col = col, lty = lty, pch = pch)
}

#' @title Where to Put the Key
#'
#' @description
#' Returns the corner of the panel a legend goes in: the one the curves leave
#' emptier. A right-skewed family puts its mass on the left, so a key fixed at
#' the top right never meets it, while a left-skewed one is covered by exactly
#' that choice.
#'
#' The side is read off the drawing. The density's center of mass across all
#' the settings is compared with the middle of the horizontal range, and the
#' key goes to the other side.
#'
#' @param y The evaluation points, a numeric vector.
#' @param dens A list of densities, one per setting, each the length of `y`.
#'   Non-finite and negative entries are treated as zero.
#'
#' @return `"topright"` or `"topleft"`, a character string of length 1.
#'   `"topright"` where no setting has any positive mass.
#'
#' @seealso [plot.continuous_distrib()] and [plot.discrete_distrib()], the
#'   callers; [plot_labels()] for what goes in the key.
#' @keywords internal
plot_legend_side <- function(y, dens) {
  w <- do.call(pmax, c(dens, list(na.rm = TRUE)))
  w[!is.finite(w) | w < 0] <- 0
  if (!any(w > 0)) return("topright")
  centre <- sum(y * w) / sum(w)
  if (centre < mean(range(y))) "topright" else "topleft"
}

#' @title The Legend Entries and the Title of a Distribution Plot
#'
#' @description
#' Splits a set of plot settings into what varies and what does not: the
#' varying parameters become one legend entry per curve, and the fixed ones are
#' stated once in the title beside the family's name. A reader then sees which
#' parameter the panel is about without counting curves, and reads the held
#' values without a second key.
#'
#' Where nothing varies there is one curve and no legend, and the title carries
#' every parameter. Where nothing is fixed the title is the family's name
#' alone. Values are rounded to three decimals.
#'
#' @param x A distribution object, read for `distrib_name` and `params`.
#' @param ps The value of [plot_settings()], read for `settings` and
#'   `varying`.
#'
#' @return A list with two components: `legend`, a character vector of length
#'   `ps$k` or `NULL` when nothing varies, and `main`, a character string of
#'   length 1 carrying a newline between the family's name and the fixed
#'   values.
#'
#' @seealso [plot_settings()], which supplies `ps`;
#'   [plot_legend_side()] for where the legend goes;
#'   [plot.continuous_distrib()] and [plot.discrete_distrib()], the callers.
#' @keywords internal
plot_labels <- function(x, ps) {
  fmt <- function(v) format(round(v, 3), trim = TRUE)
  fixed <- setdiff(x@params, ps$varying)
  if (!length(ps$varying)) {
    return(list(
      legend = NULL,
      main = sprintf("%s distribution\n(%s)", x@distrib_name,
                     paste(fixed, fmt(unlist(ps$settings[[1L]][fixed])),
                           sep = " = ", collapse = ", "))))
  }
  legend <- vapply(ps$settings, function(th) {
    paste(ps$varying, fmt(unlist(th[ps$varying])), sep = " = ",
          collapse = ", ")
  }, character(1))
  main <- if (length(fixed)) {
    sprintf("%s distribution\n(%s)", x@distrib_name,
            paste(fixed, fmt(unlist(ps$settings[[1L]][fixed])),
                  sep = " = ", collapse = ", "))
  } else {
    sprintf("%s distribution", x@distrib_name)
  }
  list(legend = legend, main = main)
}

#' @title Plot Method for Continuous Distributions
#'
#' @name plot.continuous_distrib
#'
#' @description
#' Draws the density of a continuous family over a grid of 1000 points, one
#' curve per parameter setting. A component of `theta` given as a vector asks
#' for one curve per element, so `plot(d, list(mu = 0, sigma = c(1, 2, 4)))`
#' draws three densities that share a mean and differ in scale.
#'
#' @details
#' # The two keys
#' Curves are separated by **color and line type together**, so that a printed
#' copy with no color is still readable. The parameters that vary are named in
#' a legend and those held fixed are stated once in the title, so a reader sees
#' what the panel is about without counting curves. The legend goes in
#' whichever top corner the mass leaves emptier. See
#' [plot_settings()] for the rule on lengths and [plot_keys()] for the palette.
#'
#' # The horizontal range
#' It covers every setting, not the first: from the smallest 0.5% quantile to
#' the largest 99.5% quantile over them, widened by a tenth and clamped to the
#' support. Give `xlim` to override it, which is worth doing for a heavy-tailed
#' family, whose 99.5% quantile can sit orders of magnitude past anything worth
#' drawing.
#'
#' @param x An object inheriting from `continuous_distrib`.
#' @param theta A named list or numeric vector holding one entry per parameter
#'   of `x`. A component of length \eqn{k > 1} draws \eqn{k} curves, and every
#'   component must have length 1 or that same \eqn{k}. Missing, a random
#'   parameter is drawn by [generate_random_theta()].
#' @param xlim A numeric vector of length 2, or `NULL` (the default) to compute
#'   the range from the quantiles as above.
#' @param legend Logical of length 1: draw the key when several settings are
#'   plotted. Defaults to `TRUE`. Nothing is drawn when only one is.
#' @param ... Passed to [graphics::plot()], for instance `main`, `xlab` or
#'   `ylim`. A `col` or `lty` given here wins over the palette and is recycled
#'   over the curves.
#'
#' @return `x`, invisibly. Called for the plot it draws.
#'
#' @examples
#' # One curve per value of the scale, at a shared mean.
#' plot(gaussian1_distrib(), list(mu = 0, sigma = c(1, 2, 4)))
#'
#' # A heavy-tailed family wants its own window: the 99.5% quantile of a
#' # Cauchy is two orders of magnitude past the interesting part.
#' plot(cauchy_distrib(), list(mu = 0, sigma = 1), xlim = c(-6, 6))
#'
#' @seealso [plot.discrete_distrib()] for a lattice family;
#'   [plot.distrib_fit()] to compare a fit with data;
#'   [plot_settings()], [plot_keys()] and [plot_labels()] for the pieces.
#' @importFrom graphics plot
#' @keywords internal
S7::method(plot, continuous_distrib) <- function(x, theta, xlim = NULL,
                                                 legend = TRUE, ...) {
  if (missing(theta)) {
    theta <- generate_random_theta(x)
    param_str <- paste(names(theta), round(unlist(theta), 3), sep = " = ", collapse = ", ")
    message(sprintf("Argument 'theta' is missing. Automatically generating random parameters: %s", param_str))
  }
  ps <- plot_settings(x, theta)
  labs <- plot_labels(x, ps)

  if (is.null(xlim)) {
    # the window has to hold every setting, so the quantiles are taken over
    # all of them and not over the first
    qs <- vapply(ps$settings, function(th) {
      c(suppressWarnings(distrib_quantile(x, 0.005, th)),
        suppressWarnings(distrib_quantile(x, 0.995, th)))
    }, numeric(2))
    lower_q <- min(qs[1L, ], na.rm = TRUE)
    upper_q <- max(qs[2L, ], na.rm = TRUE)

    span <- max(upper_q - lower_q, 1.0, na.rm = TRUE)

    # Clamp limits strictly to mathematical bounds to prevent evaluating out-of-domain points
    xlim <- c(
      max(lower_q - 0.1 * span, x@bounds[1], na.rm = TRUE),
      min(upper_q + 0.1 * span, x@bounds[2], na.rm = TRUE)
    )
  }

  seq_x <- seq(xlim[1], xlim[2], length.out = 1000)
  dens <- lapply(ps$settings,
                 function(th) distrib_pdf(x, seq_x, th, log = FALSE))

  dots <- list(...)
  keys <- plot_keys(ps$k, dots)
  if (is.null(dots$main)) dots$main <- labs$main
  if (is.null(dots$xlab)) dots$xlab <- "y"
  if (is.null(dots$ylab)) dots$ylab <- "Density (PDF)"
  if (is.null(dots$lwd)) dots$lwd <- 2
  if (is.null(dots$ylim)) {
    dots$ylim <- c(0, max(unlist(dens)[is.finite(unlist(dens))], na.rm = TRUE))
  }
  dots$col <- NULL
  dots$lty <- NULL

  plot_args <- c(list(x = seq_x, y = dens[[1L]], xlim = xlim, type = "n"), dots)
  do.call(graphics::plot, plot_args)
  for (i in seq_len(ps$k)) {
    graphics::lines(seq_x, dens[[i]], col = keys$col[i], lty = keys$lty[i],
                    lwd = dots$lwd)
  }
  if (isTRUE(legend) && !is.null(labs$legend)) {
    graphics::legend(plot_legend_side(seq_x, dens), legend = labs$legend,
                     col = keys$col, lty = keys$lty, lwd = dots$lwd,
                     bty = "n")
  }

  invisible(x)
}

#' @title Plot Method for Discrete Distributions
#'
#' @name plot.discrete_distrib
#'
#' @description
#' Draws the probability mass of a discrete family as stems over the integers
#' in range, one set per parameter setting. A component of `theta` given as a
#' vector asks for one setting per element, exactly as in
#' [plot.continuous_distrib()].
#'
#' @details
#' # The two keys, and why they are not the continuous ones
#' Settings are separated by **color and symbol together**, the symbol sitting
#' at the top of each stem where the point already is. The line type is not
#' used here: dashing a stem is what it would do, and a dashed stem reads as a
#' broken one, so at a support of any size the panel fills with fragments that
#' cross each other.
#'
#' Several settings are also shifted sideways by a fraction of the spacing, so
#' that equal masses at one support point stay countable instead of standing
#' one in front of another.
#'
#' @param x An object inheriting from `discrete_distrib`.
#' @param theta A named list or numeric vector holding one entry per parameter
#'   of `x`. A component of length \eqn{k > 1} draws \eqn{k} sets of stems, and
#'   every component must have length 1 or that same \eqn{k}. Missing, a random
#'   parameter is drawn by [generate_random_theta()].
#' @param xlim A numeric vector of length 2, or `NULL` (the default) to compute
#'   the range from the quantiles over every setting.
#' @param legend Logical of length 1: draw the key when several settings are
#'   plotted. Defaults to `TRUE`.
#' @param ... Passed to [graphics::plot()], for instance `main`, `xlab` or
#'   `ylim`. A `col` or `pch` given here wins over the palette and is recycled
#'   over the settings; an `lty` is accepted and has no effect on the stems.
#'
#' @return `x`, invisibly. Called for the plot it draws.
#'
#' @examples
#' # Three rates, whose masses are countable at every support point because
#' # the settings are offset sideways.
#' plot(poisson_distrib(), list(mu = c(1, 4, 10)))
#'
#' # Overdispersion at a fixed mean: the mass spreads as theta falls.
#' plot(negbin2_distrib(), list(mu = 5, theta = c(0.5, 2, 50)))
#'
#' @seealso [plot.continuous_distrib()] for a density;
#'   [plot.distrib_fit()] to compare a fit with data;
#'   [plot_settings()], [plot_keys()] and [plot_labels()] for the pieces.
#' @importFrom graphics plot segments points
#' @keywords internal
S7::method(plot, discrete_distrib) <- function(x, theta, xlim = NULL,
                                               legend = TRUE, ...) {
  if (missing(theta)) {
    theta <- generate_random_theta(x)
    param_str <- paste(names(theta), round(unlist(theta), 3), sep = " = ", collapse = ", ")
    message(sprintf("Argument 'theta' is missing. Automatically generating random parameters: %s", param_str))
  }
  ps <- plot_settings(x, theta)
  labs <- plot_labels(x, ps)

  if (is.null(xlim)) {
    qs <- vapply(ps$settings, function(th) {
      c(suppressWarnings(distrib_quantile(x, 0.005, th)),
        suppressWarnings(distrib_quantile(x, 0.995, th)))
    }, numeric(2))
    xlim <- c(
      max(floor(min(qs[1L, ], na.rm = TRUE)), x@bounds[1], na.rm = TRUE),
      min(ceiling(max(qs[2L, ], na.rm = TRUE)), x@bounds[2], na.rm = TRUE)
    )
  }

  # Generate strict integer sequence clamped to mathematical bounds
  seq_x <- seq(floor(xlim[1]), ceiling(xlim[2]), by = 1)
  seq_x <- seq_x[seq_x >= x@bounds[1] & seq_x <= x@bounds[2]]

  dens <- lapply(ps$settings,
                 function(th) distrib_pdf(x, seq_x, th, log = FALSE))

  dots <- list(...)
  keys <- plot_keys(ps$k, dots)
  if (is.null(dots$main)) dots$main <- labs$main
  if (is.null(dots$xlab)) dots$xlab <- "y"
  if (is.null(dots$ylab)) dots$ylab <- "Probability (PMF)"
  if (is.null(dots$lwd)) dots$lwd <- 2
  if (is.null(dots$ylim)) {
    dots$ylim <- c(0, max(unlist(dens)[is.finite(unlist(dens))], na.rm = TRUE))
  }
  dots$col <- NULL
  dots$lty <- NULL
  dots$pch <- NULL

  plot_args <- c(list(x = seq_x, y = dens[[1L]], xlim = xlim, type = "n"), dots)
  do.call(graphics::plot, plot_args)

  # A shift of a fraction of the unit spacing, so that equal masses at one
  # support point stay countable. The stems are drawn solid whatever the
  # setting and told apart by their symbol: a dashed stem reads as a broken
  # one, and at a support of any size the panel fills with fragments.
  shift <- if (ps$k > 1L) (seq_len(ps$k) - (ps$k + 1) / 2) * (0.55 / ps$k)
           else 0
  stem_lwd <- if (ps$k > 1L) max(1, dots$lwd - 1) else dots$lwd
  for (i in seq_len(ps$k)) {
    at <- seq_x + shift[i]
    graphics::segments(x0 = at, y0 = 0, x1 = at, y1 = dens[[i]],
                       col = keys$col[i], lwd = stem_lwd)
    graphics::points(at, dens[[i]], pch = keys$pch[i], col = keys$col[i],
                     cex = if (ps$k > 1L) 0.65 else 1)
  }
  if (isTRUE(legend) && !is.null(labs$legend)) {
    graphics::legend(plot_legend_side(seq_x, dens), legend = labs$legend,
                     col = keys$col, pch = keys$pch, bty = "n")
  }

  invisible(x)
}
