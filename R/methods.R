#' @include distrib.R generics.R ggplot.R
NULL

#' @title Print Method for `distrib` Objects
#' @name print.distrib
#'
#' @description
#' Prints the name of the distribution, its type (continuous or discrete) and
#' its dimension, followed by one row per parameter giving the name of the
#' parameter, its interpretation, the interval in which it lies and the link
#' that maps it to the unconstrained scale. The name of a wrapper combines the
#' wrapper and its parent, so `truncated(gamma2_distrib(), upper = 5)` prints
#' as `Truncated Gamma2 [upper=5]`, followed by the parameter rows of the
#' parent. The method does not show a parameter value, because a `distrib`
#' carries a parametrization and not an estimate; for estimates see
#' [print.distrib_fit()].
#'
#' @param x An object inheriting from `distrib`.
#' @param ... Unused, accepted for compatibility with [base::print()].
#'
#' @return `x`, invisibly. The function is called for the output that it writes.
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
#' Draws one parameter value per entry of `distrib@params`, and each value lies
#' inside that parameter's own interval in `params_bounds` and strictly away
#' from its ends. A bounded interval is sampled across its width, a half-line is
#' sampled from a positive draw offset from the finite end, and an unbounded
#' parameter is sampled around zero.
#'
#' The draws ignore the response entirely, so they are a fallback and never an
#' estimate. They are of order one whatever the data, so they suit a shape or a
#' probability and are the wrong order for a location or a scale on a response
#' of any size. The univariate methods of [distrib_start()] use them for the
#' starting values after the first, and the base method for all of them;
#' [check_distrib()] uses them to probe a family at an arbitrary parameter.
#'
#' @param distrib An object inheriting from `distrib`, read for `params` and
#'   `params_bounds`.
#' @param ... Unused.
#'
#' @return A named list with one numeric value per parameter, named and
#'   ordered as `distrib@params`; every value lies strictly inside its own
#'   bounds.
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
#' Turns the `theta` that a plot method receives into the list of scalar
#' settings that the method draws, one per curve. A component of length
#' \eqn{k > 1} gives \eqn{k} curves, so `list(mu = 0, sigma = c(1, 2, 4))`
#' becomes three settings that share a mean and differ in scale.
#'
#' @details
#' Every component must have length 1 or the same \eqn{k}. A length that only
#' divides \eqn{k} signals an error, because R would recycle such a
#' component without a warning.
#'
#' A vector component can be read as one setting per curve only because a plot
#' has no data against which to recycle a parameter. [distrib_pdf()] and every
#' derivative generic read a vector component as one value per observation,
#' so each setting is handed to them as scalars.
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
    stop(sprintf("'theta' lacks some parameters; the expected parameters are %s.",
                 paste(x@params, collapse = ", ")), call. = FALSE)
  }
  theta <- theta[x@params]
  lens <- vapply(theta, length, integer(1))
  if (any(lens == 0L)) {
    stop(sprintf("'theta' has empty parameters: %s",
                 paste(x@params[lens == 0L], collapse = ", ")), call. = FALSE)
  }
  k <- max(lens)
  bad <- x@params[lens != 1L & lens != k]
  if (length(bad)) {
    stop(sprintf(paste0("Every parameter must have length 1 or %d, the number ",
                        "of curves requested. %s"),
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
#' Returns the three visual keys by which a graph distinguishes its settings,
#' and each key is recycled to length \eqn{k}. A single curve is drawn in the
#' main colour of the toolkit's graphs, solid and with a filled circle.
#' Several take a qualitative palette, the six line types and six
#' distinguishable symbols, so that a printed copy with no color is still
#' readable.
#'
#' A continuous family draws with the color and the line type. A discrete
#' family draws with the color and the symbol, because a line type would dash
#' the stems and a dashed stem reads as a broken one.
#'
#' @param k How many settings are drawn, a single positive integer.
#'
#' @return A list with components `col`, `lty` and `pch`, each of length `k`.
#'
#' @seealso [plot.continuous_distrib()] and [plot.discrete_distrib()], which
#'   read different pairs of these; [plot_settings()] for `k`.
#' @keywords internal
plot_keys <- function(k) {
  if (k == 1L) return(list(col = "#1f5fa8", lty = 1L, pch = 16L))
  list(col = grDevices::hcl.colors(k, palette = "Dark 3"),
       lty = rep_len(seq_len(6L), k),
       pch = rep_len(c(16L, 17L, 15L, 18L, 1L, 2L), k))
}

#' @title The Legend Entries and the Title of a Distribution Plot
#'
#' @description
#' Splits a set of plot settings into what varies and what does not: the
#' varying parameters become one legend entry per curve, and the fixed ones are
#' stated once in the title beside the family's name. A reader then sees which
#' parameter the graph is about without counting curves, and reads the held
#' values without a second key.
#'
#' If every parameter is fixed, there is one curve and no legend, and the title
#' carries every parameter. If every parameter varies, the title is only the
#' family's name. Values are rounded to three decimals.
#'
#' @param x A distribution object, read for `distrib_name` and `params`.
#' @param ps The value of [plot_settings()], read for `settings` and
#'   `varying`.
#'
#' @return A list with two components: `legend`, a character vector of length
#'   `ps$k` or `NULL` when every parameter is fixed, and `main`, a character
#'   string of length 1 carrying a newline between the family's name and the
#'   fixed values.
#'
#' @seealso [plot_settings()], which supplies `ps`;
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

#' @title The Settings of a Univariate Graph
#'
#' @description
#' Draws random parameters when `theta` is missing, reporting them in a
#' message, and returns the settings, the legend entries and the title of the
#' graph. It is the first step of the `autoplot()` methods of
#' [plot.continuous_distrib()] and [plot.discrete_distrib()].
#'
#' @param x A univariate distribution object.
#' @param theta The caller's `theta`, possibly missing.
#'
#' @return A list with the components `ps` (the value of [plot_settings()])
#'   and `labs` (the value of [plot_labels()]).
#'
#' @keywords internal
plot_univariate_setup <- function(x, theta) {
  if (missing(theta)) {
    theta <- generate_random_theta(x)
    param_str <- paste(names(theta), round(unlist(theta), 3), sep = " = ", collapse = ", ")
    message(sprintf("Argument 'theta' is missing. The plot uses random parameters: %s", param_str))
  }
  ps <- plot_settings(x, theta)
  list(ps = ps, labs = plot_labels(x, ps))
}

#' @title Plot Method for Continuous Distributions
#'
#' @name plot.continuous_distrib
#'
#' @description
#' `autoplot()` returns a graph of the density of a continuous family over a
#' grid of 1000 points, one curve per parameter setting, and `plot()` prints
#' that graph. A component of `theta` given as a vector gives one curve per
#' element, so `plot(d, list(mu = 0, sigma = c(1, 2, 4)))` draws three
#' densities that share a mean and differ in scale.
#'
#' @details
#' # The two keys
#' Curves are separated by color and line type together, so that a printed
#' copy with no color is still readable. The parameters that vary are named in
#' a legend below the panel and those held fixed are stated once in the title,
#' so a reader sees what the graph is about without counting curves. A single
#' setting is drawn without a legend, and its title carries every parameter.
#' See [plot_settings()] for the rule on lengths and [plot_keys()] for the
#' palette.
#'
#' # The horizontal range
#' It covers every setting, not the first: from the smallest 0.5% quantile to
#' the largest 99.5% quantile over them, widened by a tenth and clamped to the
#' support. Give `xlim` to override it, which is useful for a heavy-tailed
#' family, whose 99.5% quantile can sit orders of magnitude past anything worth
#' drawing.
#'
#' The graph is drawn with \pkg{ggplot2}, which is suggested; the methods
#' signal an error that names it when it is not installed. The title, the
#' labels of the axes and the colours are changed on the returned object.
#'
#' @param x,object An object inheriting from `continuous_distrib`.
#' @param theta A named list or numeric vector holding one entry per parameter
#'   of `x`. A component of length \eqn{k > 1} draws \eqn{k} curves, and every
#'   component must have length 1 or that same \eqn{k}. When `theta` is
#'   missing, random values for all parameters are drawn by
#'   [generate_random_theta()].
#' @param xlim A numeric vector of length 2, or `NULL` (the default) to compute
#'   the range from the quantiles as above.
#' @param legend Logical of length 1: draw the legend below the panel when
#'   several settings are plotted. Defaults to `TRUE`.
#' @param ... Unused; an argument passed here signals an error.
#'
#' @return `autoplot()` returns a ggplot2 object; `plot()` prints the same
#'   object and returns it invisibly.
#'
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' # The graph has one curve per value of the scale, at a shared mean.
#' plot(gaussian1_distrib(), list(mu = 0, sigma = c(1, 2, 4)))
#'
#' # A heavy-tailed family needs its own window: the 99.5% quantile of a
#' # standard Cauchy is about ten times as far out as the window drawn here.
#' plot(cauchy_distrib(), list(mu = 0, sigma = 1), xlim = c(-6, 6))
#'
#' @seealso [plot.discrete_distrib()] for a discrete family;
#'   [plot.distrib_fit()] to compare a fit with data;
#'   [plot_settings()], [plot_keys()] and [plot_labels()] for the pieces.
#' @keywords internal
S7::method(plot, continuous_distrib) <- function(x, theta, xlim = NULL,
                                                 legend = TRUE, ...) {
  p <- autoplot_continuous_distrib(x, theta, xlim = xlim, legend = legend, ...)
  print(p)
  invisible(p)
}

#' @rdname plot.continuous_distrib
autoplot_continuous_distrib <- function(object, theta, xlim = NULL,
                                        legend = TRUE, ...) {
  check_plot_args(legend, ...)
  x <- object
  st <- plot_univariate_setup(x, theta)
  ps <- st$ps
  labs <- st$labs

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
  need_ggplot()

  seq_x <- seq(xlim[1], xlim[2], length.out = 1000)
  dens <- lapply(ps$settings,
                 function(th) distrib_pdf(x, seq_x, th, log = FALSE))

  keys <- plot_keys(ps$k)
  key <- if (is.null(labs$legend)) "density" else labs$legend
  d <- data.frame(y = rep(seq_x, ps$k), dens = unlist(dens),
                  setting = factor(rep(key, each = length(seq_x)),
                                   levels = key))
  # an infinite density at an edge of the support is left out of the curve
  d <- d[is.finite(d$dens), , drop = FALSE]
  .data <- ggplot2::.data
  ggplot2::ggplot(d, ggplot2::aes(x = .data$y, y = .data$dens,
                                  colour = .data$setting,
                                  linetype = .data$setting)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::scale_colour_manual(name = NULL,
                                 values = stats::setNames(keys$col, key)) +
    ggplot2::scale_linetype_manual(name = NULL,
                                   values = stats::setNames(keys$lty, key)) +
    ggplot2::labs(title = labs$main, x = "y", y = "density") +
    ggplot2::theme_bw() +
    plot_key_theme(legend && ps$k > 1L)
}

S7::method(autoplot, continuous_distrib) <- autoplot_continuous_distrib

#' @title Plot Method for Discrete Distributions
#'
#' @name plot.discrete_distrib
#'
#' @description
#' `autoplot()` returns a graph of the probability mass of a discrete family
#' as stems over the integers in range, one set per parameter setting, and
#' `plot()` prints that graph. A component of `theta` given as a vector gives
#' one setting per element, exactly as in [plot.continuous_distrib()].
#'
#' @details
#' # The two keys, and why they are not the continuous ones
#' Settings are separated by color and symbol together, and the symbol sits
#' at the top of each stem, where the point already is. The line type is not
#' used here, because it would dash the stems and a dashed stem reads as a
#' broken one, so at a support of any size the panel fills with fragments that
#' cross each other. The parameters that vary are named in a legend below the
#' panel and those held fixed are stated once in the title.
#'
#' Several settings are also shifted sideways by a fraction of the spacing, so
#' that equal masses at one support point stay countable instead of standing
#' one in front of another.
#'
#' The graph is drawn with \pkg{ggplot2}, which is suggested; the methods
#' signal an error that names it when it is not installed.
#'
#' @param x,object An object inheriting from `discrete_distrib`.
#' @param theta A named list or numeric vector holding one entry per parameter
#'   of `x`. A component of length \eqn{k > 1} draws \eqn{k} sets of stems, and
#'   every component must have length 1 or that same \eqn{k}. When `theta` is
#'   missing, random values for all parameters are drawn by
#'   [generate_random_theta()].
#' @param xlim A numeric vector of length 2, or `NULL` (the default) to compute
#'   the range from the quantiles over every setting.
#' @param legend Logical of length 1: draw the legend below the panel when
#'   several settings are plotted. Defaults to `TRUE`.
#' @param ... Unused; an argument passed here signals an error.
#'
#' @return `autoplot()` returns a ggplot2 object; `plot()` prints the same
#'   object and returns it invisibly.
#'
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' # The three rates have masses that are countable at every support point
#' # because the settings are offset sideways.
#' plot(poisson_distrib(), list(mu = c(1, 4, 10)))
#'
#' # Overdispersion at a fixed mean: the mass spreads as theta falls.
#' plot(negbin2_distrib(), list(mu = 5, theta = c(0.5, 2, 50)))
#'
#' @seealso [plot.continuous_distrib()] for a density;
#'   [plot.distrib_fit()] to compare a fit with data;
#'   [plot_settings()], [plot_keys()] and [plot_labels()] for the pieces.
#' @keywords internal
S7::method(plot, discrete_distrib) <- function(x, theta, xlim = NULL,
                                               legend = TRUE, ...) {
  p <- autoplot_discrete_distrib(x, theta, xlim = xlim, legend = legend, ...)
  print(p)
  invisible(p)
}

#' @rdname plot.discrete_distrib
autoplot_discrete_distrib <- function(object, theta, xlim = NULL,
                                      legend = TRUE, ...) {
  check_plot_args(legend, ...)
  x <- object
  st <- plot_univariate_setup(x, theta)
  ps <- st$ps
  labs <- st$labs

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
  need_ggplot()

  # Generate strict integer sequence clamped to mathematical bounds
  seq_x <- seq(floor(xlim[1]), ceiling(xlim[2]), by = 1)
  seq_x <- seq_x[seq_x >= x@bounds[1] & seq_x <= x@bounds[2]]

  dens <- lapply(ps$settings,
                 function(th) distrib_pdf(x, seq_x, th, log = FALSE))

  # A shift of a fraction of the unit spacing, so that equal masses at one
  # support point stay countable. The stems are drawn solid whatever the
  # setting and told apart by their symbol: a dashed stem reads as a broken
  # one, and at a support of any size the panel fills with fragments.
  shift <- if (ps$k > 1L) (seq_len(ps$k) - (ps$k + 1) / 2) * (0.55 / ps$k)
           else 0
  keys <- plot_keys(ps$k)
  key <- if (is.null(labs$legend)) "probability" else labs$legend
  n <- length(seq_x)
  d <- data.frame(at = as.vector(outer(seq_x, shift, `+`)),
                  prob = unlist(dens),
                  setting = factor(rep(key, each = n), levels = key))
  .data <- ggplot2::.data
  ggplot2::ggplot(d, ggplot2::aes(x = .data$at, y = .data$prob,
                                  colour = .data$setting)) +
    ggplot2::geom_segment(ggplot2::aes(xend = .data$at, yend = 0),
                          linewidth = if (ps$k > 1L) 0.5 else 0.7) +
    ggplot2::geom_point(ggplot2::aes(shape = .data$setting),
                        size = if (ps$k > 1L) 1.6 else 2) +
    ggplot2::scale_colour_manual(name = NULL,
                                 values = stats::setNames(keys$col, key)) +
    ggplot2::scale_shape_manual(name = NULL,
                                values = stats::setNames(keys$pch, key)) +
    ggplot2::scale_x_continuous(breaks = integer_breaks) +
    ggplot2::labs(title = labs$main, x = "y", y = "probability") +
    ggplot2::theme_bw() +
    plot_key_theme(legend && ps$k > 1L)
}

S7::method(autoplot, discrete_distrib) <- autoplot_discrete_distrib
