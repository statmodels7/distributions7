#' @include distrib.R generics.R
NULL

#' @title Print Method for `distrib` Objects
#' @name print.distrib
#' @description Custom S7 print method for objects inheriting from `distrib`.
#' @param x An object inheriting from class `"distrib"`.
#' @param ... Additional arguments (currently unused).
#' @return The object `x` invisibly.
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

#' Generate Random Parameters for `distrib` Objects
#'
#' @name generate_random_theta.distrib
#' @description
#' Generates a named list of sensible random parameters based on the mathematical domain (`params_bounds`) of the given distribution.
#' 
#' @param distrib An object inheriting from class `"distrib"`.
#' @param ... Additional arguments (currently unused).
#' @return A named list of randomly generated parameters.
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

#' The Parameter Settings a Plot Draws
#'
#' @description
#' Splits a `theta` whose components may be vectors into one setting per
#' curve, and reports which parameters vary across them.
#'
#' @details
#' A component of length \eqn{k > 1} asks for \eqn{k} curves and a component of
#' length one is held across all of them, so
#' `list(mu = 0, sigma = c(1, 2, 4))` is three settings sharing a mean.
#' Every component must therefore have length one or the same \eqn{k}; anything
#' else is rejected rather than recycled, since a length that divides \eqn{k}
#' is far more likely to be a mistake than a request.
#'
#' This meaning is available because a plot has no data to recycle against. The
#' derivative and density generics read a vector component as one value per
#' observation, which is a different question asked of the same object.
#'
#' @param x A distribution object.
#' @param theta A named list or numeric vector of parameters.
#'
#' @return A list with `settings` (one scalar `theta` per curve),
#'   `k` and `varying` (the names of the parameters that differ
#'   between settings).
#'
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

#' Colors, Line Types and Symbols for Overlaid Settings
#'
#' @description
#' The visual keys distinguishing several settings on one panel: a color, a
#' line type and a plotting symbol per setting, all cycled, so the settings are
#' told apart in color and in a printed copy that has none.
#'
#' @details
#' A continuous family is separated by color and line type, a discrete one by
#' color and symbol. Dashing a stem is what a line type would do there, and a
#' dashed stem reads as a broken one: at a support of any size the panel fills
#' with fragments that cross each other. The symbol carries the same
#' information without drawing anything extra, the point at the top of the stem
#' being already there.
#'
#' @param k The number of settings.
#' @param dots The caller's `...`; a `col`, `lty` or `pch`
#'   given there wins and is recycled over the settings.
#'
#' @return A list with `col`, `lty` and `pch`, each of length
#'   `k`.
#'
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

#' Where to Put the Key
#'
#' @description
#' The corner of the panel a legend goes in: the one the curves leave emptier.
#'
#' @details
#' A right-skewed family puts its mass on the left, so a key fixed at the top
#' right never meets it, while a left-skewed one is covered by exactly that
#' choice. The side is therefore read off the drawing: the density's center of
#' mass is compared with the middle of the horizontal range, and the key goes
#' to the other side.
#'
#' @param y The evaluation points.
#' @param dens A list of densities, one per setting.
#'
#' @return `"topright"` or `"topleft"`.
#'
#' @keywords internal
plot_legend_side <- function(y, dens) {
  w <- do.call(pmax, c(dens, list(na.rm = TRUE)))
  w[!is.finite(w) | w < 0] <- 0
  if (!any(w > 0)) return("topright")
  centre <- sum(y * w) / sum(w)
  if (centre < mean(range(y))) "topright" else "topleft"
}

#' Labels for a Plot of Several Settings
#'
#' @description
#' The legend entries and the title of a plot drawing several settings: the
#' parameters that vary go in the legend, one entry per curve, and those held
#' fixed go in the title, where they are stated once.
#'
#' @param x A distribution object.
#' @param ps The value of [plot_settings()].
#'
#' @return A list with `legend` (character, length `ps$k`, or
#'   `NULL` when nothing varies) and `main`.
#'
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

#' Plot Method for Continuous Distributions
#'
#' @name plot.continuous_distrib
#' @description
#' Visualizes the Probability Density Function (PDF) of a continuous distribution object.
#'
#' @details
#' A parameter given as a vector asks for one curve per element, so
#' `plot(d, list(mu = 0, sigma = c(1, 2, 4)))` draws three densities that
#' share a mean and differ in scale. The curves are separated by color and by
#' line type together, and the parameters that vary are named in a legend while
#' those held fixed are stated in the title. See [plot_settings()]
#' for the rule on lengths.
#'
#' The horizontal range covers every setting: it runs from the smallest
#' 0.5\% quantile to the largest 99.5\% quantile over them, widened by a tenth
#' and clamped to the support.
#'
#' @param x An object of class `"continuous_distrib"`.
#' @param theta A named list or vector of parameters matching `x@params`.
#'   A component of length \eqn{k > 1} draws \eqn{k} curves.
#' @param xlim Optional numeric vector of length 2 indicating the x-axis range.
#' @param legend Whether to draw the key when several settings are plotted.
#' @param ... Additional arguments passed to the base [plot()] function.
#'   `col` and `lty` given here win and are recycled over the curves.
#' @importFrom graphics plot
#' @keywords internal
#' @return The distribution, invisibly.
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

#' Plot Method for Discrete Distributions
#'
#' @name plot.discrete_distrib
#' @description
#' Visualizes the Probability Mass Function (PMF) of a discrete distribution object.
#'
#' @details
#' A parameter given as a vector asks for one setting per element, exactly as
#' in [plot.continuous_distrib()]. Several settings are drawn as
#' several sets of stems, shifted sideways so that one does not stand in front
#' of another, and separated by color and line type together.
#'
#' @param x An object of class `"discrete_distrib"`.
#' @param theta A named list or vector of parameters matching `x@params`.
#'   A component of length \eqn{k > 1} draws \eqn{k} sets of stems.
#' @param xlim Optional numeric vector of length 2 indicating the x-axis range.
#' @param legend Whether to draw the key when several settings are plotted.
#' @param ... Additional arguments passed to the base [plot()] function.
#'   `col` and `lty` given here win and are recycled over the settings.
#' @importFrom graphics plot segments points
#' @keywords internal
#' @return The distribution, invisibly.
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
