#' @include distrib.R generics.R

#' Print Method for `distrib` Objects
#'
#' @description Custom S7 print method for objects inheriting from `distrib`.
#' @param x An object inheriting from class \code{"distrib"}.
#' @param ... Additional arguments (currently unused).
#' @return The object \code{x} invisibly.
#' @export
S7::method(print, distrib) <- function(x, ...) {
  # Format the distribution name (Capitalize first letters)
  d_name <- gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", x@distrib_name, perl = TRUE)
  
  # Determine type based on S7 inheritance
  d_type <- if (S7::S7_inherits(x, continuous_distrib)) {
    "Continuous"
  } else if (S7::S7_inherits(x, discrete_distrib)) {
    "Discrete"
  } else {
    "Unknown"
  }
  
  cat(sprintf("Distribution: %s\n", d_name))
  cat(sprintf("Type:         %s\n", d_type))
  cat(sprintf("Dimensions:   %s\n", x@dimension))
  
  cat("\nParameters:\n")
  
  max_param_len <- max(nchar(x@params), na.rm = TRUE)
  
  for (i in seq_len(x@n_params)) {
    p_name <- x@params[i]
    
    interpretation <- if (!is.null(x@params_interpretation)) {
      paste0("(", x@params_interpretation[[p_name]], ")")
    } else {
      ""
    }
    
    link_obj <- x@link_params[[p_name]]
    link_name <- if (is.null(link_obj)) "none" else link_obj$link_name
    
    bounds <- x@params_bounds[[p_name]]
    domain_str <- paste0("(", bounds[1], ", ", bounds[2], ")")
    
    cat(sprintf(
      "  %-*s %-20s | Link: %-10s | Domain: %s\n",
      max_param_len, p_name,
      interpretation,
      link_name,
      domain_str
    ))
  }
  
  invisible(x)
}

#' Generate Random Parameters for `distrib` Objects
#'
#' @description 
#' Generates a named list of sensible random parameters based on the mathematical domain (\code{params_bounds}) of the given distribution.
#' 
#' @param distrib An object inheriting from class \code{"distrib"}.
#' @param ... Additional arguments (currently unused).
#' @return A named list of randomly generated parameters.
#' @export
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

#' Plot Method for Continuous Distributions
#'
#' @description
#' Visualizes the Probability Density Function (PDF) of a continuous distribution object.
#'
#' @param x An object of class \code{"continuous_distrib"}.
#' @param theta A named list or vector of parameters matching \code{x@params}.
#' @param xlim Optional numeric vector of length 2 indicating the x-axis range.
#' @param ... Additional arguments passed to the base \code{\link{plot}} function.
#' @importFrom graphics plot
#' @export
S7::method(plot, continuous_distrib) <- function(x, theta, xlim = NULL, ...) {
  if (missing(theta)) {
    theta <- generate_random_theta(x)
    param_str <- paste(names(theta), round(unlist(theta), 3), sep = " = ", collapse = ", ")
    message(sprintf("Argument 'theta' is missing. Automatically generating random parameters: %s", param_str))
  }
  if (is.numeric(theta) && !is.list(theta)) theta <- as.list(theta)
  if (!all(x@params %in% names(theta))) {
    stop(sprintf("Missing parameters in 'theta'. Expected: %s", paste(x@params, collapse = ", ")))
  }
  
  if (is.null(xlim)) {
    # Evaluate quantiles to find the meaningful mass, suppressing warnings for extreme tails
    lower_q <- suppressWarnings(distrib_quantile(x, 0.005, theta))
    upper_q <- suppressWarnings(distrib_quantile(x, 0.995, theta))
    
    span <- max(upper_q - lower_q, 1.0, na.rm = TRUE)
    
    # Clamp limits strictly to mathematical bounds to prevent evaluating out-of-domain points
    xlim <- c(
      max(lower_q - 0.1 * span, x@bounds[1], na.rm = TRUE),
      min(upper_q + 0.1 * span, x@bounds[2], na.rm = TRUE)
    )
  }
  
  seq_x <- seq(xlim[1], xlim[2], length.out = 1000)
  dens_y <- distrib_pdf(x, seq_x, theta, log = FALSE)
  
  dots <- list(...)
  if (is.null(dots$main)) {
    param_str <- paste(names(theta), round(unlist(theta), 2), sep = " = ", collapse = ", ")
    dots$main <- paste0(x@distrib_name, " distribution\n(", param_str, ")")
  }
  if (is.null(dots$xlab)) dots$xlab <- "y"
  if (is.null(dots$ylab)) dots$ylab <- "Density (PDF)"
  if (is.null(dots$col)) dots$col <- "black"
  if (is.null(dots$lwd)) dots$lwd <- 2
  
  plot_args <- c(list(x = seq_x, y = dens_y, xlim = xlim, type = "l"), dots)
  do.call(graphics::plot, plot_args)
  
  invisible(x)
}

#' Plot Method for Discrete Distributions
#'
#' @description
#' Visualizes the Probability Mass Function (PMF) of a discrete distribution object.
#'
#' @param x An object of class \code{"discrete_distrib"}.
#' @param theta A named list or vector of parameters matching \code{x@params}.
#' @param xlim Optional numeric vector of length 2 indicating the x-axis range.
#' @param ... Additional arguments passed to the base \code{\link{plot}} function.
#' @importFrom graphics plot segments points
#' @export
S7::method(plot, discrete_distrib) <- function(x, theta, xlim = NULL, ...) {
  if (missing(theta)) {
    theta <- generate_random_theta(x)
    param_str <- paste(names(theta), round(unlist(theta), 3), sep = " = ", collapse = ", ")
    message(sprintf("Argument 'theta' is missing. Automatically generating random parameters: %s", param_str))
  }
  if (is.numeric(theta) && !is.list(theta)) theta <- as.list(theta)
  if (!all(x@params %in% names(theta))) {
    stop(sprintf("Missing parameters in 'theta'. Expected: %s", paste(x@params, collapse = ", ")))
  }
  
  if (is.null(xlim)) {
    lower_q <- suppressWarnings(distrib_quantile(x, 0.005, theta))
    upper_q <- suppressWarnings(distrib_quantile(x, 0.995, theta))
    
    xlim <- c(
      max(floor(lower_q), x@bounds[1], na.rm = TRUE),
      min(ceiling(upper_q), x@bounds[2], na.rm = TRUE)
    )
  }
  
  # Generate strict integer sequence clamped to mathematical bounds
  seq_x <- seq(floor(xlim[1]), ceiling(xlim[2]), by = 1)
  seq_x <- seq_x[seq_x >= x@bounds[1] & seq_x <= x@bounds[2]]
  
  dens_y <- distrib_pdf(x, seq_x, theta, log = FALSE)
  
  dots <- list(...)
  if (is.null(dots$main)) {
    param_str <- paste(names(theta), round(unlist(theta), 2), sep = " = ", collapse = ", ")
    dots$main <- paste0(x@distrib_name, " distribution\n(", param_str, ")")
  }
  if (is.null(dots$xlab)) dots$xlab <- "y"
  if (is.null(dots$ylab)) dots$ylab <- "Probability (PMF)"
  if (is.null(dots$col)) dots$col <- "black"
  if (is.null(dots$lwd)) dots$lwd <- 2
  
  plot_args <- c(list(x = seq_x, y = dens_y, xlim = xlim, type = "n"), dots)
  do.call(graphics::plot, plot_args)
  
  graphics::segments(x0 = seq_x, y0 = 0, x1 = seq_x, y1 = dens_y, col = dots$col, lwd = dots$lwd)
  graphics::points(seq_x, dens_y, pch = 16, col = dots$col)
  
  invisible(x)
}