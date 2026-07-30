#' @include distrib.R generics.R utility_functions.R link_scale.R
NULL

# --- internal helpers ------------------------------------------------------

#' Map Parameters to the Link Scale
#'
#' @description
#' Applies each parameter's link, \eqn{\eta_i = g_i(\theta_i)}.
#'
#' @details
#' Optimization is carried out on \eqn{\eta}, which is unconstrained, so this is
#' how a starting value expressed in natural parameters enters the optimizer.
#' Only the first element of each parameter is taken: a fit estimates one
#' \eqn{\theta} for the whole sample.
#'
#' @param distrib An object inheriting from class \code{"distrib"}.
#' @param theta A named list of parameters on the natural scale.
#'
#' @return A numeric vector of length \code{length(distrib@params)}.
#'
#' @seealso \code{\link{fit_theta_from_eta}}, the inverse.
#' @keywords internal
fit_eta_from_theta <- function(distrib, theta) {
  params <- distrib@params
  vapply(seq_along(params), function(i) {
    linkfunctions7::linkfun(distrib@link_params[[params[i]]], theta[[i]])[1]
  }, numeric(1))
}

#' Map the Link Scale Back to Parameters
#'
#' @description
#' Applies each parameter's inverse link, \eqn{\theta_i = g_i^{-1}(\eta_i)}.
#'
#' @details
#' The inverse of \code{\link{fit_eta_from_theta}}. Because every link maps onto
#' its parameter's domain by construction, a \eqn{\theta} obtained this way is
#' admissible whatever the optimizer proposed -- which is the reason for working
#' on the link scale at all, and the reason a confidence interval built there and
#' mapped back cannot run outside the domain.
#'
#' @param distrib An object inheriting from class \code{"distrib"}.
#' @param eta A numeric vector of linear predictors, one per parameter.
#'
#' @return A named list of parameters on the natural scale.
#'
#' @seealso \code{\link{fit_eta_from_theta}}
#' @keywords internal
fit_theta_from_eta <- function(distrib, eta) {
  params <- distrib@params
  out <- lapply(seq_along(params), function(i) {
    linkfunctions7::linkinv(distrib@link_params[[params[i]]], eta[i])
  })
  names(out) <- params
  out
}

#' Jacobian of the Inverse Link at the Estimate
#'
#' @description
#' The first derivative \eqn{h_i'(\eta_i)} of each parameter's inverse link, one
#' entry per parameter.
#'
#' @details
#' This is the diagonal Jacobian the delta method needs to carry a standard error
#' from the link scale, where it is computed, to the parameter scale, where it is
#' reported.
#'
#' @param distrib An object inheriting from class \code{"distrib"}.
#' @param eta A numeric vector of linear predictors, one per parameter.
#'
#' @return A numeric vector of length \code{length(distrib@params)}.
#'
#' @seealso \code{\link{fit_distrib}}
#' @keywords internal
fit_dtheta_deta <- function(distrib, eta) {
  params <- distrib@params
  vapply(seq_along(params), function(i) {
    linkfunctions7::dlinkinv(distrib@link_params[[params[i]]], eta[i])[1]
  }, numeric(1))
}

#' Summed Score on the Link Scale
#'
#' @description
#' The gradient of the total log-likelihood with respect to \eqn{\eta}, summed
#' over observations.
#'
#' @param distrib An object inheriting from class \code{"distrib"}.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters on the natural scale.
#'
#' @return A numeric vector of length \code{length(distrib@params)}.
#'
#' @seealso \code{\link{fit_hess_matrix}}
#' @keywords internal
fit_score <- function(distrib, y, theta) {
  g <- distrib_gradient(distrib, y, theta, scale = "link")
  vapply(g, function(v) sum(v), numeric(1))
}

#' Summed Hessian on the Link Scale, as a Matrix
#'
#' @description
#' Assembles the package's named list of Hessian components into the symmetric
#' \eqn{p \times p} matrix an optimizer wants, summed over observations.
#'
#' @details
#' Derivative components are stored one per unique index pair, since the Hessian
#' is symmetric; this fills both triangles. Passing \code{expected = TRUE} gives
#' the expected Hessian, which is what makes Fisher scoring possible -- and what
#' allows a fit on a non-regular family such as the Laplace, where the observed
#' Hessian is degenerate but the information is not.
#'
#' @param distrib An object inheriting from class \code{"distrib"}.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters on the natural scale.
#' @param expected Logical; whether to use the expected Hessian.
#'
#' @return A symmetric numeric matrix with dimnames taken from the parameters.
#'
#' @seealso \code{\link{fit_score}}, \code{\link{fit_distrib}}
#' @keywords internal
fit_hess_matrix <- function(distrib, y, theta, expected) {
  params <- distrib@params
  p <- length(params)
  h <- if (expected) {
    distrib_expected_hessian(distrib, y, theta, scale = "link")
  } else {
    distrib_hessian(distrib, y, theta, scale = "link")
  }
  M <- matrix(0, p, p, dimnames = list(params, params))
  for (i in seq_len(p)) {
    for (j in i:p) {
      v <- sum(h[[paste(params[c(i, j)], collapse = "_")]])
      M[i, j] <- v
      M[j, i] <- v
    }
  }
  M
}

#' Total Log-Likelihood
#'
#' @description
#' The sum of the log-density over the observations, the objective the fit
#' maximises.
#'
#' @param distrib An object inheriting from class \code{"distrib"}.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters on the natural scale.
#'
#' @return A single number.
#'
#' @seealso \code{\link{fit_distrib}}
#' @keywords internal
fit_loglik <- function(distrib, y, theta) {
  sum(distrib_pdf(distrib, y, theta, log = TRUE))
}

#' @title S7 Class for Maximum-Likelihood Fits
#' @name distrib_fit_class
#' @description
#' Object returned by \code{\link{fit_distrib}}, holding the estimates on both the
#' parameter and the link scale together with their uncertainty.
#' @param distrib The fitted \code{distrib} object.
#' @param y The observations the fit was computed from, kept so that the fitted
#'   distribution can be compared with the data (see \code{\link{plot.distrib_fit}}).
#' @param n Number of observations.
#' @param coefficients Named estimates on the parameter scale.
#' @param se Standard errors on the parameter scale (delta method).
#' @param ci Matrix of confidence limits on the parameter scale.
#' @param eta Estimates on the link scale.
#' @param se_eta Standard errors on the link scale.
#' @param ci_eta Matrix of confidence limits on the link scale.
#' @param vcov Variance-covariance matrix on the parameter scale.
#' @param vcov_eta Variance-covariance matrix on the link scale.
#' @param loglik Maximised log-likelihood.
#' @param aic,bic Information criteria.
#' @param iterations Number of iterations used.
#' @param converged Logical convergence flag.
#' @param method Optimisation method actually used.
#' @param level Confidence level.
#' @section Methods:
#' Methods implemented for this class:
#'   \code{\link[=coef.distrib_fit]{coef()}},
#'   \code{\link[=logLik.distrib_fit]{logLik()}},
#'   \code{\link[=plot.distrib_fit]{plot()}},
#'   \code{\link[=print.distrib_fit]{print()}},
#'   \code{\link[=simulate.distrib_fit]{simulate()}},
#'   \code{\link[=vcov.distrib_fit]{vcov()}}
#'
#' @export
distrib_fit <- S7::new_class("distrib_fit",
  properties = list(
    distrib      = distrib,
    y            = S7::class_numeric,
    n            = S7::class_numeric,
    coefficients = S7::class_numeric,
    se           = S7::class_numeric,
    ci           = S7::class_any,
    eta          = S7::class_numeric,
    se_eta       = S7::class_numeric,
    ci_eta       = S7::class_any,
    vcov         = S7::class_any,
    vcov_eta     = S7::class_any,
    loglik       = S7::class_numeric,
    aic          = S7::class_numeric,
    bic          = S7::class_numeric,
    iterations   = S7::class_numeric,
    converged    = S7::class_logical,
    method       = S7::class_character,
    level        = S7::class_numeric
  )
)

#' Maximum-Likelihood Estimation
#'
#' @description
#' Fits a distribution to an i.i.d. sample by maximum likelihood. The optimisation
#' is carried out on the \strong{link (real) scale}, where the parameters are
#' unconstrained, using the analytical score and information supplied by the
#' distribution (\code{scale = "link"}). Estimates are then mapped back to the
#' parameter scale and reported with standard errors and confidence intervals.
#'
#' @param distrib An object inheriting from class \code{"distrib"}.
#' @param y A numeric vector of observations.
#' @param start Optional named list of starting values \strong{on the parameter
#'   scale}. If \code{NULL} (default) starting values are drawn with
#'   \code{\link{generate_random_theta}}, with a few random restarts on failure;
#'   supplying a sensible \code{start} makes convergence faster and more reliable.
#' @param method Optimisation method. \code{"fisher"} (default) uses Fisher
#'   scoring with the expected information, \code{"newton"} uses the observed
#'   Hessian, and \code{"bfgs"} uses \code{\link[stats]{optim}} with the analytical
#'   gradient. Fisher scoring and Newton fall back to BFGS if they fail to converge.
#' @param maxit Maximum number of iterations. Defaults to 200.
#' @param tol Convergence tolerance on the score and on the log-likelihood
#'   increment. Defaults to \code{1e-10}.
#' @param level Confidence level for the intervals. Defaults to 0.95.
#' @param n_start Number of random restarts attempted when \code{start} is
#'   \code{NULL} and the first attempt fails. Defaults to 5.
#'
#' @return An object of class \code{\link{distrib_fit}}; see its documentation for
#'   the available components. \code{coef()}, \code{vcov()} and \code{logLik()}
#'   methods are provided.
#'
#' @details
#' \strong{Why the link scale.} Optimising \eqn{\eta \in \mathbb{R}^p} rather than
#' the constrained \eqn{\theta} removes the need for box constraints: a variance
#' can never become negative, a probability never leaves \eqn{(0,1)}. The score and
#' information on that scale are obtained exactly (not numerically) through the
#' chain rule described in \code{\link{link_scale_derivatives}}.
#'
#' \strong{Standard errors.} The variance-covariance matrix on the link scale is
#' the inverse of the information at the optimum. It is mapped to the parameter
#' scale by the delta method,
#' \deqn{\widehat{\mathrm{Var}}(\hat\theta) = J\,\widehat{\mathrm{Var}}(\hat\eta)\,J, \qquad
#'       J = \mathrm{diag}\!\left(\frac{d g^{-1}}{d\eta}\Big|_{\hat\eta}\right)}
#'
#' \strong{Confidence intervals.} Intervals are built symmetrically on the link
#' scale, \eqn{\hat\eta \pm z_{1-\alpha/2}\,\mathrm{se}(\hat\eta)}, and then mapped
#' through \eqn{g^{-1}}. This guarantees that the reported limits always respect the
#' parameter's domain (a variance interval cannot contain negative values), which a
#' symmetric interval on the parameter scale would not.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' d <- gaussian_distrib()
#' y <- distrib_rng(d, 500, list(mu = 2, sigma = 3))
#' fit <- fit_distrib(d, y)
#' fit
#' coef(fit)
#' vcov(fit)
#'
#' # a bounded parameter: the interval never leaves (0, 1)
#' b <- bernoulli_distrib()
#' fit_distrib(b, rbinom(50, 1, 0.9))
#' }
#'
#' @seealso \code{\link{link_scale_derivatives}}, \code{\link{check_distrib}}
#' @importFrom stats optim qnorm setNames
#' @export
fit_distrib <- function(distrib, y, start = NULL,
                        method = c("fisher", "newton", "bfgs"),
                        maxit = 200, tol = 1e-10, level = 0.95, n_start = 5) {
  method <- match.arg(method)
  params <- distrib@params
  p <- length(params)
  n <- length(y)

  nll <- function(eta) {
    th <- fit_theta_from_eta(distrib, eta)
    # Trial points are probed all over the link scale, including places where a
    # parameter overflows and the density warns ("NaNs produced"). Those
    # warnings say nothing about the fit -- the point is simply rejected by
    # returning Inf -- so they are not passed on to the user.
    v <- suppressWarnings(tryCatch(-fit_loglik(distrib, y, th), error = function(e) Inf))
    if (!is.finite(v)) Inf else v
  }

  # --- starting values ----------------------------------------------------
  starts <- if (!is.null(start)) {
    list(fit_eta_from_theta(distrib, align_theta(distrib, start)))
  } else {
    lapply(seq_len(max(1L, n_start)), function(i) {
      fit_eta_from_theta(distrib, generate_random_theta(distrib))
    })
  }

  # --- optimisation -------------------------------------------------------
  run_bfgs <- function(eta0) {
    gr <- function(eta) {
      th <- fit_theta_from_eta(distrib, eta)
      -fit_score(distrib, y, th)
    }
    op <- stats::optim(eta0, nll, gr, method = "BFGS",
                       control = list(maxit = maxit, reltol = 1e-12))
    list(eta = op$par, converged = op$convergence == 0,
         iterations = if (is.null(op$counts[[1]])) NA_real_ else op$counts[[1]],
         method = "BFGS")
  }

  run_scoring <- function(eta0, expected) {
    eta <- eta0
    ll <- -nll(eta)
    if (!is.finite(ll)) return(NULL)
    conv <- FALSE
    it <- 0L

    for (it in seq_len(maxit)) {
      th <- fit_theta_from_eta(distrib, eta)
      U <- fit_score(distrib, y, th)
      H <- fit_hess_matrix(distrib, y, th, expected = expected)
      I <- -H                                       # information (should be PD)

      step <- tryCatch(solve(I, U), error = function(e) NULL)
      if (is.null(step) || any(!is.finite(step))) return(NULL)

      # step halving on the log-likelihood
      accepted <- FALSE
      s <- 1
      for (k in 1:30) {
        cand <- eta + s * step
        ll_new <- -nll(cand)
        if (is.finite(ll_new) && ll_new >= ll - 1e-12) {
          accepted <- TRUE
          break
        }
        s <- s / 2
      }
      if (!accepted) return(NULL)

      delta_ll <- ll_new - ll
      eta <- cand
      ll <- ll_new

      if (max(abs(U)) < tol || abs(delta_ll) < tol * (abs(ll) + tol)) {
        conv <- TRUE
        break
      }
    }

    list(eta = eta, converged = conv, iterations = it,
         method = if (expected) "Fisher scoring" else "Newton-Raphson")
  }

  res <- NULL
  for (eta0 in starts) {
    if (!is.finite(nll(eta0))) next
    # An error raised inside the optimiser must be treated like a failure to
    # converge, not propagated: at an awkward parameter value the quadrature
    # behind a numerically-approximated expected Hessian can fail outright
    # ("the integral is probably divergent"), and without this the random
    # restarts and the BFGS fallback promised below never get their turn.
    res <- tryCatch(
      switch(method,
        fisher = run_scoring(eta0, expected = TRUE),
        newton = run_scoring(eta0, expected = FALSE),
        bfgs   = run_bfgs(eta0)
      ),
      error = function(e) NULL
    )
    if (is.null(res) || !isTRUE(res$converged)) {
      alt <- tryCatch(run_bfgs(eta0), error = function(e) NULL)
      if (!is.null(alt) && isTRUE(alt$converged)) res <- alt
    }
    if (!is.null(res) && isTRUE(res$converged)) break
  }

  if (is.null(res)) {
    stop("Optimisation failed from every starting value; supply 'start'.", call. = FALSE)
  }

  # --- inference at the optimum -------------------------------------------
  eta_hat <- res$eta
  theta_hat <- fit_theta_from_eta(distrib, eta_hat)

  # The estimates stand on their own; if the information cannot be evaluated at
  # the optimum the fit is still returned, with a missing variance matrix rather
  # than an error that throws the estimates away too.
  I_eta <- tryCatch(
    -fit_hess_matrix(distrib, y, theta_hat, expected = !identical(method, "newton")),
    error = function(e) matrix(NA_real_, p, p, dimnames = list(params, params))
  )
  V_eta <- tryCatch(solve(I_eta), error = function(e) matrix(NA_real_, p, p))
  dimnames(V_eta) <- list(params, params)

  J <- fit_dtheta_deta(distrib, eta_hat)
  V_theta <- diag(J, p, p) %*% V_eta %*% diag(J, p, p)
  dimnames(V_theta) <- list(params, params)

  se_eta <- sqrt(pmax(diag(V_eta), 0))
  se_theta <- sqrt(pmax(diag(V_theta), 0))

  z <- stats::qnorm(1 - (1 - level) / 2)
  ci_eta <- cbind(lower = eta_hat - z * se_eta, upper = eta_hat + z * se_eta)
  rownames(ci_eta) <- params

  # map the link-scale interval through g^{-1}; sort in case the link decreases
  ci_theta <- t(vapply(seq_len(p), function(i) {
    lk <- distrib@link_params[[params[i]]]
    # na.last keeps the pair at length two: sort() drops NAs by default, and an
    # information matrix that could not be evaluated makes both ends missing.
    sort(c(linkfunctions7::linkinv(lk, ci_eta[i, 1]),
           linkfunctions7::linkinv(lk, ci_eta[i, 2])), na.last = TRUE)
  }, numeric(2)))
  dimnames(ci_theta) <- list(params, c("lower", "upper"))

  ll_hat <- fit_loglik(distrib, y, theta_hat)
  coefs <- stats::setNames(vapply(theta_hat, function(v) v[1], numeric(1)), params)

  distrib_fit(
    distrib = distrib, y = y, n = n,
    coefficients = coefs,
    se = stats::setNames(se_theta, params),
    ci = ci_theta,
    eta = stats::setNames(eta_hat, params),
    se_eta = stats::setNames(se_eta, params),
    ci_eta = ci_eta,
    vcov = V_theta, vcov_eta = V_eta,
    loglik = ll_hat,
    aic = -2 * ll_hat + 2 * p,
    bic = -2 * ll_hat + log(n) * p,
    iterations = res$iterations,
    converged = isTRUE(res$converged),
    method = res$method,
    level = level
  )
}

#' Print Method for Maximum-Likelihood Fits
#'
#' @name print.distrib_fit
#' @param x A \code{\link{distrib_fit}} object.
#' @param digits Number of significant digits. Defaults to 4.
#' @param ... Unused.
#' @return \code{x}, invisibly.
S7::method(print, distrib_fit) <- function(x, digits = 4, ...) {
  pct <- format(100 * x@level, trim = TRUE)
  lo <- paste0(format((1 - x@level) / 2 * 100, trim = TRUE), "%")
  hi <- paste0(format((1 + x@level) / 2 * 100, trim = TRUE), "%")

  cat("Maximum-likelihood fit: ", x@distrib@distrib_name, "\n", sep = "")
  cat("Observations: ", x@n,
      "   Log-likelihood: ", format(x@loglik, digits = digits),
      "   AIC: ", format(x@aic, digits = digits),
      "   BIC: ", format(x@bic, digits = digits), "\n", sep = "")
  cat("Method: ", x@method,
      if (x@converged) sprintf(" (converged in %s iterations)", x@iterations)
      else " (DID NOT CONVERGE)", "\n\n", sep = "")

  tab <- cbind(Estimate = x@coefficients, `Std. Error` = x@se, x@ci)
  colnames(tab) <- c("Estimate", "Std. Error", lo, hi)
  cat("Parameter scale (", pct, "% CI mapped from the link scale):\n", sep = "")
  print(round(tab, digits))

  links <- vapply(x@distrib@params,
                  function(p) x@distrib@link_params[[p]]@link_name, character(1))
  tab_eta <- cbind(Estimate = x@eta, `Std. Error` = x@se_eta)
  cat("\nLink scale (", paste(links, collapse = ", "), "):\n", sep = "")
  print(round(tab_eta, digits))

  invisible(x)
}

#' Extract Estimates from a Maximum-Likelihood Fit
#'
#' @name coef.distrib_fit
#' @param object A \code{\link{distrib_fit}} object.
#' @param scale Either \code{"parameter"} (default) or \code{"link"}.
#' @param ... Unused.
#' @return A named numeric vector of estimates.
#' @importFrom stats coef
S7::method(coef, distrib_fit) <- function(object, scale = c("parameter", "link"), ...) {
  scale <- match.arg(scale)
  if (scale == "link") object@eta else object@coefficients
}

#' Variance-Covariance Matrix of a Maximum-Likelihood Fit
#'
#' @name vcov.distrib_fit
#' @param object A \code{\link{distrib_fit}} object.
#' @param scale Either \code{"parameter"} (default) or \code{"link"}.
#' @param ... Unused.
#' @return A variance-covariance matrix.
#' @importFrom stats vcov
S7::method(vcov, distrib_fit) <- function(object, scale = c("parameter", "link"), ...) {
  scale <- match.arg(scale)
  if (scale == "link") object@vcov_eta else object@vcov
}

#' Log-Likelihood of a Maximum-Likelihood Fit
#'
#' @name logLik.distrib_fit
#' @param object A \code{\link{distrib_fit}} object.
#' @param ... Unused.
#' @return An object of class \code{logLik}.
#' @importFrom stats logLik
S7::method(logLik, distrib_fit) <- function(object, ...) {
  structure(object@loglik,
            df = length(object@coefficients),
            nobs = object@n,
            class = "logLik")
}

#' Simulate from a Fitted Distribution
#'
#' @name simulate.distrib_fit
#'
#' @description
#' Draws new samples from the fitted distribution, evaluated at the maximum
#' likelihood estimates. Each replicate has the same length as the data the model
#' was fitted to, which makes the result directly comparable with the observations
#' and suitable for a parametric bootstrap or a posterior-predictive style check.
#'
#' @param object A \code{\link{distrib_fit}} object.
#' @param nsim Number of replicates to draw. Defaults to 1.
#' @param seed Optional seed. If supplied it is used to initialise the generator,
#'   and the state of \code{.Random.seed} in effect before the call is restored
#'   afterwards, so that simulating does not disturb the calling stream. The seed
#'   actually used is attached to the result as the \code{"seed"} attribute.
#' @param ... Unused.
#'
#' @return A data frame with \code{nsim} columns named \code{sim_1}, ...,
#'   \code{sim_nsim}, each holding one replicate of \code{object@n} draws.
#'
#' @seealso \code{\link{fit_distrib}}, \code{\link{plot.distrib_fit}}
#'
#' @examples
#' set.seed(1)
#' y <- rnorm(200, 3, 2)
#' fit <- fit_distrib(gaussian_distrib(), y)
#'
#' sims <- simulate(fit, 20, seed = 42)
#' dim(sims)
#'
#' # a parametric bootstrap of any statistic
#' quantile(vapply(sims, median, numeric(1)), c(0.025, 0.975))
#'
#' @importFrom stats simulate
S7::method(simulate, distrib_fit) <- function(object, nsim = 1, seed = NULL, ...) {
  nsim <- as.integer(nsim)
  if (length(nsim) != 1L || is.na(nsim) || nsim < 1L) {
    stop("'nsim' must be a single positive integer.", call. = FALSE)
  }

  # The seed protocol of stats::simulate: honour `seed`, leave the caller's
  # random stream exactly as it was, and report what was used.
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    stats::runif(1)
  }
  saved <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (is.null(seed)) {
    used <- saved
  } else {
    set.seed(seed)
    used <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", saved, envir = .GlobalEnv), add = TRUE)
  }

  theta <- as.list(object@coefficients)
  out <- lapply(seq_len(nsim), function(i) distrib_rng(object@distrib, object@n, theta))
  names(out) <- paste0("sim_", seq_len(nsim))

  structure(as.data.frame(out), seed = used)
}

#' Plot a Fitted Distribution against the Data
#'
#' @name plot.distrib_fit
#'
#' @description
#' Compares the fitted distribution with the sample it was estimated from. For a
#' continuous distribution the observations are summarised by a kernel density
#' estimate, with the fitted density drawn on top and a rug of the data
#' underneath. For a discrete one the observed relative frequencies are drawn as
#' bars with the fitted probability mass overlaid, since a kernel density would
#' misrepresent a lattice-valued sample.
#'
#' @param x A \code{\link{distrib_fit}} object.
#' @param n_grid Number of points at which the fitted density is evaluated
#'   (continuous distributions only). Defaults to 512.
#' @param rug Logical; draw a rug of the observations. Defaults to \code{TRUE}
#'   when there are at most 2000 of them.
#' @param legend Logical; add a legend. Defaults to \code{TRUE}.
#' @param col_fit,col_data Colours of the fitted curve and of the empirical
#'   summary.
#' @param ... Further arguments passed to \code{\link[graphics]{plot}}, for
#'   instance \code{main}, \code{xlab} or \code{xlim}.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{fit_distrib}}, \code{\link{simulate.distrib_fit}}
#'
#' @examples
#' set.seed(1)
#' y <- rgamma(300, shape = 4, rate = 2)
#' fit <- fit_distrib(gamma_distrib(), y)
#' plot(fit)
#'
#' @importFrom graphics lines legend rug barplot points
S7::method(plot, distrib_fit) <- function(x, n_grid = 512, rug = NULL,
                                          legend = TRUE,
                                          col_fit = "#B22222", col_data = "#4682B4",
                                          ...) {
  y <- x@y
  if (!length(y)) {
    stop("This fit carries no data to plot.", call. = FALSE)
  }
  theta <- as.list(x@coefficients)
  d <- x@distrib
  dots <- list(...)
  if (is.null(dots$main)) {
    dots$main <- paste0("Fitted ", d@distrib_name, " and observed data")
  }
  if (is.null(dots$xlab)) dots$xlab <- "y"

  if (S7::S7_inherits(d, discrete_distrib)) {
    ks <- seq(max(min(y), d@bounds[1]), min(max(y), d@bounds[2]))
    obs <- as.numeric(table(factor(y, levels = ks))) / length(y)
    fitted <- distrib_pdf(d, ks, theta)

    if (is.null(dots$ylab)) dots$ylab <- "probability"
    mids <- do.call(graphics::barplot,
      c(list(obs, names.arg = ks, col = col_data, border = NA,
             ylim = c(0, max(obs, fitted) * 1.1)), dots))
    graphics::points(mids, fitted, col = col_fit, pch = 16)
    graphics::lines(mids, fitted, col = col_fit, lwd = 2)
    if (isTRUE(legend)) {
      graphics::legend("topright", bty = "n",
        legend = c("observed frequency", "fitted pmf"),
        fill = c(col_data, NA), border = c(NA, NA),
        lty = c(NA, 1), lwd = c(NA, 2), col = c(NA, col_fit))
    }
    return(invisible(x))
  }

  dens <- stats::density(y)
  # Evaluate the fit strictly inside the support: a density is often unbounded or
  # undefined exactly at the edge, and the kernel estimate spills across it.
  b <- d@bounds
  nudge <- function(v, into) if (is.finite(v)) v + into * 1e-8 * max(1, abs(v)) else v
  lo <- max(dens$x[1], nudge(b[1], 1))
  hi <- min(dens$x[length(dens$x)], nudge(b[2], -1))
  grid <- seq(lo, hi, length.out = n_grid)
  fitted <- distrib_pdf(d, grid, theta)
  fitted[!is.finite(fitted)] <- NA_real_

  if (is.null(dots$ylab)) dots$ylab <- "density"
  if (is.null(dots$ylim)) dots$ylim <- c(0, max(c(dens$y, fitted), na.rm = TRUE) * 1.05)
  if (is.null(dots$xlim)) dots$xlim <- range(dens$x)

  do.call(graphics::plot,
    c(list(dens$x, dens$y, type = "l", col = col_data, lwd = 2), dots))
  graphics::lines(grid, fitted, col = col_fit, lwd = 2)

  if (is.null(rug)) rug <- length(y) <= 2000
  if (isTRUE(rug)) graphics::rug(y, col = col_data)

  if (isTRUE(legend)) {
    graphics::legend("topright", bty = "n", lwd = 2,
      col = c(col_data, col_fit),
      legend = c("kernel density", paste("fitted", d@distrib_name)))
  }

  invisible(x)
}
