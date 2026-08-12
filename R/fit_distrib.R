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
#' @param approx How the expectation is approximated when the distribution has
#'   no closed form for it. Ignored when it has one, and when \code{expected}
#'   is \code{FALSE}.
#' @param nsim Monte Carlo sample size, used when \code{approx = "mc"}.
#'
#' @return A symmetric numeric matrix with dimnames taken from the parameters.
#'
#' @seealso \code{\link{fit_score}}, \code{\link{fit_distrib}}
#' @keywords internal
fit_hess_matrix <- function(distrib, y, theta, expected,
                            approx = "bartlett", nsim = 10000) {
  params <- distrib@params
  p <- length(params)
  h <- if (expected) {
    distrib_expected_hessian(distrib, y, theta,
      scale = "link", approx = approx, nsim = nsim
    )
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
#' maximizes.
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
#' @param loglik Maximized log-likelihood.
#' @param aic,bic Information criteria.
#' @param iterations Number of iterations used.
#' @param converged Logical convergence flag.
#' @param method Optimization method actually used.
#' @param criterion Which stopping rule ended the run, as optimizers7 reports it.
#' @param note Any remark the optimizer attached to the run.
#' @param counts How many times the objective and its gradient were evaluated.
#' @param score The max-norm of the score \strong{per observation} at the
#'   reported optimum, which is the quantity the stopping rule tested.
#' @param elapsed Seconds spent optimizing, summed over every starting value
#'   and every fallback attempted.
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
#' @return An object of class \code{distrib_fit}.
#'
#' @examples
#' set.seed(1)
#' y <- distrib_rng(gaussian1_distrib(), 200, list(mu = 1, sigma = 2))
#' fit <- fit_distrib(gaussian1_distrib(), y)
#' S7::S7_inherits(fit, distrib_fit)
#' coef(fit)
#' logLik(fit)
#'
#' @seealso \code{\link{fit_distrib}}
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
    criterion    = S7::class_character,
    note         = S7::class_character,
    counts       = S7::class_any,
    score        = S7::class_numeric,
    elapsed      = S7::class_numeric,
    level        = S7::class_numeric
  )
)

#' Maximum-Likelihood Estimation
#'
#' @description
#' Fits a distribution to an i.i.d. sample by maximum likelihood. The optimization
#' is carried out on the \strong{link (real) scale}, where the parameters are
#' unconstrained, using the analytical score and information supplied by the
#' distribution (\code{scale = "link"}). Estimates are then mapped back to the
#' parameter scale and reported with standard errors and confidence intervals.
#'
#' @param distrib An object inheriting from class \code{"distrib"}.
#' @param y A numeric vector of observations.
#' @param start Optional named list of starting values \strong{on the parameter
#'   scale}. If \code{NULL} (default) they come from
#'   \code{\link{distrib_start}}, which lets a family compute them from the
#'   data; families that do not say otherwise fall back to random draws, with
#'   restarts.
#' @param method How to optimize. One argument, taking one of three things:
#'   \itemize{
#'     \item \code{\link{fisher_scoring}()}, the default --- Newton's method
#'       with the \strong{expected} information in place of the Hessian, the
#'       object carrying how that information is to be obtained when the family
#'       has no closed form for it;
#'     \item an optimizer object from \pkg{optimizers7}, used as given and
#'       receiving the analytical gradient and the \strong{observed} Hessian,
#'       so that \code{method = lbfgs(criterion = crit_grad(1e-12))} selects
#'       both the algorithm and the stopping rule;
#'     \item one of the strings \code{"fisher"}, \code{"newton"} or
#'       \code{"bfgs"}, kept as short names for the three ready-made
#'       strategies. The first two fall back to BFGS if they fail to converge;
#'       an optimizer the caller chose is never silently replaced.
#'   }
#'   The iteration limit and the stopping rule belong to the method and are
#'   set there: on an optimizer object through its own \code{maxit} and
#'   \code{criterion}, on \code{\link{fisher_scoring}()} through the same two
#'   arguments, and otherwise left at the defaults of
#'   \code{\link[optimizers7]{crit_grad}} and of the optimizer. The objective
#'   handed to the optimizer is the \emph{mean} negative log-likelihood, so a
#'   tolerance on its gradient is a tolerance on the score \strong{per
#'   observation} whatever the sample size, and
#'   \code{\link[optimizers7]{crit_grad}} documents why its default sits where
#'   it does.
#' @param level Confidence level for the intervals. Defaults to 0.95.
#' @param n_start How many starting values to ask \code{\link{distrib_start}}
#'   for when \code{start} is \code{NULL}. Defaults to 5. A family that returns
#'   its own estimate returns one and ignores this.
#'
#' @return An object of class \code{\link{distrib_fit}}; see its documentation for
#'   the available components. \code{coef()}, \code{vcov()} and \code{logLik()}
#'   methods are provided.
#'
#' @details
#' \strong{Why the link scale.} Optimizing \eqn{\eta \in \mathbb{R}^p} rather than
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
#' d <- gaussian1_distrib()
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
#' @seealso \code{\link{link_scale_derivatives}}, \code{\link{check_distrib}},
#'   \code{\link[optimizers7]{minimize}}
#' @importFrom stats qnorm setNames
#' @export
fit_distrib <- function(distrib, y, start = NULL,
                        method = fisher_scoring(),
                        level = 0.95, n_start = 5) {
  # 'start' comes before 'method' in the signature, so an optimizer passed
  # positionally lands in it. What the caller then sees, several frames down,
  # is align_theta() refusing to coerce an S7 object to a list -- an error
  # that names neither the argument nor the mistake. The check is here, where
  # both are known.
  if (!is.null(start) && !is.list(start)) {
    hint <- if (S7::S7_inherits(start, optimizers7::optimizer) ||
                S7::S7_inherits(start, FisherScoring)) {
      "\n  It looks like an optimizer: pass it as 'method = ', since 'start'\n  is the third argument and 'method' the fourth."
    } else {
      ""
    }
    stop(sprintf(paste0(
      "'start' must be NULL or a named list of parameter values, and it is ",
      "'%s'.%s"
    ), paste(class(start), collapse = "/"), hint), call. = FALSE)
  }

  # One argument says how to optimize, and it takes one of three things: a
  # fisher_scoring() specification, an optimizers7 optimizer, or the name of
  # one of the three ready-made strategies. How the expected information is to
  # be approximated is a property of Fisher scoring and lives on that object,
  # not among fit_distrib()'s own arguments, where it would sit next to
  # optimizers that never look at it.
  optimizer <- NULL
  approx <- "bartlett"
  nsim <- 10000
  fs <- NULL

  if (S7::S7_inherits(method, FisherScoring)) {
    fs <- method
    method <- "fisher"
    approx <- fs@approx
    nsim <- fs@nsim
  } else if (S7::S7_inherits(method, optimizers7::optimizer)) {
    optimizer <- method
    method <- "custom"
    # A stopping rule the optimizer cannot evaluate is a mistake in the call,
    # not a numerical failure, and it is raised here so that it reaches the
    # caller with its own explanation rather than being caught by the restart
    # loop below and reported as a fit that never converged.
    optimizers7::check_criterion(optimizer)
  } else {
    method <- match.arg(method, c("fisher", "newton", "bfgs"))
  }

  # A strategy chosen where it will be ignored is a mistake in the call rather
  # than a harmless redundancy: silently accepting it is how a caller comes to
  # believe a fit used a method it did not.
  if (!is.null(fs) && !identical(fs@approx, "bartlett") &&
      has_exact_expected_hessian(distrib)) {
    stop(sprintf(paste0(
      "'%s' computes its expected information in closed form, so the 'approx'\n",
      "  of fisher_scoring() would be ignored. Use fisher_scoring() with no\n",
      "  arguments: the fit will take the exact expression."
    ), distrib@distrib_name), call. = FALSE)
  }

  params <- distrib@params
  p <- length(params)
  # The row count for a multivariate response, its length otherwise: n is the
  # number of OBSERVATIONS, which is what BIC and the printed summary mean.
  n <- n_obs(distrib, y)

  # The optimizer is handed the MEAN negative log-likelihood, and its gradient
  # and Hessian divided by n with it. Scaling an objective by a positive
  # constant leaves the maximum and every Newton step where they were --
  # H^{-1}g is unchanged when both are divided by n -- so nothing about the
  # path changes; what changes is what a stopping rule means. A criterion on
  # the gradient of the SUMMED log-likelihood tests a quantity whose attainable
  # floor grows with the sample, so the same fit met a bound of 1e-10 on one
  # platform and missed it on another. The score per observation is of order
  # one whatever n is, so a tolerance on it means the same thing everywhere,
  # and it means it for a criterion the CALLER supplies as much as for the
  # default -- which scaling the tolerance instead of the objective would not
  # have done.
  #
  # The scaling settles what a threshold MEANS but not how small one can be.
  # A line search accepts a step only when the objective decreases by a
  # definite amount, and near the maximum that decrease, of order
  # |U/n|^2 / (2 lambda), sinks below the rounding of the objective itself.
  # The search then rejects every step whatever the direction, so the smallest
  # gradient a run can reach is around sqrt(2 lambda eps |l/n|). Measured over
  # several families, methods and samples it is usually near 1e-15 and reaches
  # 1e-8, which is why the default tolerance is 1e-6 and not tighter.
  #
  # The reported quantities are the ordinary ones: loglik, the information and
  # the standard errors are all recomputed below from the unscaled likelihood.
  scale_n <- max(1, n)

  nll <- function(eta) {
    th <- fit_theta_from_eta(distrib, eta)
    # Trial points are probed all over the link scale, including places where a
    # parameter overflows and the density warns ("NaNs produced"). Those
    # warnings say nothing about the fit -- the point is simply rejected by
    # returning Inf -- so they are not passed on to the user.
    v <- suppressWarnings(tryCatch(-fit_loglik(distrib, y, th), error = function(e) Inf))
    if (!is.finite(v)) Inf else v / scale_n
  }

  # --- starting values ----------------------------------------------------
  # distrib_start() lets a family compute a start from the DATA. The default
  # draws at random from the parameter domains, as before, but a family that
  # knows its own estimator says so: the four-dimensional gaussian of the iris
  # measurements never converges from the origin and converges in one iteration
  # from the sample mean and covariance.
  starts <- if (!is.null(start)) {
    list(fit_eta_from_theta(distrib, align_theta(distrib, start)))
  } else {
    th0 <- distrib_start(distrib, y, n_start = max(1L, n_start))
    if (!is.list(th0) || !length(th0)) {
      stop("distrib_start() must return a non-empty list of parameter lists.",
           call. = FALSE)
    }
    lapply(th0, function(th) {
      fit_eta_from_theta(distrib, align_theta(distrib, th))
    })
  }

  # --- optimization -------------------------------------------------------
  # The objective, its gradient and its Hessian are those of the NEGATIVE
  # log-likelihood, since optimizers7 minimizes. Fisher scoring is Newton's
  # method with the expected information supplied in place of the observed
  # Hessian, so the two named strategies differ only in that argument.
  nll_gr <- function(eta) {
    th <- fit_theta_from_eta(distrib, eta)
    -fit_score(distrib, y, th) / scale_n
  }
  nll_he <- function(expected) {
    function(eta) {
      th <- fit_theta_from_eta(distrib, eta)
      -fit_hess_matrix(distrib, y, th,
        expected = expected, approx = approx, nsim = nsim
      ) / scale_n
    }
  }

  # How the run stops and how long it may take are properties of the method,
  # so they are read off the method and nowhere else. Whatever the caller did
  # not set stays at the optimizer's own default, which is the only place the
  # constant is written down: a copy of it here would keep its old value the
  # day the default moved. The objective the rule sees is the mean negative
  # log-likelihood, so a tolerance on its gradient is a tolerance on the score
  # per observation without the rule having to know n.
  opt_args <- list()
  if (!is.null(fs)) {
    if (!is.null(fs@criterion)) opt_args$criterion <- fs@criterion
    if (!is.null(fs@maxit)) opt_args$maxit <- fs@maxit
  }

  # Time is accumulated over every attempt rather than taken from the run that
  # is kept: what a caller wants to know is what the fit cost, and the restarts
  # and the fallback are part of that.
  spent <- 0

  run <- function(opt, eta0, he, label) {
    r <- optimizers7::minimize(opt, fn = nll, par = eta0, gr = nll_gr, he = he)
    if (length(r@elapsed) && is.finite(r@elapsed)) spent <<- spent + r@elapsed
    # The optimizer's gradient is that of -l(eta)/n, so its max-norm IS the
    # score per observation the stopping rule tested. Keeping it means the
    # object can say how close to stationary it ended, which is the one number
    # a run that did not converge is worth reading for.
    list(eta = r@par, converged = isTRUE(r@converged),
         iterations = r@iterations, method = label,
         value = r@value, criterion_met = r@criterion_met,
         message = r@message, counts = r@counts,
         score = if (length(r@gradient)) max(abs(r@gradient)) else NA_real_)
  }

  # The fallback inherits the rule and the budget of the method it stands in
  # for: a run reported under a stopping rule the caller did not ask for would
  # be a different fit under the same name.
  run_bfgs <- function(eta0) {
    run(do.call(optimizers7::bfgs, opt_args), eta0, NULL, "BFGS")
  }

  run_chosen <- function(eta0) {
    switch(method,
      fisher = run(do.call(optimizers7::newton, opt_args),
                   eta0, nll_he(TRUE), "Fisher scoring"),
      newton = run(do.call(optimizers7::newton, opt_args),
                   eta0, nll_he(FALSE), "Newton-Raphson"),
      bfgs   = run_bfgs(eta0),
      custom = run(optimizer, eta0, nll_he(FALSE), optimizer@name)
    )
  }

  # A discrete family's likelihood is a product of probabilities, so its
  # logarithm cannot be positive. That is a fact about the model rather than a
  # tolerance, and it is the one test that separates a run which has found a
  # better fit from one which has left the region where the mass function is
  # computable: the beta-binomial's two beta functions cancel to the last digit
  # once its shapes pass about 1e15, and every mass then comes back as one, so
  # the run reports a log-likelihood of zero and wins every comparison against
  # a real fit. The mass function is written to avoid that (see
  # betabinom_log_mass), and this stands behind it for any family whose mass
  # breaks down somewhere its author did not foresee.
  is_discrete <- S7::S7_inherits(distrib, discrete_distrib)
  impossible <- 0L

  res <- NULL
  for (eta0 in starts) {
    if (!is.finite(nll(eta0))) next
    # An error raised inside the optimizer must be treated like a failure to
    # converge, not propagated: at an awkward parameter value the quadrature
    # behind a numerically-approximated expected Hessian can fail outright
    # ("the integral is probably divergent"), and without this the random
    # restarts and the BFGS fallback promised below never get their turn.
    this <- tryCatch(run_chosen(eta0), error = function(e) NULL)
    # An explicitly chosen optimizer is not silently replaced; the fallback
    # belongs to the two named second-order strategies, which can fail on a
    # Hessian the distribution cannot supply at that point.
    if (method %in% c("fisher", "newton") &&
        (is.null(this) || !isTRUE(this$converged))) {
      alt <- tryCatch(run_bfgs(eta0), error = function(e) NULL)
      # The fallback is preferred when it converges, and kept as a last resort
      # when the chosen method raised: discarding a run that reached a point
      # merely because it did not satisfy the stopping rule reports "failed
      # from every starting value" for a fit that exists. A family whose score
      # contains a finite difference cannot drive the gradient below the
      # difference's own error, so this is reachable on correct code.
      if (!is.null(alt) && (isTRUE(alt$converged) || is.null(this))) this <- alt
    }
    # The run is discarded before it is compared, not ranked below the others:
    # the objective is what breaks ties, and a number that is not a
    # log-likelihood wins every tie it is allowed to enter.
    if (!is.null(this) && is_discrete && -this$value * scale_n > 1e-8) {
      impossible <- impossible + 1L
      this <- NULL
    }
    # Keep the BEST run, not the last one. Several starting values exist
    # precisely because one of them may end somewhere poor, and overwriting on
    # every pass reported whichever start happened to come last -- which on a
    # four-dimensional gaussian meant a random start's answer instead of the
    # maximum likelihood estimate. Convergence comes first and the objective
    # breaks ties: a converged run is what the fit promises, and a lower value
    # reached without meeting the stopping rule does not replace one.
    if (!is.null(this)) {
      better <- is.null(res) ||
        (isTRUE(this$converged) && !isTRUE(res$converged)) ||
        (isTRUE(this$converged) == isTRUE(res$converged) &&
           this$value < res$value)
      if (better) res <- this
    }
    if (!is.null(res) && isTRUE(res$converged)) break
  }

  if (is.null(res)) {
    if (impossible > 0L) {
      stop("Every run reached a positive log-likelihood, which a discrete ",
           "family cannot have: its mass function is breaking down at the ",
           "parameters reached. Supply 'start', or bound the parameters.",
           call. = FALSE)
    }
    stop("Optimization failed from every starting value; supply 'start'.", call. = FALSE)
  }

  # --- inference at the optimum -------------------------------------------
  eta_hat <- res$eta
  theta_hat <- fit_theta_from_eta(distrib, eta_hat)

  # Which information the standard errors are read off. The expected one is
  # taken when the fit itself used it, which is Fisher scoring, or when the
  # family writes it out and it therefore costs one evaluation; otherwise the
  # observed one, which every family has.
  #
  # The condition used to be `!identical(method, "newton")`, a test against a
  # STRING, and `method` has accepted an optimizer OBJECT since the
  # delegation to optimizers7. An object is normalized to "custom" above, so
  # `newton()` failed that test and the expected information was assembled
  # anyway -- by quadrature, for a family that has no closed form. Measured on
  # a user-defined Gompertz: `method = "newton"` fits in 0.15 s and
  # `method = newton()`, the same algorithm on the same data, had not
  # returned after five minutes. The `tryCatch` below does not help, a
  # quadrature that fails to converge raising nothing.
  use_expected <- identical(method, "fisher") ||
    (!identical(method, "newton") && has_exact_expected_hessian(distrib))

  # The estimates stand on their own; if the information cannot be evaluated at
  # the optimum the fit is still returned, with a missing variance matrix rather
  # than an error that throws the estimates away too.
  I_eta <- tryCatch(
    -fit_hess_matrix(distrib, y, theta_hat, expected = use_expected),
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
    criterion = if (is.null(res$criterion_met)) "" else res$criterion_met,
    note = if (is.null(res$message)) "" else res$message,
    counts = res$counts,
    score = if (is.null(res$score)) NA_real_ else res$score,
    elapsed = spent,
    level = level
  )
}

#' Render a Duration With a Unit Matched to Its Size
#'
#' @description
#' Formats a time in seconds using milliseconds below a second, seconds below a
#' minute, and minutes and seconds above, so that the figure a fit prints stays
#' readable whatever the sample size.
#'
#' @param sec A single non-negative number of seconds.
#'
#' @return A character string.
#'
#' @seealso \code{\link{print.distrib_fit}}
#' @keywords internal
fit_format_elapsed <- function(sec) {
  if (!length(sec) || !is.finite(sec)) return(NA_character_)
  if (sec < 1)  return(sprintf("%.3g ms", sec * 1000))
  if (sec < 60) return(sprintf("%.3g s", sec))
  sprintf("%d min %02d s", as.integer(sec %/% 60), as.integer(round(sec %% 60)))
}

#' Print Method for Maximum-Likelihood Fits
#'
#' @name print.distrib_fit
#' @param x A \code{\link{distrib_fit}} object.
#' @param digits Number of significant digits. Defaults to 4.
#' @param ... Unused.
#' @return \code{x}, invisibly.
S7::method(print, distrib_fit) <- function(x, digits = 4, ...) {
  lo <- paste0(format((1 - x@level) / 2 * 100, trim = TRUE), "%")
  hi <- paste0(format((1 + x@level) / 2 * 100, trim = TRUE), "%")
  mv <- S7::S7_inherits(x@distrib, multivariate_distrib)

  cat("Maximum-likelihood fit: ", x@distrib@distrib_name, "\n", sep = "")
  cat("Observations: ", x@n,
      "   Log-likelihood: ", format(x@loglik, digits = digits),
      "   AIC: ", format(x@aic, digits = digits),
      "   BIC: ", format(x@bic, digits = digits), "\n", sep = "")

  # What the optimizer did, in the shape optimizers7 reports it. The stopping
  # rule that ended the run is the thing that says what "converged" means here,
  # and a run that stopped without meeting one is worth reading for the same
  # reason.
  cat("Method: ", x@method, "   iterations: ", x@iterations, sep = "")
  if (!is.null(x@counts) && all(c("f", "g") %in% names(x@counts))) {
    cat("   evaluations: f ", x@counts[["f"]], ", g ", x@counts[["g"]], sep = "")
  }
  # Printed whenever it was measured, zero included: on a coarse clock a fast
  # fit really does report 0, and that is the reading rather than the absence
  # of one.
  if (length(x@elapsed) && is.finite(x@elapsed)) {
    cat("   time: ", fit_format_elapsed(x@elapsed), sep = "")
  }
  cat("\n")
  crit <- if (length(x@criterion) && nzchar(x@criterion)) x@criterion else "no rule reported"
  cat(if (x@converged) "Converged: yes (" else "Converged: NO (", crit, ")\n", sep = "")
  # A run that did not converge is worth reading for one number: how close to
  # stationary it ended. Printing it beside the failure saves the reader from
  # recomputing the score to find out whether the point is usable.
  if (!x@converged && length(x@score) && is.finite(x@score)) {
    cat("Score per observation at the reported point: ",
        format(x@score, digits = 3), "\n", sep = "")
  }
  if (length(x@note) && nzchar(x@note)) cat("Note: ", x@note, "\n", sep = "")
  cat("\n")

  # --- the estimates -------------------------------------------------------
  # For a multivariate family the parameters are coordinates of a covariance
  # structure, and nobody reads a Cholesky coordinate. What goes here is what
  # the model is about, which mv_summary() assembles from mv_derived(): the
  # location, and whatever the structure's matrix decomposes into. A structure
  # with fewer free values produces fewer quantities, so the block is as small
  # as the model is.
  if (mv) {
    tab_mv <- tryCatch(mv_summary(x), error = function(e) NULL)
  } else {
    tab_mv <- NULL
  }

  if (!is.null(tab_mv)) {
    loc <- tryCatch(
      {
        th <- as.list(x@coefficients)
        v <- mv_location(x@distrib, th)
        p <- x@distrib@n_dim
        m <- cbind(
          Estimate = unname(v), `Std. Error` = unname(x@se[seq_len(p)]),
          unname(x@ci[seq_len(p), , drop = FALSE])
        )
        rownames(m) <- names(v)
        colnames(m) <- c("Estimate", "Std. Error", lo, hi)
        m
      },
      error = function(e) NULL
    )
    if (!is.null(loc)) {
      cat("Location:\n")
      print(round(loc, digits))
      cat("\n")
    }
    names(tab_mv)[3:4] <- c(lo, hi)
    blocks <- attr(tab_mv, "block")
    for (b in unique(blocks)) {
      cat(b, ":\n", sep = "")
      print(round(tab_mv[blocks == b, , drop = FALSE], digits))
      cat("\n")
    }
    # Any parameter the structure does not account for -- the degrees of
    # freedom of a t, for instance -- still belongs in the table.
    extra <- setdiff(
      x@distrib@params,
      c(paste0("mu", seq_len(x@distrib@n_dim)),
        grep("^(sigma|omega)_", x@distrib@params, value = TRUE))
    )
    if (length(extra)) {
      tab_x <- cbind(Estimate = x@coefficients[extra],
                     `Std. Error` = x@se[extra], x@ci[extra, , drop = FALSE])
      colnames(tab_x) <- c("Estimate", "Std. Error", lo, hi)
      cat("Other parameters:\n")
      print(round(tab_x, digits))
      cat("\n")
    }
  } else {
    tab <- cbind(Estimate = x@coefficients, `Std. Error` = x@se, x@ci)
    colnames(tab) <- c("Estimate", "Std. Error", lo, hi)
    cat("Parameter scale:\n")
    print(round(tab, digits))
    cat("\n")
  }

  # --- the scale the fit was actually computed on --------------------------
  # The link-scale interval is the one computed; everything above is its image.
  # Showing it makes the mapping visible. When every link is the identity the
  # two tables would be the same numbers twice, so the block is dropped and the
  # reason said once -- which is what a multivariate gaussian gives, its
  # constraint living in the structure rather than in a link.
  links <- vapply(x@distrib@params,
                  function(pp) x@distrib@link_params[[pp]]@link_name, character(1))
  if (all(links == "identity")) {
    cat("The link scale is the parameter scale: every link is the identity.\n")
  } else {
    tab_eta <- cbind(Estimate = x@eta, `Std. Error` = x@se_eta, x@ci_eta)
    colnames(tab_eta) <- c("Estimate", "Std. Error", lo, hi)
    cat("Link scale:\n")
    print(round(tab_eta, digits))
  }

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

#' Confidence Intervals for a Maximum-Likelihood Fit
#'
#' @name confint.distrib_fit
#' @description
#' Returns Wald intervals. They are built symmetrically on the link scale, where
#' the parameters are unconstrained, and mapped through \eqn{g^{-1}} when the
#' parameter scale is requested, so that a limit can never leave the parameter's
#' domain. The two ends are sorted after mapping, because a link need not be
#' increasing.
#'
#' @param object A \code{\link{distrib_fit}} object.
#' @param parm Parameters to report, given by name or position. Defaults to all.
#' @param level Confidence level. Defaults to the level the fit was computed at,
#'   and any other value is obtained from the stored estimates and standard
#'   errors without refitting.
#' @param scale Either \code{"parameter"} (default) or \code{"link"}.
#' @param ... Unused.
#' @return A two-column matrix of confidence limits, one row per parameter.
#' @importFrom stats confint qnorm
S7::method(confint, distrib_fit) <- function(object, parm, level = object@level,
                                             scale = c("parameter", "link"), ...) {
  scale <- match.arg(scale)
  params <- object@distrib@params

  z <- stats::qnorm(1 - (1 - level) / 2)
  ci_eta <- cbind(lower = object@eta - z * object@se_eta,
                  upper = object@eta + z * object@se_eta)
  rownames(ci_eta) <- params

  out <- if (scale == "link") {
    ci_eta
  } else {
    m <- t(vapply(seq_along(params), function(i) {
      lk <- object@distrib@link_params[[params[i]]]
      sort(c(linkfunctions7::linkinv(lk, ci_eta[i, 1]),
             linkfunctions7::linkinv(lk, ci_eta[i, 2])), na.last = TRUE)
    }, numeric(2)))
    dimnames(m) <- list(params, c("lower", "upper"))
    m
  }

  if (!missing(parm)) {
    idx <- if (is.character(parm)) match(parm, params) else as.integer(parm)
    if (anyNA(idx) || any(idx < 1L | idx > length(params))) {
      stop("'parm' does not name a parameter of this fit.", call. = FALSE)
    }
    out <- out[idx, , drop = FALSE]
  }
  colnames(out) <- paste0(
    format(c(1 - level, 1 + level) / 2 * 100, trim = TRUE), "%")
  out
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
#' @param seed Optional seed. If supplied it is used to initialize the generator,
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
#' fit <- fit_distrib(gaussian1_distrib(), y)
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

  # The seed protocol of stats::simulate: honor `seed`, leave the caller's
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
#' continuous distribution the observations are summarized by a kernel density
#' estimate, with the fitted density drawn on top and a rug of the data
#' underneath. For a discrete one the observed relative frequencies are drawn as
#' bars with the fitted probability mass overlaid, since a kernel density would
#' misrepresent a discrete sample.
#'
#' @param x A \code{\link{distrib_fit}} object.
#' @param n_grid Number of points at which the fitted density is evaluated
#'   (continuous distributions only). Defaults to 512.
#' @param rug Logical; draw a rug of the observations. Defaults to \code{TRUE}
#'   when there are at most 2000 of them.
#' @param legend Logical; add a legend. Defaults to \code{TRUE}.
#' @param col_fit,col_data Colors of the fitted curve and of the empirical
#'   summary.
#' @param mv_which For a multivariate fit, which coordinates to show. Defaults
#'   to all of them, and at most three are drawn: above that the panel matrix
#'   stops being readable.
#' @param ... Further arguments passed to \code{\link[graphics]{plot}}, for
#'   instance \code{main}, \code{xlab} or \code{xlim}.
#'
#' @details
#' A univariate fit is drawn as the fitted density over a histogram or, for a
#' discrete family, over the observed proportions. A multivariate one is drawn
#' as a panel matrix: the fitted marginal density and a kernel estimate of the
#' data on the diagonal, the fitted contours over the observations below it,
#' and the fitted correlation above.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{fit_distrib}}, \code{\link{simulate.distrib_fit}},
#'   \code{\link{plot.multivariate_distrib}}
#'
#' @examples
#' set.seed(1)
#' y <- rgamma(300, shape = 4, rate = 2)
#' fit <- fit_distrib(gamma2_distrib(), y)
#' plot(fit)
#'
#' @importFrom graphics lines legend rug barplot points
S7::method(plot, distrib_fit) <- function(x, n_grid = 512, rug = NULL,
                                          legend = TRUE,
                                          col_fit = "#B22222", col_data = "#4682B4",
                                          mv_which = NULL,
                                          ...) {
  y <- x@y
  if (!length(y)) {
    stop("This fit carries no data to plot.", call. = FALSE)
  }
  theta <- as.list(x@coefficients)
  d <- x@distrib

  # A fitted multivariate density is a panel matrix rather than a curve: the
  # fitted marginal against a kernel estimate on the diagonal, and the fitted
  # contours over the observations off it. The kernel estimate is the
  # comparison that does not assume the model.
  if (S7::S7_inherits(d, multivariate_distrib)) {
    mv_pairs_panels(d, align_theta(d, theta), which = mv_which,
      n_grid = 80L, col_fit = col_fit,
      data = as_mv_matrix(d, y), col_data = col_data
    )
    return(invisible(x))
  }
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
