#' @include distrib.R generics.R utility_functions.R link_scale.R
NULL

# --- internal helpers ------------------------------------------------------

#' @title Parameters Carried to the Link Scale
#'
#' @description
#' Applies each parameter's own link, \eqn{\eta_i = g_i(\theta_i)}, and returns
#' the resulting vector of linear predictors. A starting value expressed in
#' natural parameters enters the optimizer through this map, because the
#' optimizer works on \eqn{\eta}, which is unconstrained.
#'
#' Only the first element of each component of `theta` is read. A fit estimates
#' one \eqn{\theta} for the whole sample, so a component of length \eqn{n} is a
#' vector that the caller built for a density evaluation, not \eqn{n} separate
#' parameters: `list(mu = c(1, 99, 99), sigma = 2)` gives the same \eqn{\eta} as
#' `list(mu = 1, sigma = 2)`.
#'
#' @param distrib An object inheriting from `distrib`, supplying `params` and
#'   `link_params`.
#' @param theta A named list of parameters on the parameter scale, in any
#'   order; the components are read by position in `distrib@params`, so it must
#'   already have been through [align_theta()]. A value outside the parameter's
#'   domain reaches the link, which returns `NaN` or a non-finite value there.
#'
#' @return An unnamed numeric vector of length `length(distrib@params)`, one
#'   linear predictor per parameter, in the order given by `distrib@params`.
#'
#' @seealso [fit_theta_from_eta()], the inverse; [fit_dtheta_deta()] for the
#'   Jacobian between the two scales; [fit_distrib()], the caller.
#' @keywords internal
fit_eta_from_theta <- function(distrib, theta) {
  params <- distrib@params
  vapply(seq_along(params), function(i) {
    linkfunctions7::linkfun(distrib@link_params[[params[i]]], theta[[i]])[1]
  }, numeric(1))
}

#' @title The Link Scale Carried Back to Parameters
#'
#' @description
#' Applies each parameter's inverse link, \eqn{\theta_i = g_i^{-1}(\eta_i)},
#' and returns the parameters as the named list that every generic in the
#' package expects. It is the inverse of [fit_eta_from_theta()].
#'
#' Every link maps onto the open interior of its parameter's domain, and
#' `linkfunctions7` clamps the result strictly inside when the arithmetic
#' saturates, so a \eqn{\theta} obtained this way is admissible whatever the
#' optimizer proposed. On a Gaussian scale carried by the logarithm,
#' \eqn{\eta = -800} gives \eqn{\sigma = 1.9\times 10^{-77}} and
#' \eqn{\eta = 800} gives \eqn{1.8\times 10^{308}}: both values are extreme
#' and both lie inside \eqn{(0, \infty)}, so the density can be evaluated and
#' the point rejected on its likelihood. For the same reason a confidence
#' interval built on the link scale and mapped back cannot run outside the
#' domain.
#'
#' @param distrib An object inheriting from `distrib`, supplying `params` and
#'   `link_params`.
#' @param eta A numeric vector of linear predictors, one per parameter, in the
#'   order given by `distrib@params`. A non-finite entry propagates to the
#'   parameter that it names.
#'
#' @return A named list of length `length(distrib@params)`, in which each
#'   component is a number on the parameter scale, named and ordered as
#'   `distrib@params`.
#'
#' @seealso [fit_eta_from_theta()], the inverse; [fit_dtheta_deta()] for the
#'   first derivative of the same map; [linkfunctions7::linkinv()].
#' @keywords internal
fit_theta_from_eta <- function(distrib, eta) {
  params <- distrib@params
  out <- lapply(seq_along(params), function(i) {
    linkfunctions7::linkinv(distrib@link_params[[params[i]]], eta[i])
  })
  names(out) <- params
  out
}

#' @title Jacobian of the Inverse Link at the Estimate
#'
#' @description
#' Returns \eqn{h_i'(\eta_i) = dg_i^{-1}/d\eta_i}, one entry per parameter,
#' evaluated at the supplied linear predictors. Because each parameter carries
#' its own scalar link the Jacobian of \eqn{\theta} in \eqn{\eta} is diagonal,
#' and this is its diagonal.
#'
#' The delta method uses it to carry a variance matrix from the link
#' scale, where the fit computes it, to the parameter scale, where it is
#' reported: \eqn{\widehat{\mathrm{Var}}(\hat\theta) = J V J} with
#' \eqn{J = \mathrm{diag}(h')}. On a Gaussian at \eqn{\eta = (1.5, \log 2.5)},
#' with the identity on \eqn{\mu} and the logarithm on \eqn{\sigma}, it returns
#' `c(1, 2.5)`.
#'
#' @param distrib An object inheriting from `distrib`, supplying `params` and
#'   `link_params`.
#' @param eta A numeric vector of linear predictors, one per parameter, in the
#'   order of `distrib@params`. Only `eta[i]` is read for the i-th
#'   parameter, so entries beyond the parameter count are ignored.
#'
#' @return An unnamed numeric vector of length `length(distrib@params)`, in the
#'   order of `distrib@params`.
#'
#' @seealso [fit_theta_from_eta()], the map that this function differentiates;
#'   [fit_distrib()], which uses it for the standard errors and nothing else;
#'   [linkfunctions7::dlinkinv()].
#' @keywords internal
fit_dtheta_deta <- function(distrib, eta) {
  params <- distrib@params
  vapply(seq_along(params), function(i) {
    linkfunctions7::dlinkinv(distrib@link_params[[params[i]]], eta[i])[1]
  }, numeric(1))
}

#' @title Summed Score on the Link Scale
#'
#' @description
#' Returns \eqn{\partial \ell / \partial \eta} for the whole sample: the
#' package's per-observation gradient on the link scale, summed over
#' observations, as a plain numeric vector in parameter order. At the maximum
#' likelihood estimate every entry is zero, up to the tolerance of the
#' stopping rule that ended the fit.
#'
#' @param distrib An object inheriting from `distrib`.
#' @param y A numeric vector of observations, or the response matrix of a
#'   multivariate family.
#' @param theta A named list of parameters on the parameter scale, aligned to
#'   `distrib@params`.
#' @param threads How many threads the family's compiled kernels may use, as an
#'   integer count. Defaults to `1L`, which takes the sequential path. The sum
#'   does not depend on the count.
#'
#' @return An unnamed numeric vector of length `length(distrib@params)`, in the
#'   order of `distrib@params`.
#'
#' @seealso [fit_hess_matrix()] for the second-order counterpart;
#'   [distrib_gradient()], which supplies the per-observation components;
#'   [fit_distrib()], the caller.
#' @keywords internal
fit_score <- function(distrib, y, theta, threads = 1L) {
  g <- distrib_gradient(distrib, y, theta, scale = "link", threads = threads)
  vapply(g, function(v) sum(v), numeric(1))
}

#' @title Summed Hessian on the Link Scale, as a Matrix
#'
#' @description
#' Assembles the package's named list of second-derivative components into the
#' symmetric \eqn{p \times p} matrix that an optimizer requires, summed over
#' observations and on the link scale. Components are stored one per unordered
#' index pair, so this fills both triangles from the one value.
#'
#' With `expected = TRUE` it assembles the expected Hessian instead, which turns
#' Newton's method into Fisher scoring and makes a fit possible on a family
#' whose observed curvature is unusable. On a Laplace at
#' \eqn{\sigma = 1} with 400 draws the observed matrix is
#' \eqn{\bigl(\begin{smallmatrix} 0 & 2\\ 2 & -392\end{smallmatrix}\bigr)}:
#' \eqn{\partial^2\ell/\partial\mu^2} is zero almost everywhere, so the
#' determinant is \eqn{-(\sum_i \mathrm{sign}(y_i - \mu))^2/\sigma^2}, which is
#' negative unless the signs balance exactly, and the matrix is indefinite.
#' The expected information at the same point is \eqn{-400 I}.
#'
#' @param distrib An object inheriting from `distrib`.
#' @param y A numeric vector of observations, or the response matrix of a
#'   multivariate family.
#' @param theta A named list of parameters on the parameter scale, aligned to
#'   `distrib@params`.
#' @param expected Logical of length 1, with no default. `TRUE` assembles the
#'   expected Hessian and `FALSE` the observed one.
#' @param method `NULL` (the default) or an object inheriting from
#'   [ExpectMethod] that selects how the expected information is obtained:
#'   `NULL` takes the closed form of the family where it has one and
#'   [expect_opg()] otherwise, and an object is used for every family; see
#'   [expected_derivative_methods]. It is read only when `expected` is `TRUE`.
#' @param threads How many threads the family's compiled kernels may use, as an
#'   integer count. Defaults to `1L`. The matrix does not depend on the count.
#'
#' @return A symmetric numeric matrix of dimension
#'   `length(distrib@params)`, whose dimnames are both `distrib@params`.
#'
#' @seealso [fit_score()] for the first-order counterpart;
#'   [distrib_hessian()] and [distrib_expected_hessian()], which supply the
#'   components; [fisher_scoring()], where `method` is set.
#' @keywords internal
fit_hess_matrix <- function(distrib, y, theta, expected,
                            method = NULL, threads = 1L) {
  params <- distrib@params
  p <- length(params)
  h <- if (expected) {
    distrib_expected_hessian(distrib, y, theta,
      scale = "link", method = method, threads = threads
    )
  } else {
    distrib_hessian(distrib, y, theta, scale = "link", threads = threads)
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

#' Error Message for a Method That Is Not an Object
#'
#' @description
#' Builds the message of [fit_distrib()] for a `method` that is neither a
#' [FisherScoring] object nor an optimizer object.
#'
#' @param method The rejected value.
#'
#' @return A single character string.
#'
#' @seealso [fit_distrib()], [fisher_scoring()]
#' @keywords internal
method_error <- function(method) {
  sprintf(paste0(
    "'method' must be fisher_scoring() or an optimizers7 optimizer object ",
    "such as optimizers7::newton() or optimizers7::bfgs(), and it is '%s'."),
    paste(class(method), collapse = "/"))
}

#' @title Total Log-Likelihood
#'
#' @description
#' Sums the log-density over the observations, \eqn{\sum_i \log f(y_i;
#' \theta)}, and returns the one number that the fit maximizes. Every
#' observation is read at the same \eqn{\theta}, which is the i.i.d. assumption
#' that the fitting layer makes.
#'
#' The value is not divided by \eqn{n}. [fit_distrib()] scales it when it
#' hands the objective to an optimizer, and recomputes it unscaled here for the
#' `loglik`, `aic` and `bic` that a fit reports.
#'
#' @param distrib An object inheriting from `distrib`.
#' @param y A numeric vector of observations, or the response matrix of a
#'   multivariate family.
#' @param theta A named list of parameters on the parameter scale, aligned to
#'   `distrib@params`.
#'
#' @return A single number. It is `-Inf` when any observation has zero density,
#'   and `NaN` when the density itself is not computable at `theta`.
#'
#' @seealso [distrib_pdf()], which supplies the terms;
#'   [logLik.distrib_fit()] for the value that a fit carries.
#' @keywords internal
fit_loglik <- function(distrib, y, theta) {
  sum(distrib_pdf(distrib, y, theta, log = TRUE))
}

#' @title S7 Class for Maximum-Likelihood Fits
#' @name distrib_fit_class
#'
#' @description
#' This is the class of the object that [fit_distrib()] returns. It holds the
#' estimates on both the parameter scale and the link scale, their standard
#' errors and confidence limits, the maximized log-likelihood with the two
#' information criteria built on it, and a record of what the optimizer did.
#'
#' Both scales are kept because the fit is computed on one and read on the
#' other. The variance matrix is the inverse of the observed information at
#' \eqn{\hat\eta}; `vcov` is its image under the delta method and `ci` is
#' `ci_eta` mapped through \eqn{g^{-1}}, so a limit on the parameter scale
#' respects the parameter's domain by construction.
#'
#' This page documents the raw S7 constructor. It validates nothing and is not
#' the way to build a fit; call [fit_distrib()].
#'
#' @param distrib The fitted `distrib` object, carrying the parametrization and
#'   the links in which the estimates are expressed.
#' @param y The observations from which the fit was computed, kept so that the
#'   fitted distribution can be compared with the data by [plot.distrib_fit()]
#'   and so that the information can be recomputed by [vcov.distrib_fit()].
#' @param n The number of observations: the row count for a multivariate
#'   response and the length otherwise. `bic` and the printed header both
#'   read it.
#' @param coefficients Named numeric estimates on the parameter scale, one per
#'   parameter, in the order given by `distrib@params`.
#' @param se Named standard errors on the parameter scale, from the delta
#'   method applied to the inverse of the observed information. An entry is `NaN` where the corresponding variance came back
#'   negative and `NA` where it is missing.
#' @param ci A two-column numeric matrix of confidence limits on the parameter
#'   scale, one row per parameter, columns `lower` and `upper`.
#' @param eta Named estimates on the link scale, \eqn{\hat\eta = g(\hat\theta)}.
#' @param se_eta Named standard errors on the link scale, the square roots of
#'   the diagonal of `vcov_eta`. These are the standard errors that the fit
#'   computes; `se` is derived from them. Another information matrix is
#'   selected after the fit through [vcov.distrib_fit()].
#' @param ci_eta A two-column numeric matrix of confidence limits on the link
#'   scale, symmetric about `eta`.
#' @param vcov The variance-covariance matrix on the parameter scale, with both
#'   dimnames set to the parameter names.
#' @param aic,bic Information criteria built on the unscaled log-likelihood,
#'   \eqn{-2\ell + 2p} and \eqn{-2\ell + p\log n} respectively, with \eqn{p}
#'   the number of estimated parameters.
#' @param vcov_eta The variance-covariance matrix on the link scale, the
#'   inverse of the information at \eqn{\hat\eta}. Every entry is `NA` when the
#'   information could not be evaluated or inverted there; the estimates stand
#'   in that case and only the uncertainty is missing.
#' @param loglik The maximized log-likelihood, summed over observations and
#'   not divided by \eqn{n}.
#' @param iterations How many iterations the run that was kept took.
#' @param converged Logical of length 1: whether a stopping rule confirmed
#'   convergence. A run that exhausted its iteration budget is not converged
#'   whatever point it reached.
#' @param method The optimization method actually used, which is not always the
#'   one requested: `"Fisher scoring"` and `"Newton-Raphson"` fall back to
#'   `"BFGS"` when they fail, and this records which one produced the estimates.
#' @param criterion Which stopping rule ended the run, in the words in which
#'   \pkg{optimizers7} reports it, such as
#'   `"gradient (max-norm) < 1e-06"`. It is empty when no rule fired.
#' @param note Any remark that the optimizer attached to the run, such as
#'   `"the line search found no acceptable step"`. It is empty when there is
#'   no remark.
#' @param counts A named list of evaluation counts with components `f`, `g` and
#'   `h`: how many times the objective, its gradient and its Hessian were
#'   evaluated in the run that was kept.
#' @param score The max-norm of the score per observation at the reported
#'   optimum. This is the quantity that the stopping rule tested, so it
#'   measures how close to stationary a run ended, and it is the number to read
#'   on a fit that did not converge.
#' @param elapsed Seconds spent optimizing, summed over every starting value
#'   and every fallback attempted, not just the run that was kept.
#' @param level The confidence level at which `ci` and `ci_eta` were built.
#'
#' @section Methods:
#' The following methods are registered on this class:
#'   [`coef()`][coef.distrib_fit],
#'   [`confint()`][confint.distrib_fit],
#'   [`logLik()`][logLik.distrib_fit],
#'   [`plot()`][plot.distrib_fit],
#'   [`print()`][print.distrib_fit],
#'   [`simulate()`][simulate.distrib_fit],
#'   [`vcov()`][vcov.distrib_fit]
#'
#' `coef()`, `vcov()` and `confint()` each take a `scale` argument and report
#' either scale from the stored components.
#'
#' @return An S7 object of class `distrib_fit`, with the properties above.
#'
#' @examples
#' set.seed(1)
#' d <- gaussian1_distrib()
#' y <- distrib_rng(d, 200, list(mu = 1, sigma = 2))
#' fit <- fit_distrib(d, y)
#' S7::S7_inherits(fit, distrib_fit)
#'
#' # The two scales are stored together, and eta is the link of the estimate.
#' rbind(parameter = coef(fit), link = coef(fit, scale = "link"))
#' all.equal(fit@eta[["sigma"]], log(fit@coefficients[["sigma"]]))
#'
#' # The parameter-scale interval is the link-scale one mapped through g^-1,
#' # so the lower limit of a scale stays positive whatever the sample.
#' fit@ci
#' all.equal(fit@ci[["sigma", "lower"]], exp(fit@ci_eta[["sigma", "lower"]]))
#'
#' # The fit records what the optimizer did. 'score' is the max-norm of the
#' # score per observation at the point reported, which the rule tested.
#' c(iterations = fit@iterations, converged = fit@converged, score = fit@score)
#' fit@criterion
#'
#' @seealso [fit_distrib()], which builds one;
#'   [print.distrib_fit()] for the printed layout;
#'   [confint.distrib_fit()] to recompute an interval at another level.
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

#' @title Maximum-Likelihood Estimation for a Distribution
#'
#' @description
#' Fits a `distrib` object to an i.i.d. sample by maximum likelihood and
#' returns a [distrib_fit()] carrying the estimates, their standard errors and
#' confidence limits on both scales. The optimization runs on the link
#' scale, where the parameters are unconstrained, and it is driven by the
#' analytical score and information that the family supplies at
#' `scale = "link"`; the estimates are then mapped back and reported on the
#' parameter scale.
#'
#' On a Gaussian the answer is the closed-form estimate to the printed digit:
#' \eqn{\hat\mu = \bar y}, \eqn{\hat\sigma^2 = \frac{1}{n}\sum (y_i - \bar
#' y)^2}, with \eqn{\mathrm{se}(\hat\mu) = \hat\sigma/\sqrt{n}} and
#' \eqn{\mathrm{se}(\hat\sigma) = \hat\sigma/\sqrt{2n}}.
#'
#' @param distrib An object inheriting from `distrib`: a univariate family, a
#'   multivariate one, or any wrapper of either.
#' @param y A numeric vector of observations, or an \eqn{n \times p} matrix for
#'   a multivariate family. Every observation is read at the same \eqn{\theta}.
#' @param start Optional named list of starting values on the parameter
#'   scale. `NULL`, the default, takes the starting values from
#'   [distrib_start()], which computes them from the data for a family that
#'   registers a method and otherwise draws them from the parameter domains,
#'   with the restarts below. A value that is
#'   neither `NULL` nor a list signals an error naming the argument, because
#'   `start` sits before `method` in the signature and an optimizer passed
#'   positionally lands here.
#' @param method An object that selects how to optimize. The argument takes one
#'   of two kinds of object:
#'
#'   - [fisher_scoring()], the default: Newton's method with the expected
#'     information in place of the observed Hessian, the object carrying the
#'     method ([ExpectMethod]) that obtains that information. With
#'     `method = NULL` the closed form of the family is used where there is
#'     one, and [expect_opg()] otherwise; a method given there is used for
#'     every family. Fisher scoring falls back to BFGS if it fails to
#'     converge;
#'   - an optimizer object from \pkg{optimizers7}, used as given and receiving
#'     the analytical gradient and the observed Hessian, so that
#'     `method = optimizers7::lbfgs(criterion = optimizers7::crit_grad(1e-12))`
#'     selects the algorithm and the stopping rule together. A stopping rule
#'     that the optimizer cannot evaluate is rejected here, where the message
#'     can name it. An optimizer that the caller chose is never replaced.
#'
#'   Any other value, a character string included, signals an error.
#'
#'   The iteration limit and the stopping rule belong to the method and are set
#'   there, on an optimizer through its own `maxit` and `criterion` and on
#'   [fisher_scoring()] through the same two arguments. For
#'   [fisher_scoring()] the rule is [optimizers7::crit_grad()] at its default
#'   tolerance when the method does not set a rule; an optimizer object keeps
#'   its own default.
#' @param level Confidence level for `ci` and `ci_eta`, a single number in
#'   \eqn{(0, 1)}; any other value signals an error. Defaults to 0.95. Any
#'   other level is available afterwards from [confint.distrib_fit()] without
#'   refitting.
#' @param n_start How many starting values to request from [distrib_start()] when
#'   `start` is `NULL`. Defaults to 5. A family that returns its own estimate
#'   returns one and ignores this. It is ignored entirely when `start` is
#'   given.
#' @param threads How many threads the fit may use inside this R process, a
#'   single whole number (see [numericals7::local_threads()]). The default, `1`,
#'   is sequential and takes the sequential code path. The count reaches the
#'   family's compiled per-observation kernels as an argument; the result does
#'   not depend on it, bit for bit, because every parallel region decomposes
#'   its work over the elements of its output and never splits a reduction.
#'
#' @return An S7 object of class [distrib_fit()]. Its `coefficients`, `se` and
#'   `ci` are on the parameter scale and its `eta`, `se_eta` and `ci_eta` on
#'   the link scale; `loglik`, `aic` and `bic` are computed from the unscaled
#'   log-likelihood at the estimate, and `converged`, `score`, `iterations`,
#'   `method`, `criterion` and `counts` record what the optimizer did. See that
#'   page for every component.
#'
#'   Signals an error when every starting value fails to produce a usable run,
#'   and a separate one when every run of a discrete family reached a positive
#'   log-likelihood, which is impossible for a product of probabilities and
#'   indicates that the mass function has broken down at the parameters
#'   reached.
#'
#' @details
#' # Why the link scale
#' Optimizing \eqn{\eta \in \mathbb{R}^p} in place of the constrained
#' \eqn{\theta} removes the need for box constraints: a scale cannot become
#' negative and a probability cannot leave \eqn{(0,1)}, because every link maps
#' onto the interior of its parameter's domain. The score and the information
#' on that scale follow from those on the parameter scale through the chain
#' rule of [link_scale_derivatives()], and they are exact wherever the
#' family's own derivatives are.
#'
#' # Scaling of the objective
#' The optimizer receives \eqn{-\ell(\eta)/n}, and the gradient and
#' the Hessian are divided by \eqn{n} as well. Scaling by a positive constant
#' moves neither the maximum nor any Newton step, since \eqn{H^{-1}g} is
#' unchanged when both are divided by \eqn{n}. The scaling changes the meaning
#' of a threshold: a tolerance on the gradient of this objective is a tolerance
#' on the score per observation whatever the sample size, so the same rule means
#' the same thing at \eqn{n = 10} and at \eqn{n = 10^7}. `loglik`, the
#' information and every standard error are recomputed unscaled at the optimum.
#'
#' # Standard errors and intervals
#' The variance matrix on the link scale is the inverse of the observed
#' information at \eqn{\hat\eta}, the negative Hessian of the log-likelihood,
#' for every family and every `method`. The matrix with which the optimizer
#' ended is not used, and another matrix is selected after the fit with the
#' `information` argument of [vcov.distrib_fit()], [print.distrib_fit()] and
#' [confint.distrib_fit()] (see [information_objects]). The
#' delta method carries the matrix to the parameter scale,
#' \deqn{\widehat{\mathrm{Var}}(\hat\theta) = J\,\widehat{\mathrm{Var}}(\hat\eta)\,J,
#'       \qquad J = \mathrm{diag}\!\left(\frac{dg^{-1}}{d\eta}\Big|_{\hat\eta}\right).}
#' Intervals are built symmetrically on the link scale, \eqn{\hat\eta \pm
#' z_{1-\alpha/2}\,\mathrm{se}(\hat\eta)}, and mapped through \eqn{g^{-1}},
#' sorting the pair in case the link decreases. The limits therefore respect
#' the parameter's domain, which a symmetric interval on the parameter scale
#' would not.
#'
#' # Restarts, the fallback and the tie-break
#' Each starting value is tried in turn, and the search stops once two runs
#' that converged have reached the best objective found so far, to within
#' \eqn{10^{-8}} times the larger of 1 and the absolute value of the mean
#' negative log-likelihood. A single converged run is not taken as the answer
#' because it may be a degenerate stationary point.
#' Where fewer than two runs converge every start is tried and the best is
#' kept. Fisher scoring falls back to BFGS from the same starting value when
#' it fails; an optimizer object that the caller supplied does not. Among the runs that finish, a run that converged beats one that
#' did not converge and the objective breaks ties, so the fit reports the best
#' run and not the last. A run of a discrete family whose log-likelihood came
#' back positive is discarded before any comparison.
#'
#' @examples
#' set.seed(1)
#' d <- gaussian1_distrib()
#' y <- distrib_rng(d, 500, list(mu = 2, sigma = 3))
#' fit <- fit_distrib(d, y)
#' fit
#'
#' # The Gaussian MLE is closed form, and the fit reaches it.
#' all.equal(coef(fit)[["mu"]], mean(y))
#' all.equal(coef(fit)[["sigma"]], sqrt(mean((y - mean(y))^2)))
#'
#' # So are the standard errors: sigma/sqrt(n) and sigma/sqrt(2n).
#' s <- coef(fit)[["sigma"]]
#' all.equal(unname(fit@se), c(s / sqrt(500), s / sqrt(2 * 500)))
#'
#' # A bounded parameter: the interval is built on the link scale and mapped
#' # back, so it cannot contain a probability outside (0, 1).
#' b <- bernoulli_distrib()
#' fb <- fit_distrib(b, rbinom(50, 1, 0.9))
#' rbind(link = confint(fb, scale = "link"), parameter = confint(fb))
#'
#' # A non-regular family. The Laplace's observed curvature in the location is
#' # zero almost everywhere, so Newton's method has nothing to invert; Fisher
#' # scoring uses the information instead and reaches the closed-form estimates.
#' yl <- distrib_rng(laplace_distrib(), 400, list(mu = 0, sigma = 1))
#' fl <- fit_distrib(laplace_distrib(), yl)
#' c(fitted = coef(fl)[["mu"]], median = median(yl))
#' c(fitted = coef(fl)[["sigma"]], mad = mean(abs(yl - median(yl))))
#'
#' @seealso [distrib_fit()] for the object returned;
#'   [fisher_scoring()] for the default method;
#'   [distrib_start()] for where the starting values come from;
#'   [link_scale_derivatives()] for the chain rule that the score uses;
#'   [check_distrib()] to validate a family before fitting it;
#'   [optimizers7::minimize()], which runs the search.
#' @importFrom stats qnorm setNames
#' @export
fit_distrib <- function(distrib, y, start = NULL,
                        method = fisher_scoring(),
                        level = 0.95, n_start = 5,
                        threads = 1) {
  # The count is read once and passed DOWN as an argument; the process-level
  # RcppParallel setting is sized here and restored when this frame exits,
  # so a fit never leaves it moved for the code that runs after it. At
  # threads = 1 neither call touches anything.
  tc <- numericals7::local_threads(threads)
  # The response is a numeric vector, or a numeric matrix for a multivariate
  # family. Anything else -- the object cens() returns is the case met in
  # practice -- reached the first coercion and stopped with "cannot coerce
  # type 'object' to vector of type 'double'".
  if (!is.numeric(y) && !is.logical(y)) {
    cens_hint <- if (any(grepl("censor", tolower(class(y))))) {
      "\n  A censored response has no likelihood in fit_distrib(): fit the observed values."
    } else {
      ""
    }
    stop(sprintf("'y' must be a numeric vector, and it is '%s'.%s",
                 paste(class(y), collapse = "/"), cens_hint), call. = FALSE)
  }
  # A level outside (0, 1) reached qnorm() at the end of the fit and came back
  # as NaN intervals with only qnorm()'s warning.
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      level <= 0 || level >= 1) {
    stop("'level' must be a single number in (0, 1).", call. = FALSE)
  }
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

  # One argument says how to optimize, and it takes one of two objects: a
  # fisher_scoring() specification or an optimizers7 optimizer. How the
  # expected information is to be approximated is a property of Fisher scoring
  # and lives on that object, not among fit_distrib()'s own arguments, where it
  # would sit next to optimizers that never look at it. The code below tells
  # the two apart by class: `fs` is the specification, or NULL for an
  # optimizer.
  optimizer <- NULL
  info_method <- NULL
  fs <- NULL

  if (S7::S7_inherits(method, FisherScoring)) {
    fs <- method
    info_method <- fs@method
  } else if (S7::S7_inherits(method, optimizers7::optimizer)) {
    optimizer <- method
    # A stopping rule the optimizer cannot evaluate is a mistake in the call,
    # not a numerical failure, and it is raised here so that it reaches the
    # caller with its own explanation rather than being caught by the restart
    # loop below and reported as a fit that never converged.
    optimizers7::check_criterion(optimizer)
    # An optimizer class with no minimize() method of its own fails at every
    # start, and the restart loop below reported that as "Optimization failed
    # from every starting value", which names neither the class nor the
    # missing method.
    if (is.null(tryCatch(S7::method(optimizers7::minimize, S7::S7_class(optimizer)),
                         error = function(e) NULL))) {
      stop(sprintf(paste0(
        "'method' is an optimizer of class '%s', which has no minimize() ",
        "method.\n  Register one with S7::method(optimizers7::minimize, ",
        "<class>) before fitting."), class(optimizer)[[1L]]), call. = FALSE)
    }
  } else {
    stop(method_error(method), call. = FALSE)
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
    -fit_score(distrib, y, th, threads = tc) / scale_n
  }
  nll_he <- function(expected) {
    function(eta) {
      th <- fit_theta_from_eta(distrib, eta)
      -fit_hess_matrix(distrib, y, th,
        expected = expected, method = info_method, threads = tc
      ) / scale_n
    }
  }

  # How the run stops and how long it may take are properties of the method,
  # so they are read off the method and nowhere else. A budget the caller did
  # not set stays at the optimizer's own, which is the only place that constant
  # is written down. The objective the rule sees is the mean negative
  # log-likelihood, so a tolerance on its gradient is a tolerance on the score
  # per observation without the rule having to know n.
  #
  # ⚠️ THE RULE IS NAMED HERE AND THE CONSTANT IS NOT. From optimizers7 0.6.0
  # the gradient methods default to a disjunction that also stops on a stalled
  # objective, and this function cannot use it: the loop below reads
  # `converged` as the signal to try another start and to fall back to BFGS,
  # so a rule that reports success at a stall turns a multi-start search into a
  # single-start one that keeps the stall. Measured on
  # folded(gaussian1_distrib()) at n = 3000 with mu = 1.2, sigma = 2: under the
  # wider rule Fisher scoring reports convergence after 141 iterations at a
  # score of 0.57, with mu = 0.103 and sigma = 3.572 and a log-likelihood 413
  # units below what the same call reaches under this one. A maximum likelihood
  # fit promises a stationary point, so it asks for the rule that tests one;
  # `crit_grad()` carries no number, and its tolerance stays optimizers7's.
  opt_args <- list(criterion = optimizers7::crit_grad())
  if (!is.null(fs)) {
    if (!is.null(fs@criterion)) opt_args$criterion <- fs@criterion
    if (!is.null(fs@maxit)) opt_args$maxit <- fs@maxit
  }

  # Time is accumulated over every attempt rather than taken from the run that
  # is kept: what a caller wants to know is what the fit cost, and the restarts
  # and the fallback are part of that.
  spent <- 0

  run <- function(opt, eta0, he, label) {
    r <- optimizers7::minimize(opt, fn = nll, par = eta0, gr = nll_gr, hess = he)
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
    if (!is.null(fs)) {
      run(do.call(optimizers7::newton, opt_args),
          eta0, nll_he(TRUE), "Fisher scoring")
    } else {
      run(optimizer, eta0, nll_he(FALSE), optimizer@name)
    }
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
  conv_values <- numeric(0)
  agree_tol <- function(v) 1e-8 * max(1, abs(v))
  for (eta0 in starts) {
    if (!is.finite(nll(eta0))) next
    # An error raised inside the optimizer must be treated like a failure to
    # converge, not propagated: at an awkward parameter value the quadrature
    # behind a numerically-approximated expected Hessian can fail outright
    # ("the integral is probably divergent"), and without this the random
    # restarts and the BFGS fallback promised below never get their turn.
    this <- tryCatch(run_chosen(eta0), error = function(e) NULL)
    # An explicitly chosen optimizer is not silently replaced; the fallback
    # belongs to Fisher scoring, which can fail on a Hessian the distribution
    # cannot supply at that point.
    if (!is.null(fs) &&
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
           this$value < res$value - agree_tol(res$value))
      if (better) res <- this
      if (isTRUE(this$converged)) conv_values <- c(conv_values, this$value)
    }
    # The search stops once TWO converged runs have reached the best objective
    # found, not at the first converged run. A single converged run can be a
    # degenerate stationary point: on zero_adjusted(pig1_distrib()) the first
    # start converged at a log-likelihood of -2851.73, where nine of the ten
    # starts reach -2754.36. Running every start costs 6.5 times the first-
    # converged rule over the univariate families; confirming the best by a
    # second run costs one more fit wherever the first converged run is already
    # the maximum. The objective is -l/n, so the agreement is per observation.
    # A later run replaces the kept one only by MORE than that tolerance: two
    # runs at one maximum differ in the last bits, and letting a random start
    # win on rounding made the reported estimate depend on the random stream
    # at 1e-9 relative, which is what the thread-count twin test caught.
    if (!is.null(res) && isTRUE(res$converged) &&
        sum(abs(conv_values - res$value) <= agree_tol(res$value)) >= 2L) break
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

  # The reported variance matrix is the inverse of the default information at
  # the optimum (default_information()): the observed information when the
  # log-likelihood is twice differentiable in every parameter, and the
  # expected information of the family otherwise, for every method. The matrix
  # that the optimizer ended with is not used: a scoring step may be driven by
  # any positive definite matrix, the score being exact, so Fisher scoring on a
  # family with no closed form runs on the outer product of the observed scores
  # (see expected_by_opg()), which has the same fixed point and is a noisier
  # estimate of the information. Another matrix is selected after the fit, by
  # the `information` argument of vcov().
  #
  # The estimates stand on their own; if the information cannot be evaluated at
  # the optimum the fit is still returned, with a missing variance matrix rather
  # than an error that throws the estimates away too.
  info0 <- default_information(distrib)
  I_eta <- tryCatch(
    -fit_hess_matrix(distrib, y, theta_hat,
                     expected = S7::S7_inherits(info0, ExpectedInfo),
                     method = if (S7::S7_inherits(info0, ExpectedInfo)) info0@method),
    error = function(e) matrix(NA_real_, p, p, dimnames = list(params, params))
  )
  V_eta <- tryCatch(solve(I_eta), error = function(e) matrix(NA_real_, p, p))
  dimnames(V_eta) <- list(params, params)
  inf <- fit_inference(distrib, eta_hat, V_eta, level)

  ll_hat <- fit_loglik(distrib, y, theta_hat)
  coefs <- stats::setNames(vapply(theta_hat, function(v) v[1], numeric(1)), params)

  distrib_fit(
    distrib = distrib, y = y, n = n,
    coefficients = coefs,
    se = inf$se,
    ci = inf$ci,
    eta = stats::setNames(eta_hat, params),
    se_eta = inf$se_eta,
    ci_eta = inf$ci_eta,
    vcov = inf$vcov, vcov_eta = inf$vcov_eta,
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

#' Standard Errors and Intervals From a Link-Scale Variance Matrix
#'
#' @description
#' Carries a variance matrix of the estimates on the link scale to the
#' quantities that a fit stores and reports: the variance matrix on the
#' parameter scale by the delta method, the standard errors on both scales and
#' the Wald intervals, which are symmetric on the link scale and mapped through
#' the inverse links.
#'
#' @param distrib The fitted `distrib` object.
#' @param eta The estimates on the link scale.
#' @param V_eta The variance matrix on the link scale, named by the
#'   parameters; entries are `NA` where it could not be computed.
#' @param level The confidence level of the intervals.
#'
#' @return A list with the components `vcov`, `vcov_eta`, `se`, `se_eta`, `ci`
#'   and `ci_eta`, laid out as the properties of the same names in
#'   [distrib_fit()].
#'
#' @seealso [fit_distrib()], which stores the result for the observed
#'   information; [vcov.distrib_fit()], which recomputes it for another matrix.
#' @keywords internal
fit_inference <- function(distrib, eta, V_eta, level) {
  params <- distrib@params
  p <- length(params)
  J <- fit_dtheta_deta(distrib, eta)
  V_theta <- diag(J, p, p) %*% V_eta %*% diag(J, p, p)
  dimnames(V_theta) <- list(params, params)

  # a negative variance gives NaN, not a standard error of zero
  var_eta <- diag(V_eta)
  var_theta <- diag(V_theta)
  var_eta[!is.na(var_eta) & var_eta < 0] <- NaN
  var_theta[!is.na(var_theta) & var_theta < 0] <- NaN
  se_eta <- sqrt(var_eta)
  se_theta <- sqrt(var_theta)

  z <- stats::qnorm(1 - (1 - level) / 2)
  ci_eta <- cbind(lower = eta - z * se_eta, upper = eta + z * se_eta)
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

  list(vcov = V_theta, vcov_eta = V_eta,
       se = stats::setNames(se_theta, params),
       se_eta = stats::setNames(se_eta, params),
       ci = ci_theta, ci_eta = ci_eta)
}

#' Standard Errors and Intervals of a Fit for a Chosen Information
#'
#' @description
#' Returns the quantities of [fit_inference()] for the information selected by
#' an [Information] object, or by the default of the family when `information`
#' is `NULL`. The fit stores the quantities of the default information
#' ([default_information()]), and they are returned when the selected
#' information is the default; otherwise the link-scale variance matrix is
#' recomputed by [vcov.distrib_fit()] and carried through [fit_inference()].
#'
#' @param object A [distrib_fit()] object.
#' @param information `NULL`, an [observed_info()] or an [expected_info()]
#'   object.
#'
#' @return A list as returned by [fit_inference()].
#'
#' @seealso [vcov.distrib_fit()], [print.distrib_fit()],
#'   [confint.distrib_fit()], [mv_summary()]
#' @keywords internal
fit_report <- function(object, information = NULL) {
  information <- resolve_information(object@distrib, information)
  if (same_information(information, default_information(object@distrib))) {
    return(list(vcov = object@vcov, vcov_eta = object@vcov_eta,
                se = object@se, se_eta = object@se_eta,
                ci = object@ci, ci_eta = object@ci_eta))
  }
  V_eta <- vcov(object, scale = "link", information = information)
  fit_inference(object@distrib, object@eta, V_eta, object@level)
}

#' @title A Duration Rendered With a Unit Matched to Its Size
#'
#' @description
#' Formats a time in seconds for the line that [print.distrib_fit()] shows,
#' choosing the unit from the size: milliseconds below a second, seconds below
#' a minute, and minutes with seconds above. A fit of a few hundred
#' observations then reads `40 ms` where one of ten million reads
#' `2 min 07 s`, and in neither case does the reader have to count zeros.
#'
#' @param sec A single non-negative number of seconds. `NA`, a non-finite
#'   value and a zero-length vector all give `NA_character_`; the value is not
#'   otherwise validated, and a negative number is formatted as it stands.
#'
#' @return A character string of length 1: `"0.4 ms"`, `"40 ms"`, `"1.5 s"`,
#'   `"59.4 s"`, `"1 min 01 s"`. Milliseconds, and seconds below a minute, are
#'   given to three significant digits; above a minute the seconds are rounded
#'   to the nearest second. The unit is chosen on the rounded value, so that
#'   rounding never carries a field past its limit.
#'
#' @seealso [print.distrib_fit()], the only caller; [fit_distrib()], which
#'   accumulates the figure over every starting value and every fallback.
#' @keywords internal
fit_format_elapsed <- function(sec) {
  if (!length(sec) || !is.finite(sec)) return(NA_character_)
  # each unit is chosen on the value as it will be printed, so that rounding
  # cannot produce "1e+03 ms", "60 s" or "1 min 60 s"
  if (signif(sec * 1000, 3) < 1000) return(sprintf("%.3g ms", sec * 1000))
  if (signif(sec, 3) < 60) return(sprintf("%.3g s", sec))
  s <- round(sec)
  sprintf("%d min %02d s", as.integer(s %/% 60), as.integer(s %% 60))
}

#' @title Print Method for Maximum-Likelihood Fits
#'
#' @name print.distrib_fit
#'
#' @description
#' Shows a fit in four blocks: the family and the sample size with the
#' log-likelihood and the two information criteria; what the optimizer did;
#' the estimates on the parameter scale; and the same on the link scale, so
#' that the interval that the fit built is visible beside its image. When
#' every link is the identity the last block is replaced by a one-line note.
#'
#' The optimizer line names the method that produced the estimates, which is
#' not always the one requested, and the convergence line names the stopping
#' rule that ended the run. Without that rule `converged` does not identify
#' the test: it records that some test was met. A run that did not converge
#' also prints the score per observation at the point where it stopped, and
#' this number indicates whether the point is usable.
#'
#' For a multivariate fit the coordinates of the covariance structure are
#' replaced by the quantities in which the model is written (the location, the
#' standard deviations and the correlations that [mv_summary()] derives),
#' because a log-Cholesky coordinate has no direct interpretation. Any
#' parameter that the structure does not cover, such as the degrees of
#' freedom of a \eqn{t}, is
#' printed after them.
#'
#' @param x A [distrib_fit()] object.
#' @param digits Number of digits to print. Defaults to 4. The header figures
#'   (log-likelihood, AIC, BIC) are formatted to this many significant digits,
#'   and the tables are rounded to this many decimal places by
#'   [base::round()].
#' @param information `NULL` (the default), an [observed_info()] object or an
#'   [expected_info()] object, which selects the information matrix whose
#'   inverse gives the standard errors and the intervals shown; see
#'   [vcov.distrib_fit()]. `NULL` selects the default of the family (see
#'   [information_objects]). A line under the header names the matrix whenever
#'   it is not the observed information.
#' @param ... Unused, accepted for compatibility with [base::print()].
#'
#' @return `x`, invisibly. The function is called for the output that it writes.
#'
#' @examples
#' set.seed(1)
#' d <- gaussian1_distrib()
#' y <- distrib_rng(d, 500, list(mu = 2, sigma = 3))
#' print(fit_distrib(d, y))
#'
#' # Standard errors from the expected information instead.
#' print(fit_distrib(d, y), information = expected_info())
#'
#' # The interval is built in the link-scale block. On a Gamma written
#' # in the mean and the variance both parameters carry a log link, so both
#' # lower limits are positive by construction.
#' g <- fit_distrib(gamma2_distrib(), rgamma(300, shape = 4, rate = 2))
#' print(g, digits = 3)
#'
#' @seealso [fit_distrib()] for the object;
#'   [confint.distrib_fit()] to recompute the intervals at another level;
#'   [plot.distrib_fit()] to compare the fit with the data;
#'   [mv_summary()] for the multivariate table.
S7::method(print, distrib_fit) <- function(x, digits = 4,
                                           information = NULL, ...) {
  information <- resolve_information(x@distrib, information)
  rep_x <- fit_report(x, information)
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
  if (S7::S7_inherits(information, ExpectedInfo)) {
    cat("Standard errors from the expected information",
        if (is.null(information@method)) "" else
          paste0(" (", sub("^Expectation by ", "", method_label(information@method)), ")"),
        ".\n", sep = "")
  }
  cat("\n")

  # --- the estimates -------------------------------------------------------
  # For a multivariate family the parameters are coordinates of a covariance
  # structure, and nobody reads a Cholesky coordinate. What goes here is what
  # the model is about, which mv_summary() assembles from mv_derived(): the
  # location, and whatever the structure's matrix decomposes into. A structure
  # with fewer free values produces fewer quantities, so the block is as small
  # as the model is.
  if (mv) {
    tab_mv <- tryCatch(mv_summary(x, information = information),
                       error = function(e) NULL)
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
          Estimate = unname(v), `Std. Error` = unname(rep_x$se[seq_len(p)]),
          unname(rep_x$ci[seq_len(p), , drop = FALSE])
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
                     `Std. Error` = rep_x$se[extra],
                     rep_x$ci[extra, , drop = FALSE])
      colnames(tab_x) <- c("Estimate", "Std. Error", lo, hi)
      cat("Other parameters:\n")
      print(round(tab_x, digits))
      cat("\n")
    }
  } else {
    tab <- cbind(Estimate = x@coefficients, `Std. Error` = rep_x$se, rep_x$ci)
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
    tab_eta <- cbind(Estimate = x@eta, `Std. Error` = rep_x$se_eta, rep_x$ci_eta)
    colnames(tab_eta) <- c("Estimate", "Std. Error", lo, hi)
    cat("Link scale:\n")
    print(round(tab_eta, digits))
  }

  invisible(x)
}

#' @title Estimates From a Maximum-Likelihood Fit
#'
#' @name coef.distrib_fit
#'
#' @description
#' Returns the maximum likelihood estimates, on either of the two scales that
#' the fit carries. The default is the parameter scale, \eqn{\hat\theta}, the
#' scale on which the family is interpreted; `scale = "link"` gives
#' \eqn{\hat\eta = g(\hat\theta)}, the point that the optimizer reached.
#'
#' The fit does not recompute either vector. Both were stored at the optimum,
#' and each is the image of the other under the family's links, so
#' `coef(fit, "link")` is `g(coef(fit))` component by component.
#'
#' @param object A [distrib_fit()] object.
#' @param scale Either `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]; anything else signals an error.
#' @param ... Unused, accepted for compatibility with [stats::coef()].
#'
#' @return A named numeric vector of length `length(object@distrib@params)`,
#'   named and ordered as `object@distrib@params`.
#'
#' @examples
#' set.seed(1)
#' d <- gaussian1_distrib()
#' fit <- fit_distrib(d, distrib_rng(d, 400, list(mu = 1, sigma = 2)))
#'
#' coef(fit)
#' coef(fit, scale = "link")
#'
#' # The scale carries a log link and the location the identity, so the link
#' # scale is the logarithm of one estimate and the other unchanged.
#' all.equal(coef(fit, "link")[["sigma"]], log(coef(fit)[["sigma"]]))
#' all.equal(coef(fit, "link")[["mu"]], coef(fit)[["mu"]])
#'
#' @seealso [vcov.distrib_fit()] and [confint.distrib_fit()], which take the
#'   same `scale`; [fit_distrib()] for the fit itself.
#' @importFrom stats coef
S7::method(coef, distrib_fit) <- function(object, scale = c("parameter", "link"), ...) {
  scale <- match.arg(scale)
  if (scale == "link") object@eta else object@coefficients
}

#' @title Variance-Covariance Matrix of a Maximum-Likelihood Fit
#'
#' @name vcov.distrib_fit
#'
#' @description
#' Returns the estimated variance matrix of the estimates, on either scale.
#' The matrix is computed on the link scale as the inverse of an information
#' matrix at \eqn{\hat\eta}. The parameter-scale matrix is its image under the
#' delta method,
#' \deqn{\widehat{\mathrm{Var}}(\hat\theta) = J\,\widehat{\mathrm{Var}}(\hat\eta)\,J,
#'       \qquad J = \mathrm{diag}\!\left(\frac{dg^{-1}}{d\eta}\Big|_{\hat\eta}\right),}
#' and the Jacobian is diagonal because each parameter carries its own scalar
#' link.
#'
#' Every entry is `NA` when the information could not be evaluated or inverted
#' at the optimum. The estimates stand in that case and only the uncertainty
#' is missing.
#'
#' @details
#' # Which information
#'
#' The `information` argument selects the matrix that is inverted, after the
#' fit and independently of the `method` that computed the estimates. It is
#' `NULL`, [observed_info()] or [expected_info()] (see [information_objects]):
#'
#' - `NULL`, the default: [observed_info()] when every parameter of the family
#'   has [params_order()] of at least 2, and [expected_info()] with
#'   `method = NULL` otherwise, for the whole matrix. A family with a kink in a
#'   parameter (the Laplace, the elastic net) has a singular or undefined
#'   observed information in that parameter, and the closed-form expected
#'   information gives finite standard errors.
#' - [observed_info()]: the negative Hessian of the log-likelihood at the
#'   data. It exists for every family, and it describes the curvature of the
#'   likelihood in the sample at hand.
#' - [expected_info()]: the Fisher information. With `method = NULL` the
#'   closed form of the family is used where there is one, and [expect_opg()]
#'   otherwise. An object in `method` is used for every family, one of
#'   [expect_opg()], [expect_integrate()] and [expect_mc()] (see
#'   [expected_derivative_methods]).
#'
#' The two matrices agree asymptotically and differ in a sample. The expected
#' information is averaged over the model, and the observed information is the
#' curvature of the likelihood in the sample at hand. For a family with no
#' closed-form expected information, [expect_opg()] is the outer product of the
#' observed scores, which costs one gradient and is a noisier estimate of the
#' information than the observed Hessian; `expect_integrate(identity = "score")`
#' evaluates the expectation itself, a sum over the support for a discrete
#' family and a quadrature for a continuous one, and is dearer.
#'
#' The fit stores the matrix of the default information and returns it as
#' stored. Any other information is recomputed at the optimum, so it costs one
#' evaluation of that matrix and no refit.
#'
#' @param object A [distrib_fit()] object.
#' @param scale Either `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]; anything else signals an error.
#' @param information `NULL` (the default), [observed_info()] or
#'   [expected_info()], the information matrix to invert; `NULL` selects the
#'   default of the family described above. [expected_info()] carries the
#'   `method` that obtains the expectation. Any other value signals an error.
#'   See [information_objects].
#' @param ... Unused, accepted for compatibility with [stats::vcov()].
#'
#' @return A symmetric numeric matrix of dimension
#'   `length(object@distrib@params)`, with both dimnames set to the parameter
#'   names. With the default `information` its diagonal is the square of
#'   `object@se` on the parameter scale and of `object@se_eta` on the link
#'   scale.
#'
#' @examples
#' set.seed(1)
#' d <- gaussian1_distrib()
#' fit <- fit_distrib(d, distrib_rng(d, 400, list(mu = 1, sigma = 2)))
#'
#' vcov(fit)
#' sqrt(diag(vcov(fit)))          # the standard errors the fit reports
#' all.equal(sqrt(diag(vcov(fit))), fit@se)
#'
#' # The parameter-scale matrix is the link-scale one under the delta method.
#' J <- diag(c(1, coef(fit)[["sigma"]]))   # d theta / d eta: identity, then exp
#' all.equal(unname(vcov(fit)), J %*% vcov(fit, "link") %*% J)
#'
#' # The location and the scale of a Gaussian are orthogonal, so the
#' # off-diagonal entry is zero rather than merely small.
#' vcov(fit)[["mu", "sigma"]]
#'
#' # The two informations are different estimators of one quantity: close on
#' # a well-specified fit of this size, and not identical.
#' sqrt(diag(vcov(fit, information = observed_info())))
#' sqrt(diag(vcov(fit, information = expected_info())))
#'
#' @seealso [coef.distrib_fit()] and [confint.distrib_fit()], which take the
#'   same `scale`; [information_objects] for the two information objects;
#'   [distrib_expected_hessian()] for the information itself, and
#'   [expected_by_opg()] for what `method = expect_opg()` computes.
#' @importFrom stats vcov
S7::method(vcov, distrib_fit) <- function(object, scale = c("parameter", "link"),
                                          information = NULL, ...) {
  scale <- match.arg(scale)
  information <- resolve_information(object@distrib, information)
  if (same_information(information, default_information(object@distrib))) {
    return(if (scale == "link") object@vcov_eta else object@vcov)
  }

  # Recomputed at the optimum. Everything needed is on the fit -- the family,
  # the sample and the estimates -- so no refit is involved and the cost is
  # one evaluation of the information.
  distrib <- object@distrib
  params <- distrib@params
  p <- length(params)
  theta_hat <- object@coefficients
  V_eta <- tryCatch({
    exp_info <- S7::S7_inherits(information, ExpectedInfo)
    I_eta <- -fit_hess_matrix(distrib, object@y, as.list(theta_hat),
                              expected = exp_info,
                              method = if (exp_info) information@method)
    solve(I_eta)
  }, error = function(e) matrix(NA_real_, p, p))
  dimnames(V_eta) <- list(params, params)
  if (scale == "link") return(V_eta)

  J <- fit_dtheta_deta(distrib, object@eta)
  V_theta <- diag(J, p, p) %*% V_eta %*% diag(J, p, p)
  dimnames(V_theta) <- list(params, params)
  V_theta
}

#' @title Confidence Intervals for a Maximum-Likelihood Fit
#'
#' @name confint.distrib_fit
#'
#' @description
#' Returns Wald intervals for the estimated parameters. They are built
#' symmetrically on the link scale, \eqn{\hat\eta \pm z_{1-\alpha/2}\,
#' \mathrm{se}(\hat\eta)}, and mapped through \eqn{g^{-1}} when the parameter
#' scale is requested, so a limit cannot leave the parameter's domain: a scale
#' has a positive lower limit and a probability stays inside \eqn{(0,1)}. The
#' two ends are sorted after mapping, because a link need not be increasing.
#'
#' Any level can be computed from the stored estimate and standard error
#' without refitting, so an interval at 99% follows from a fit computed at 95%
#' at the cost of one quantile.
#'
#' @param object A [distrib_fit()] object.
#' @param parm Which parameters to report, by name or by position. If missing,
#'   which is the default, all parameters are reported. A name that is not a
#'   parameter of this fit, or a position outside the range, signals an error
#'   naming the argument.
#' @param level Confidence level, a single number in \eqn{(0, 1)}. Defaults to
#'   the level at which the fit was computed, `object@level`.
#' @param scale Either `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]; anything else signals an error. It is the fourth
#'   argument, so name it: `confint(fit, "sigma", "link")` passes `"link"`
#'   as `level`, which then signals an error.
#' @param information `NULL` (the default), an [observed_info()] object or an
#'   [expected_info()] object, which selects the information matrix whose
#'   inverse gives the standard errors that the interval is built on; see
#'   [vcov.distrib_fit()].
#' @param ... Unused, accepted for compatibility with [stats::confint()].
#'
#' @return A numeric matrix with one row per requested parameter and two
#'   columns, named for the two tail probabilities as percentages
#'   (`"2.5%"` and `"97.5%"` at the default level). Row names are the
#'   parameter names. Both entries are `NA` where the standard error is.
#'
#' @examples
#' set.seed(1)
#' d <- gaussian1_distrib()
#' fit <- fit_distrib(d, distrib_rng(d, 400, list(mu = 1, sigma = 2)))
#'
#' confint(fit)
#' confint(fit, level = 0.99)     # no refit: wider, from the same estimates
#' confint(fit, "sigma")
#'
#' # The interval is built on the link scale and mapped back, so the lower
#' # limit of a scale is exp() of a real number and cannot be negative.
#' confint(fit, "sigma", scale = "link")
#' lo_eta <- confint(fit, "sigma", scale = "link")[1]
#' all.equal(confint(fit, "sigma")[1], exp(lo_eta), check.attributes = FALSE)
#'
#' # A probability near the boundary is fitted from 47 successes in 50 trials.
#' fb <- fit_distrib(bernoulli_distrib(), rep(0:1, c(3, 47)))
#' confint(fb)                    # inside (0, 1) by construction
#'
#' @seealso [coef.distrib_fit()] and [vcov.distrib_fit()], which take the same
#'   `scale`; [fit_distrib()], whose `level` sets the default here.
#' @importFrom stats confint qnorm
S7::method(confint, distrib_fit) <- function(object, parm, level = object@level,
                                             scale = c("parameter", "link"),
                                             information = NULL, ...) {
  scale <- match.arg(scale)
  se_eta <- fit_report(object, information)$se_eta
  params <- object@distrib@params
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      level <= 0 || level >= 1) {
    stop("'level' must be a single number in (0, 1).", call. = FALSE)
  }

  z <- stats::qnorm(1 - (1 - level) / 2)
  ci_eta <- cbind(lower = object@eta - z * se_eta,
                  upper = object@eta + z * se_eta)
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

#' @title Log-Likelihood of a Maximum-Likelihood Fit
#'
#' @name logLik.distrib_fit
#'
#' @description
#' Returns the maximized log-likelihood as a `logLik` object, so that
#' [stats::AIC()], [stats::BIC()] and anything else in \pkg{stats} that reads
#' one can be applied to a fit. The value is summed over observations and is
#' not divided by \eqn{n}, whatever scaling the optimizer used.
#'
#' The `df` attribute is the number of estimated parameters, which for these
#' fits is every parameter of the family, because a fit estimates all of them
#' and holds none fixed. The `nobs` attribute is `object@n`, the row count for a
#' multivariate response and the length otherwise.
#'
#' @param object A [distrib_fit()] object.
#' @param ... Unused, accepted for compatibility with [stats::logLik()].
#'
#' @return An object of class `logLik`: a single number carrying the attributes
#'   `df` (the parameter count) and `nobs` (the observation count).
#'
#' @examples
#' set.seed(1)
#' d <- gaussian1_distrib()
#' y <- distrib_rng(d, 400, list(mu = 1, sigma = 2))
#' fit <- fit_distrib(d, y)
#'
#' logLik(fit)
#' c(df = attr(logLik(fit), "df"), nobs = attr(logLik(fit), "nobs"))
#'
#' # The criteria that the fit reports are the ones that stats builds from this.
#' all.equal(AIC(logLik(fit)), fit@aic)
#' all.equal(BIC(logLik(fit)), fit@bic)
#'
#' # Two families are compared on the same data. The Student t has one more
#' # parameter, and the AIC shows whether the tails are worth it.
#' ft <- fit_distrib(student_t1_distrib(), y)
#' c(gaussian = AIC(logLik(fit)), student_t = AIC(logLik(ft)))
#'
#' @seealso [fit_distrib()] for the fit; [stats::AIC()] and [stats::BIC()],
#'   which read this.
#' @importFrom stats logLik
S7::method(logLik, distrib_fit) <- function(object, ...) {
  structure(object@loglik,
            df = length(object@coefficients),
            nobs = object@n,
            class = "logLik")
}

#' @title Simulate From a Fitted Distribution
#'
#' @name simulate.distrib_fit
#'
#' @description
#' Draws new samples from the fitted distribution, evaluated at the maximum
#' likelihood estimates and ignoring their uncertainty. Each replicate has the
#' same length as the data from which the fit was computed, so a replicate is
#' directly comparable with the observations: this is the draw that a
#' parametric bootstrap and a posterior-predictive style check both need.
#'
#' The draws come from [distrib_rng()], so a family with no closed-form
#' generator is sampled by the package's own fallback and the cost is that
#' fallback's.
#'
#' @param object A [distrib_fit()] object.
#' @param nsim Number of replicates to draw, a single positive integer.
#'   Defaults to 1.
#' @param seed Optional seed. When supplied it initializes the generator, and
#'   the `.Random.seed` in effect before the call is restored afterwards, so
#'   simulating does not disturb the caller's stream. When `NULL`, the default,
#'   the caller's stream is used and advanced. Either way the seed actually in
#'   force is attached to the result as its `"seed"` attribute, so a run can
#'   be reproduced after the fact.
#' @param ... Unused, accepted for compatibility with [stats::simulate()].
#'
#' @return For a univariate fit, a data frame of `object@n` rows and `nsim`
#'   columns named `sim_1` to `sim_nsim`, in which each column is one
#'   replicate. For a
#'   multivariate fit each replicate is a matrix, and the data frame has
#'   `object@n` rows and `nsim * p` columns named `sim_1.v1`, `sim_1.v2`, and
#'   so on. The `"seed"` attribute carries the generator state described
#'   above.
#'
#' @examples
#' set.seed(1)
#' y <- rnorm(200, 3, 2)
#' fit <- fit_distrib(gaussian1_distrib(), y)
#'
#' sims <- simulate(fit, 20, seed = 42)
#' dim(sims)
#'
#' # The replicates give a parametric bootstrap of any statistic, here the
#' # median.
#' quantile(vapply(sims, median, numeric(1)), c(0.025, 0.975))
#'
#' # The seed argument leaves the caller's stream where it found it.
#' set.seed(7); before <- runif(1)
#' set.seed(7); invisible(simulate(fit, 2, seed = 42)); after <- runif(1)
#' all.equal(before, after)
#'
#' @seealso [fit_distrib()] for the fit; [distrib_rng()] for the generator;
#'   [plot.distrib_fit()] to compare the fit with the data from which it was
#'   computed.
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

#' @title Plot a Fitted Distribution Against the Data
#'
#' @name plot.distrib_fit
#'
#' @description
#' Compares the fitted distribution with the sample from which it was
#' estimated, choosing the comparison from the family's support. A continuous
#' family is drawn as the fitted density over a kernel estimate of the data,
#' with a rug of the observations underneath. A discrete one is drawn as the
#' observed relative frequencies in bars with the fitted mass overlaid,
#' because a kernel estimate misrepresents a sample from a discrete
#' distribution.
#'
#' A multivariate fit is drawn as a panel matrix: the fitted marginal density
#' and a kernel estimate of the data on the diagonal, the fitted contours over
#' the observations below it, the fitted correlation above.
#'
#' @param x A [distrib_fit()] object.
#' @param n_grid Number of points at which the fitted density is evaluated, a
#'   single positive integer. Defaults to 512. Read for a continuous family
#'   only; a discrete one is evaluated on its support.
#' @param rug Logical of length 1: draw a rug of the observations. Defaults to
#'   `TRUE` when there are at most 2000 of them, above which the rug becomes a
#'   solid band that conveys no information.
#' @param legend Logical of length 1: add a legend. Defaults to `TRUE`.
#' @param col_fit,col_data Colors of the fitted curve and of the empirical
#'   summary, in any form that [grDevices::col2rgb()] accepts.
#' @param mv_which For a multivariate fit, which coordinates to show, by index.
#'   Defaults to all of them. At most three coordinates can be drawn, because
#'   a larger panel matrix is not readable, so a fit with more than three
#'   coordinates signals an error unless `mv_which` selects at most three.
#'   The argument is ignored for a univariate fit.
#' @param ... Further arguments passed to [graphics::plot()], such as `main`,
#'   `xlab` or `xlim`.
#'
#' @return `x`, invisibly. The function is called for the plot that it draws.
#'
#' @examples
#' set.seed(1)
#' y <- rgamma(300, shape = 4, rate = 2)
#' fit <- fit_distrib(gamma2_distrib(), y)
#' plot(fit)
#'
#' # A count family is compared on its support, not through a kernel estimate.
#' fp <- fit_distrib(poisson_distrib(), rpois(300, 4))
#' plot(fp, main = "Poisson fit")
#'
#' @seealso [fit_distrib()] for the fit; [simulate.distrib_fit()] for draws
#'   from it; [plot.multivariate_distrib()], the same panel matrix drawn from
#'   a distribution with no data beside it.
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
      which_arg = "mv_which",
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
    lab <- c("kernel density", paste("fitted", d@distrib_name))
    if (isTRUE(rug)) {
      graphics::legend("topright", bty = "n",
        col = c(col_data, col_fit, col_data),
        lwd = c(2, 2, NA), pch = c(NA, NA, 124),
        legend = c(lab, "observations"))
    } else {
      graphics::legend("topright", bty = "n", lwd = 2,
        col = c(col_data, col_fit), legend = lab)
    }
  }

  invisible(x)
}
