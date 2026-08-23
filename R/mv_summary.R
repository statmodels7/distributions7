#' @include multivariate.R mvgaussian_distrib.R mvstudent_t_distrib.R
NULL

#' Interpretable Quantities of a Multivariate Distribution
#'
#' @description
#' Returns the quantities a reader of a multivariate fit actually wants ---
#' standard deviations, correlations, and whatever else the family's matrix
#' parameter means --- together with the Jacobian needed to carry standard
#' errors onto them.
#'
#' @details
#' The free values of a \pkg{parameters7} structure are coordinates chosen so
#' that an optimizer can move freely; they are not quantities anyone reads.
#' The logarithm of the third diagonal entry of a Cholesky factor has an
#' estimate and a standard error, and neither answers a question. This generic
#' names the quantities that do, and supplies
#' \eqn{\partial g/\partial\theta} so that [mv_summary()] can apply
#' the delta method to them.
#'
#' Each quantity also declares the scale its confidence interval should be
#' built on, exactly as [fit_distrib()] builds a univariate interval
#' on the link scale and maps it back: a standard deviation is intervalled on
#' the log scale so that the interval cannot reach zero, a correlation on
#' Fisher's \eqn{z = \mathrm{artanh}(\rho)} so that it cannot leave
#' \eqn{(-1, 1)}, and an unconstrained quantity on its own scale.
#'
#' The default method, registered on [multivariate_distrib()],
#' returns the distinct entries of the matrix [mv_sigma()] produces,
#' named after the coordinates they belong to, with a Jacobian obtained by one
#' central difference. A family whose matrix is not a covariance therefore
#' still reports something on its original scale, which is better than
#' reporting a Cholesky coordinate.
#'
#' @param distrib An object inheriting from class
#'   [multivariate_distrib()].
#' @param theta A named list or vector of parameters.
#' @param ... Passed to methods.
#'
#' @return A list with
#'   \describe{
#'     \item{`value`}{a named numeric vector of the quantities;}
#'     \item{`jacobian`}{a matrix with one row per quantity and one column
#'       per parameter of `distrib`;}
#'     \item{`transform`}{a character vector, one of `"identity"`,
#'       `"log"` or `"atanh"` per quantity, naming the scale its
#'       interval is built on;}
#'     \item{`block`}{a character vector labeling the group each quantity
#'       belongs to, used to lay the printed summary out.}
#'   }
#'
#' @seealso [mv_summary()], [mv_sigma()]
#'
#' @examples
#' d <- mvgaussian_distrib(2)
#' theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
#'               sigma_L2.1 = 0.5)
#' mv_derived(d, theta)$value
#'
#' @export
mv_derived <- S7::new_generic("mv_derived", "distrib",
  function(distrib, theta, ...) {
    theta <- align_theta(distrib, theta)
    S7::S7_dispatch()
  }
)


#' Distinct Entries of a Symmetric Matrix, and Their Labels
#'
#' @description
#' Returns the row and column indices of the lower triangle including the
#' diagonal, together with a label for each.
#'
#' @param p The side of the matrix.
#' @param prefix The string the labels start with.
#'
#' @return A list with `i`, `j` and `name`.
#'
#' @keywords internal
mv_entry_index <- function(p, prefix) {
  ij <- which(lower.tri(matrix(0, p, p), diag = TRUE), arr.ind = TRUE)
  ij <- ij[order(ij[, "col"], ij[, "row"]), , drop = FALSE]
  list(
    i = ij[, "row"], j = ij[, "col"],
    name = sprintf("%s_v%d_v%d", prefix, ij[, "row"], ij[, "col"])
  )
}


#' The Quantities the Matrix Parameter Is About
#'
#' @description
#' The block a \pkg{parameters7} family declares through
#' [parameters7::param_readable()], with its Jacobian widened from
#' the free vector to the whole parameter vector of the distribution.
#'
#' @details
#' The free values of the structure occupy a contiguous stretch of
#' `distrib@params`, after the means and before anything the family adds
#' of its own, so widening the Jacobian is placing its columns in that stretch
#' and leaving the rest at zero: the quantities depend on no other parameter.
#' A family that declares nothing yields `NULL` and the summary is what
#' it was.
#'
#' @param distrib A [multivariate_distrib()] object.
#' @param theta A named list of parameters.
#'
#' @return A list in the shape of [mv_derived()], or `NULL`.
#'
#' @keywords internal
mv_param_block <- function(distrib, theta) {
  # A family written outside the package inherits from multivariate_distrib
  # directly and need not be built on a parameters7 structure at all, so the
  # property is asked for rather than assumed.
  if (!"param" %in% S7::prop_names(distrib)) return(NULL)
  s <- distrib@param
  p <- distrib@n_dim
  v <- mv_flat_theta(distrib, theta)
  r <- parameters7::param_readable(s, v[p + seq_len(s@n_free)])
  if (is.null(r)) return(NULL)

  nm <- names(r$value)
  jac <- matrix(0, length(nm), distrib@n_params,
                dimnames = list(nm, distrib@params))
  jac[, p + seq_len(s@n_free)] <- r$jacobian
  # Which matrix the family describes is the distribution's business, not the
  # structure's: the same structure carries a covariance on one side of a
  # model and a precision on the other.
  label <- if (isTRUE(distrib@inverted)) {
    paste(r$label, "(precision)")
  } else {
    r$label
  }
  list(value = r$value, jacobian = jac, transform = r$transform,
       block = stats::setNames(rep(label, length(nm)), nm))
}


#' Append One Block of Derived Quantities to Another
#'
#' @description
#' Concatenates two lists in the shape of [mv_derived()], returning
#' the first unchanged when the second is `NULL`.
#'
#' @param out A list as described in [mv_derived()].
#' @param extra A list of the same shape, or `NULL`.
#'
#' @return A list as described in [mv_derived()].
#'
#' @seealso [mv_param_block()]
#'
#' @keywords internal
mv_append_block <- function(out, extra) {
  if (is.null(extra)) return(out)
  list(
    value = c(out$value, extra$value),
    jacobian = rbind(out$jacobian, extra$jacobian),
    transform = c(out$transform, extra$transform),
    block = c(out$block, extra$block)
  )
}


#' @title Matrix Entries as the Default Interpretable Quantities
#' @name mv_derived.multivariate_distrib
#' @description
#' The distinct entries of the matrix [mv_sigma()] returns, with a
#' Jacobian from one central difference in each parameter. This is what a
#' family gets when it says nothing more specific: the matrix on its own scale,
#' named after the coordinates, rather than the matrix parameter's coordinates.
#' @param distrib A [multivariate_distrib()] object.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A list as described in [mv_derived()].
#' @keywords internal
S7::method(mv_derived, multivariate_distrib) <- function(distrib, theta, ...) {
  p <- distrib@n_dim
  v <- mv_flat_theta(distrib, theta)
  idx <- mv_entry_index(p, "sigma")

  at <- function(x) {
    m <- mv_sigma(distrib, as.list(stats::setNames(x, distrib@params)))
    m[cbind(idx$i, idx$j)]
  }
  val <- at(v)
  jac <- matrix(0, length(val), length(v),
                dimnames = list(idx$name, distrib@params))
  for (k in seq_along(v)) {
    h <- 1e-5 * max(1, abs(v[k]))
    up <- dn <- v
    up[k] <- v[k] + h
    dn[k] <- v[k] - h
    jac[, k] <- (at(up) - at(dn)) / (2 * h)
  }
  mv_append_block(list(
    value = stats::setNames(val, idx$name),
    jacobian = jac,
    transform = stats::setNames(
      ifelse(idx$i == idx$j, "log", "identity"), idx$name
    ),
    block = stats::setNames(
      ifelse(idx$i == idx$j, "Variances", "Covariances"), idx$name
    )
  ), mv_param_block(distrib, theta))
}


#' Standard Deviations and Correlations of a Structured Matrix
#'
#' @description
#' Turns a symmetric positive definite matrix and the derivatives of its
#' entries into the standard deviations and correlations it decomposes into,
#' together with their Jacobian.
#'
#' @details
#' Writing \eqn{\Sigma = D R D} with \eqn{D} the diagonal of standard
#' deviations, the two readings are \eqn{s_j = \sqrt{\Sigma_{jj}}} and
#' \eqn{\rho_{jk} = \Sigma_{jk}/(s_j s_k)}. Differentiating,
#' \deqn{\frac{\partial s_j}{\partial \eta_k} = \frac{A_k[j,j]}{2 s_j},
#'       \qquad
#'       \frac{\partial \rho_{jk}}{\partial \eta_l}
#'         = \frac{A_l[j,k]}{s_j s_k}
#'           - \frac{\rho_{jk}}{2}
#'             \left(\frac{A_l[j,j]}{\Sigma_{jj}}
#'                 + \frac{A_l[k,k]}{\Sigma_{kk}}\right),}
#' with \eqn{A_k = \partial\Sigma/\partial\eta_k}.
#'
#' @param sigma The matrix.
#' @param a A list of derivative matrices, one per parameter, in the order of
#'   the distribution's parameters. Entries may be `NULL` for parameters
#'   the matrix does not depend on.
#' @param params The parameter names, used to label the Jacobian's columns.
#' @param sd_label The label the diagonal quantities are named with.
#' @param cor_label The label the off-diagonal quantities are named with.
#' @param sd_block,cor_block The headings the two groups print under.
#'
#' @return A list as described in [mv_derived()].
#'
#' @keywords internal
mv_sd_cor <- function(sigma, a, params, sd_label = "sd", cor_label = "cor",
                      sd_block = "Standard deviations",
                      cor_block = "Correlations") {
  p <- ncol(sigma)
  s <- sqrt(diag(sigma))
  nm_sd <- sprintf("%s_v%d", sd_label, seq_len(p))

  pairs <- if (p > 1L) {
    utils::combn(p, 2L)
  } else {
    matrix(integer(0), nrow = 2L, ncol = 0L)
  }
  nm_cor <- if (ncol(pairs)) {
    sprintf("%s_v%d_v%d", cor_label, pairs[1, ], pairs[2, ])
  } else {
    character(0)
  }
  rho <- if (ncol(pairs)) {
    vapply(seq_len(ncol(pairs)), function(m) {
      j <- pairs[1, m]
      k <- pairs[2, m]
      sigma[j, k] / (s[j] * s[k])
    }, numeric(1))
  } else {
    numeric(0)
  }

  nm <- c(nm_sd, nm_cor)
  jac <- matrix(0, length(nm), length(params), dimnames = list(nm, params))
  for (l in seq_along(params)) {
    al <- a[[l]]
    if (is.null(al)) next
    jac[seq_len(p), l] <- diag(al) / (2 * s)
    if (ncol(pairs)) {
      jac[p + seq_len(ncol(pairs)), l] <- vapply(
        seq_len(ncol(pairs)), function(m) {
          j <- pairs[1, m]
          k <- pairs[2, m]
          al[j, k] / (s[j] * s[k]) -
            (rho[m] / 2) * (al[j, j] / sigma[j, j] + al[k, k] / sigma[k, k])
        }, numeric(1)
      )
    }
  }

  list(
    value = stats::setNames(c(s, rho), nm),
    jacobian = jac,
    transform = stats::setNames(
      c(rep("log", p), rep("atanh", length(rho))), nm
    ),
    block = stats::setNames(
      c(rep(sd_block, p), rep(cor_block, length(rho))), nm
    )
  )
}


#' Derivatives of the Covariance with Respect to Every Parameter
#'
#' @description
#' Returns \eqn{\partial\Sigma/\partial\theta_k} for each parameter of a
#' multivariate distribution built on a \pkg{parameters7} structure, as a list
#' aligned with `distrib@params` and `NULL` where the covariance does
#' not depend on the parameter.
#'
#' @details
#' The mean components and, for a Student \eqn{t}, the degrees of freedom leave
#' the matrix alone, so those entries are `NULL` and cost nothing. When
#' the matrix parameter parametrizes the precision the chain rule of an inverse
#' applies, \eqn{\partial\Sigma/\partial\eta_k = -\Sigma A_k \Sigma}.
#'
#' @param distrib A distribution carrying a `param` property.
#' @param theta A named list of parameters, already aligned.
#' @param n_before How many parameters precede the matrix parameter's free values.
#'
#' @return A list of matrices and `NULL`s, of length
#'   `distrib@n_params`.
#'
#' @keywords internal
mv_sigma_derivs <- function(distrib, theta, n_before) {
  s <- distrib@param
  v <- mv_flat_theta(distrib, theta)
  eta <- v[n_before + seq_len(s@n_free)]
  d <- parameters7::param_d1(s, eta)
  inverted <- isTRUE(S7::prop_exists(distrib, "inverted")) &&
    isTRUE(distrib@inverted)
  if (inverted) {
    sigma <- mv_sigma(distrib, theta)
    d <- lapply(d, function(ak) -(sigma %*% ak %*% sigma))
  }
  out <- vector("list", distrib@n_params)
  out[n_before + seq_len(s@n_free)] <- lapply(d, unname)
  out
}


#' @title Standard Deviations and Correlations of a Multivariate Gaussian
#' @name mv_derived.MvGaussianDistrib
#' @description
#' The standard deviations and correlations of the response, whichever side the
#' structure parametrizes. A precision structure additionally reports the
#' **partial** correlations, which are what it describes directly:
#' \eqn{-\Omega_{jk}/\sqrt{\Omega_{jj}\Omega_{kk}}} is the correlation of two
#' coordinates given all the others, and it is zero exactly where the precision
#' has a zero.
#' @param distrib A [MvGaussianDistrib()] object.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A list as described in [mv_derived()].
#' @keywords internal
S7::method(mv_derived, MvGaussianDistrib) <- function(distrib, theta, ...) {
  p <- distrib@n_dim
  sigma <- unname(mv_sigma(distrib, theta))
  a <- mv_sigma_derivs(distrib, theta, n_before = p)
  out <- mv_sd_cor(sigma, a, distrib@params)

  if (!isTRUE(distrib@inverted)) {
    return(mv_append_block(out, mv_param_block(distrib, theta)))
  }

  # The precision's own reading. Omega and its derivatives are the matrix parameter's
  # matrix directly, with no inversion, so they are taken from it rather than
  # from the covariance computed above.
  s <- distrib@param
  v <- mv_flat_theta(distrib, theta)
  eta <- v[p + seq_len(s@n_free)]
  omega <- unname(parameters7::param_value(s, eta))
  aw <- vector("list", distrib@n_params)
  aw[p + seq_len(s@n_free)] <- lapply(
    parameters7::param_d1(s, eta), unname
  )

  # A partial correlation is minus the correlation of the precision, so its
  # derivatives are those of mv_sd_cor() with the sign flipped. The diagonal
  # quantities that helper produces are square roots and are dropped: what a
  # precision's diagonal means is the CONDITIONAL VARIANCE
  # \eqn{1/\Omega_{jj} = \mathrm{Var}(Y_j \mid Y_{-j})}, whose ratio to the
  # marginal variance is \eqn{1 - R_j^2} for the regression of that coordinate
  # on all the others. The variance is the quantity with that reading, not its
  # square root.
  pc <- mv_sd_cor(omega, aw, distrib@params)
  # In two dimensions there is nothing to condition on, so the partial
  # correlation IS the correlation and printing it again would be noise.
  keep <- grepl("^cor_", names(pc$value)) & p >= 3L
  cvar <- 1 / diag(omega)
  nm_cvar <- sprintf("cvar_v%d", seq_len(p))
  jac_cvar <- matrix(0, p, distrib@n_params,
                     dimnames = list(nm_cvar, distrib@params))
  for (l in seq_along(aw)) {
    if (is.null(aw[[l]])) next
    # d(1/w_jj)/d eta = -A[j,j] / w_jj^2
    jac_cvar[, l] <- -diag(aw[[l]]) / diag(omega)^2
  }

  mv_append_block(list(
    value = c(out$value,
              stats::setNames(cvar, nm_cvar),
              stats::setNames(-pc$value[keep],
                              sub("^cor_", "pcor_", names(pc$value)[keep]))),
    jacobian = rbind(
      out$jacobian,
      jac_cvar,
      `rownames<-`(-pc$jacobian[keep, , drop = FALSE],
                   sub("^cor_", "pcor_", names(pc$value)[keep]))
    ),
    transform = c(out$transform,
                  stats::setNames(rep("log", p), nm_cvar),
                  stats::setNames(rep("atanh", sum(keep)),
                                  sub("^cor_", "pcor_", names(pc$value)[keep]))),
    block = c(out$block,
              stats::setNames(rep("Conditional variances", p), nm_cvar),
              stats::setNames(rep("Partial correlations", sum(keep)),
                              sub("^cor_", "pcor_", names(pc$value)[keep])))
  ), mv_param_block(distrib, theta))
}


#' @title Scale Standard Deviations and Correlations of a Multivariate t
#' @name mv_derived.MvStudentTDistrib
#' @description
#' The square roots of the diagonal of the **scale** matrix and the
#' correlations it implies. The correlations are those of the response as well,
#' since the covariance is \eqn{\nu\Sigma/(\nu-2)} and a positive multiple of a
#' matrix leaves its correlations alone; the diagonal quantities are not
#' standard deviations of the response and are named to say so.
#' @param distrib A [MvStudentTDistrib()] object.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A list as described in [mv_derived()].
#' @keywords internal
S7::method(mv_derived, MvStudentTDistrib) <- function(distrib, theta, ...) {
  p <- distrib@n_dim
  sigma <- unname(mv_sigma(distrib, theta))
  a <- mv_sigma_derivs(distrib, theta, n_before = p)
  mv_append_block(
    mv_sd_cor(sigma, a, distrib@params,
              sd_label = "scale_sd", sd_block = "Scale standard deviations"),
    mv_param_block(distrib, theta)
  )
}


#' Interpretable Estimates of a Multivariate Fit
#'
#' @description
#' Reports the standard deviations, correlations and whatever else
#' [mv_derived()] declares for the fitted distribution, each with its
#' standard error and confidence interval.
#'
#' @details
#' A multivariate fit estimates the free values of a \pkg{parameters7}
#' structure, and those are coordinates rather than quantities: the estimate
#' and standard error of `sigma_log_L2` answer no question anybody asked.
#' This function carries the fit's variance matrix onto the quantities that do,
#' by the delta method,
#' \deqn{\widehat{\mathrm{Var}}\{g(\hat\theta)\}
#'   = J \, \widehat{\mathrm{Var}}(\hat\theta) \, J^\top,
#'   \qquad J = \partial g/\partial\theta,}
#' with \eqn{J} taken from [mv_derived()], in closed form for the
#' families that ship with the package.
#'
#' Each interval is built on the scale the quantity declares and mapped back,
#' which is the same discipline [fit_distrib()] applies to a
#' univariate parameter: a standard deviation is intervalled on the log scale
#' and a correlation on Fisher's \eqn{z}, so neither interval can leave the set
#' its quantity lives in. An interval on the raw scale would routinely put a
#' correlation above one.
#'
#' @param object A [distrib_fit()] of a multivariate distribution.
#' @param level The confidence level. Defaults to the fit's own.
#'
#' @return A data frame with one row per quantity and the columns
#'   `Estimate`, `Std. Error` and the two confidence limits, carrying
#'   the attribute `"block"` that names the group each row belongs to.
#'
#' @seealso [mv_derived()], [confint.distrib_fit()]
#'
#' @examples
#' set.seed(1)
#' d <- mvgaussian_distrib(2)
#' y <- distrib_rng(d, 500, list(mu1 = 0, mu2 = 1, sigma_log_L1 = 0,
#'                               sigma_log_L2 = 0, sigma_L2.1 = 0.7))
#' fit <- fit_distrib(d, y)
#'
#' # the standard deviations and the correlation, which is what one reads
#' mv_summary(fit)
#'
#' @export
mv_summary <- function(object, level = object@level) {
  if (!S7::S7_inherits(object, distrib_fit)) {
    stop("'object' must be a distrib_fit.", call. = FALSE)
  }
  d <- object@distrib
  if (!S7::S7_inherits(d, multivariate_distrib)) {
    stop(paste0(
      "mv_summary() is for a multivariate fit. For a univariate one the\n",
      "  parameters are already the interpretable quantities; use confint()."
    ), call. = FALSE)
  }
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    stop("'level' must be a single number in (0, 1).", call. = FALSE)
  }

  der <- mv_derived(d, as.list(object@coefficients))
  v <- object@vcov
  se <- sqrt(pmax(diag(der$jacobian %*% v %*% t(der$jacobian)), 0))

  z <- stats::qnorm(1 - (1 - level) / 2)
  est <- der$value
  fwd <- list(identity = identity, log = log, atanh = atanh,
              logit = stats::qlogis)
  inv <- list(identity = identity, log = exp, atanh = tanh,
              logit = stats::plogis)

  lo <- hi <- numeric(length(est))
  for (i in seq_along(est)) {
    tr <- der$transform[[i]]
    # The derivative of the transform is what carries the standard error onto
    # the scale the interval is built on; the interval is then mapped back.
    d_tr <- switch(tr,
      identity = 1,
      log = 1 / est[i],
      atanh = 1 / (1 - est[i]^2),
      logit = 1 / (est[i] * (1 - est[i]))
    )
    e_tr <- fwd[[tr]](est[i])
    s_tr <- se[i] * abs(d_tr)
    lo[i] <- inv[[tr]](e_tr - z * s_tr)
    hi[i] <- inv[[tr]](e_tr + z * s_tr)
  }

  out <- data.frame(
    Estimate = as.numeric(est),
    `Std. Error` = as.numeric(se),
    lower = lo,
    upper = hi,
    row.names = names(est),
    check.names = FALSE
  )
  names(out)[3:4] <- paste0(
    format(100 * c((1 - level) / 2, 1 - (1 - level) / 2), trim = TRUE), "%"
  )
  attr(out, "block") <- unname(der$block)
  out
}
