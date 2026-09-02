#' @include multivariate.R mvgaussian_distrib.R mvstudent_t_distrib.R
NULL

#' @title Interpretable Quantities of a Multivariate Distribution
#'
#' @description
#' Returns the quantities a reader of a multivariate fit wants: standard
#' deviations, correlations, and whatever else the family's matrix
#' parametrization is about. It returns with them the Jacobian that carries
#' standard errors onto those quantities, and the scale each one's confidence
#' interval should be built on.
#'
#' @details
#' # Why coordinates are not quantities
#'
#' The free values of a \pkg{parameters7} parametrization are coordinates
#' chosen so that an optimizer can move freely. The logarithm of the second
#' diagonal entry of a Cholesky factor has an estimate and a standard error,
#' and neither answers a question anybody asked. This generic names the
#' quantities that do and supplies \eqn{\partial g/\partial\theta}, so
#' [mv_summary()] can apply the delta method to them.
#'
#' # The interval scale each quantity declares
#'
#' Each quantity carries the scale its interval should be built on, exactly as
#' [fit_distrib()] builds a univariate interval on the link scale and maps it
#' back. A standard deviation is intervalled on the log scale, so the interval
#' cannot reach zero; a correlation on Fisher's
#' \eqn{z = \mathrm{artanh}(\rho)}, so it cannot leave \eqn{(-1, 1)}; an
#' unconstrained quantity on its own scale. On the raw scale a correlation's
#' interval routinely exceeds one.
#'
#' # What the default method reports
#'
#' The method registered on [multivariate_distrib()] returns the distinct
#' entries of the matrix [mv_sigma()] produces, named `sigma_v1_v1` and so on,
#' with a Jacobian from one central difference per parameter. A family whose
#' matrix is not a covariance still reports it on its own scale, which beats
#' reporting a Cholesky coordinate. The two shipped families override with
#' closed-form Jacobians.
#'
#' @param distrib An object inheriting from [multivariate_distrib()].
#' @param theta A named list or vector of parameters, each component a single
#'   number. Aligned by the generic before dispatch, so any order is accepted.
#' @param ... Passed to methods. No shipped method reads it.
#'
#' @return A named list with
#'   \describe{
#'     \item{`value`}{a named numeric vector of the quantities;}
#'     \item{`jacobian`}{a numeric matrix with one row per quantity and one
#'       column per parameter of `distrib`, in `distrib@params` order;}
#'     \item{`transform`}{a character vector, one of `"identity"`, `"log"`,
#'       `"atanh"` or `"logit"` per quantity, naming the scale its interval is
#'       built on;}
#'     \item{`block`}{a character vector labeling the group each quantity
#'       belongs to, by which the printed summary is laid out.}
#'   }
#'   All four are named after the quantities and are the same length.
#'
#' @section Notation:
#' \eqn{\theta} is the parameter vector, \eqn{g} the map to the reported
#' quantities, \eqn{\Sigma} the matrix the family carries and \eqn{\rho} a
#' correlation.
#'
#' @seealso [mv_summary()], which applies the delta method to this,
#'   [mv_derived.MvGaussianDistrib()] and [mv_derived.MvStudentTDistrib()] for
#'   the two closed-form methods, and [mv_sigma()] for the matrix.
#'
#' @examples
#' d <- mvgaussian1_distrib(2)
#' theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
#'               sigma_L2.1 = 0.5)
#' der <- mv_derived(d, theta)
#' der$value
#'
#' # Each quantity declares the scale its interval is built on.
#' der$transform
#' der$block
#'
#' # The Jacobian is exact for this family, against a numerical one.
#' g <- function(v) {
#'   t2 <- as.list(v); names(t2) <- d@params
#'   mv_derived(d, t2)$value
#' }
#' max(abs(der$jacobian - numDeriv::jacobian(g, unlist(theta))))
#'
#' # The mean parameters do not enter any of them, so those columns are zero.
#' der$jacobian[, c("mu1", "mu2")]
#'
#' # A structured matrix reports its own quantities as a further block.
#' a <- mvgaussian1_distrib(3, parameters7::ar1(3))
#' mv_derived(a, as.list(stats::setNames(c(0, 0, 0, 0.1, 0.3), a@params)))$value
#'
#' @export
mv_derived <- S7::new_generic("mv_derived", "distrib",
  function(distrib, theta, ...) {
    theta <- align_theta(distrib, theta)
    S7::S7_dispatch()
  }
)


#' @title Distinct Entries of a Symmetric Matrix, and Their Labels
#'
#' @description
#' Returns the row and column indices of the lower triangle of a
#' \eqn{p \times p} matrix, diagonal included, together with a label per entry.
#' A symmetric matrix has `p * (p + 1) / 2` distinct entries, and reporting all
#' \eqn{p^2} would give every off-diagonal quantity twice.
#'
#' @param p The side of the matrix, a single positive whole number.
#' @param prefix The label prefix, a single string. Entry \eqn{(i, j)} is named
#'   `prefix_vi_vj`, so `prefix = "cor"` gives `cor_v2_v1`.
#'
#' @return A named list with `i` and `j`, integer vectors of length
#'   `p * (p + 1) / 2` holding the row and column of each entry in
#'   column-major order over the lower triangle, and `name`, the character
#'   vector of labels.
#'
#' @seealso [mv_sd_cor()] and [mv_derived.multivariate_distrib()], which label
#'   their quantities with this.
#'
#' @examples
#' distributions7:::mv_entry_index(3, "cor")
#'
#' # The count is p(p + 1) / 2, the distinct entries of a symmetric matrix.
#' c(got = length(distributions7:::mv_entry_index(4, "sigma")$i),
#'   expected = 4 * 5 / 2)
#'
#' # Every index pair is on or below the diagonal.
#' idx <- distributions7:::mv_entry_index(4, "sigma")
#' all(idx$i >= idx$j)
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


#' @title The Quantities the Matrix Parametrization Is About
#'
#' @description
#' Returns the block a \pkg{parameters7} family declares through
#' `param_readable()`, widened to the distribution's own parameter vector: the
#' parametrization's Jacobian columns are placed in the stretch of
#' `distrib@params` its free values occupy, and every other column is zero. An
#' AR(1) covariance, for instance, is about a scale and a correlation, and
#' those are what a reader wants beside the standard deviations.
#'
#' @details
#' A family that declares nothing returns `NULL`, and the summary is then
#' whatever [mv_derived()] produced without it. The label says `(precision)`
#' when the distribution is inverted, because the parametrization does not know
#' which side of the model it was handed to.
#'
#' A multivariate family written OUTSIDE this package may have no `param`
#' property at all, so the property is asked for with `S7::prop_names()` rather
#' than assumed.
#'
#' @param distrib A [multivariate_distrib()] object.
#' @param theta A named list of parameters, already aligned.
#'
#' @return `NULL` where the parametrization declares nothing, or a named list
#'   in the shape of [mv_derived()]'s return: `value`, `jacobian`, `transform`
#'   and `block`, with the Jacobian as wide as `distrib@params`.
#'
#' @seealso [parameters7::param_readable()] for what a parametrization
#'   declares, [mv_append_block()] for how the result is joined on, and
#'   [mv_derived.MvGaussianDistrib()] for the consumer.
#'
#' @examples
#' # A log-Cholesky covariance declares nothing: its free values are
#' # coordinates and nothing more.
#' d <- mvgaussian1_distrib(2)
#' theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
#'               sigma_L2.1 = 0.5)
#' th <- distributions7:::align_theta(d, theta)
#' is.null(distributions7:::mv_param_block(d, th))
#'
#' # An AR(1) covariance is about a scale and a correlation, and says so.
#' a <- mvgaussian1_distrib(3, parameters7::ar1(3))
#' th <- as.list(stats::setNames(c(0, 0, 0, 0.1, 0.3), a@params))
#' pb <- distributions7:::mv_param_block(a, distributions7:::align_theta(a, th))
#' pb$value
#' pb$transform
#' pb$block
#'
#' # The Jacobian is as wide as the distribution's parameter vector, with the
#' # mean columns zero.
#' dim(pb$jacobian)
#' pb$jacobian[, 1:3]
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


#' @title Append One Block of Derived Quantities to Another
#'
#' @description
#' Concatenates two results in the shape of [mv_derived()], stacking the
#' `value`, `transform` and `block` vectors and row-binding the two Jacobians.
#' It is what joins a parametrization's own quantities, from
#' [mv_param_block()], onto the standard deviations and correlations.
#'
#' @param out A named list with `value`, `jacobian`, `transform` and `block`,
#'   as [mv_derived()] returns.
#' @param extra A second such list, or `NULL`, in which case `out` is returned
#'   unchanged. Its Jacobian must have the same number of columns as `out`'s.
#'
#' @return A named list of the same four components, with the rows of `extra`
#'   after those of `out`.
#'
#' @seealso [mv_param_block()] for the usual second argument and
#'   [mv_derived.MvGaussianDistrib()] for the caller.
#'
#' @examples
#' d <- mvgaussian1_distrib(2)
#' theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
#'               sigma_L2.1 = 0.5)
#' der <- mv_derived(d, theta)
#'
#' # A NULL second argument leaves the first alone.
#' identical(distributions7:::mv_append_block(der, NULL), der)
#'
#' # Otherwise the rows stack and the Jacobian keeps its width.
#' both <- distributions7:::mv_append_block(der, der)
#' c(rows = length(both$value), cols = ncol(both$jacobian))
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
#'
#' @description
#' Returns the distinct entries of the matrix [mv_sigma()] produces, named
#' `sigma_vi_vj` after the coordinates they belong to, with a Jacobian from one
#' central difference in each parameter. This is the fallback for a family that
#' registers no method of its own: reporting the matrix on its own scale is
#' worth more to a reader than reporting a Cholesky coordinate, even where the
#' entries are not the standard deviations and correlations a closed-form
#' method would give.
#'
#' @details
#' The diagonal entries are variances and are intervalled on the log scale; the
#' off-diagonal ones are unconstrained given the diagonal and are intervalled on
#' their own. That is weaker than the closed-form methods, whose correlations
#' ride Fisher's \eqn{z} and so cannot leave \eqn{(-1, 1)}.
#'
#' @param distrib An object inheriting from [multivariate_distrib()] that has
#'   no method of its own.
#' @param theta A named list of parameters, already aligned by the generic.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list with `value`, `jacobian`, `transform` and `block`, as
#'   [mv_derived()] documents. There are `p * (p + 1) / 2` quantities, blocked
#'   as `"Variances"` and `"Covariances"`.
#'
#' @seealso [mv_derived.MvGaussianDistrib()] and
#'   [mv_derived.MvStudentTDistrib()] for the two closed-form methods,
#'   [mv_entry_index()] for the labeling, and [mv_derived()] for the generic.
#'
#' @examples
#' d <- mvgaussian1_distrib(2)
#' theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
#'               sigma_L2.1 = 0.5)
#'
#' # The gaussian overrides, so reach the default method directly.
#' base <- S7::method(mv_derived, multivariate_distrib)
#' d0 <- base(d, distributions7:::align_theta(d, theta))
#' d0$value
#' d0$transform
#'
#' # The entries really are the matrix, read off its lower triangle.
#' mv_sigma(d, theta)
#'
#' # Against the closed-form method, which reports standard deviations and a
#' # correlation instead of variances and a covariance.
#' mv_derived(d, theta)$value
#'
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


#' @title Standard Deviations and Correlations of a Structured Matrix
#'
#' @description
#' Turns a symmetric positive definite matrix and the derivatives of its
#' entries into the standard deviations and correlations a reader wants, with
#' the Jacobian written out. Writing \eqn{\Sigma = DRD} with \eqn{D} the
#' diagonal of standard deviations,
#' \deqn{\frac{\partial s_j}{\partial\theta_k} = \frac{A_k[j,j]}{2 s_j},
#'   \qquad
#'   \frac{\partial \rho_{jl}}{\partial\theta_k}
#'     = \frac{A_k[j,l]}{s_j s_l}
#'       - \frac{\rho_{jl}}{2}\left(\frac{A_k[j,j]}{\Sigma_{jj}}
#'         + \frac{A_k[l,l]}{\Sigma_{ll}}\right),}
#' with \eqn{A_k = \partial\Sigma/\partial\theta_k}. Nothing new is computed:
#' the arrays come from the parametrization the caller already has.
#'
#' @param sigma A \eqn{p \times p} symmetric positive definite numeric matrix.
#' @param a A list of \eqn{p \times p} numeric matrices, one per parameter of
#'   the distribution and in `params` order, holding
#'   \eqn{\partial\Sigma/\partial\theta_k}. A parameter the matrix does not
#'   depend on contributes a matrix of zeros.
#' @param params The distribution's parameter names, which become the column
#'   names of the Jacobian.
#' @param sd_label The prefix for the diagonal quantities. Defaults to `"sd"`;
#'   the Student t passes `"scale_sd"`, its diagonal quantities not being
#'   standard deviations of the response.
#' @param cor_label The prefix for the off-diagonal quantities. Defaults to
#'   `"cor"`.
#' @param sd_block,cor_block The block labels the printed summary groups by.
#'
#' @return A named list with `value`, `jacobian`, `transform` and `block`, as
#'   [mv_derived()] documents. There are \eqn{p} diagonal quantities on the log
#'   scale and \eqn{p(p-1)/2} off-diagonal ones on Fisher's \eqn{z}.
#'
#' @section Notation:
#' \eqn{\Sigma} is the matrix, \eqn{D} the diagonal of its square roots,
#' \eqn{R} the correlation matrix, \eqn{s_j} a standard deviation,
#' \eqn{\rho_{jl}} a correlation, \eqn{\theta} the parameter vector and
#' \eqn{A_k = \partial\Sigma/\partial\theta_k}.
#'
#' @seealso [mv_sigma_derivs()] for the arrays it takes,
#'   [mv_derived.MvGaussianDistrib()] and [mv_derived.MvStudentTDistrib()] for
#'   the two callers, and [mv_entry_index()] for the labeling.
#'
#' @examples
#' d <- mvgaussian1_distrib(2)
#' theta <- distributions7:::align_theta(
#'   d, list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
#'           sigma_L2.1 = 0.5))
#' S <- mv_sigma(d, theta)
#' a <- distributions7:::mv_sigma_derivs(d, theta, 2)
#'
#' sc <- distributions7:::mv_sd_cor(S, a, d@params)
#' sc$value
#' sc$transform
#'
#' # The values really are the square roots of the diagonal and the
#' # correlations off it.
#' c(sqrt(diag(S)), S[2, 1] / sqrt(S[1, 1] * S[2, 2]))
#'
#' # The Jacobian against a numerical one taken through mv_sigma().
#' g <- function(v) {
#'   t2 <- as.list(v); names(t2) <- d@params
#'   M <- mv_sigma(d, t2)
#'   c(sqrt(diag(M)), M[2, 1] / sqrt(M[1, 1] * M[2, 2]))
#' }
#' max(abs(sc$jacobian - numDeriv::jacobian(g, unlist(theta))))
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


#' @title Derivatives of the Covariance with Respect to Every Parameter
#'
#' @description
#' Returns \eqn{\partial\Sigma/\partial\theta_k} for EVERY parameter of a
#' multivariate distribution, in `distrib@params` order, so that the list can
#' be indexed by parameter position without any bookkeeping at the call site.
#' The parametrization supplies the derivatives of its own free values; the
#' location components and, for a Student \eqn{t}, the degrees of freedom leave
#' the matrix alone and contribute a matrix of zeros.
#'
#' @details
#' Where the parametrization carries the PRECISION the arrays are carried onto
#' the covariance by \eqn{\partial\Sigma/\partial\eta_k = -\Sigma A_k \Sigma},
#' so the result is always in the covariance whichever side the model is
#' written on.
#'
#' @param distrib A [multivariate_distrib()] object with a `param` property.
#' @param theta A named list of parameters, already aligned.
#' @param n_before The number of parameters before the matrix parametrization's
#'   free values, which is \eqn{p} for both shipped families. Those positions
#'   get zero matrices.
#'
#' @return A list of \eqn{p \times p} numeric matrices, as long as
#'   `distrib@params`. Every entry outside the parametrization's stretch is a
#'   matrix of zeros.
#'
#' @section Notation:
#' \eqn{\Sigma} is the covariance or scale matrix, \eqn{\eta} the free vector of
#' the matrix parametrization, \eqn{\theta} the distribution's parameter vector
#' and \eqn{A_k} the parametrization's own derivative array.
#'
#' @seealso [mv_sd_cor()], the consumer, and [mvg_pieces()], which does the
#'   same conversion for the gaussian's own methods.
#'
#' @examples
#' d <- mvgaussian1_distrib(2)
#' theta <- distributions7:::align_theta(
#'   d, list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
#'           sigma_L2.1 = 0.5))
#' a <- distributions7:::mv_sigma_derivs(d, theta, 2)
#'
#' # One entry per parameter, and the two mean entries are zero.
#' length(a) == d@n_params
#' a[[1]]
#'
#' # The third entry is the derivative in the first free value, against a
#' # difference of mv_sigma().
#' h <- 1e-5
#' tp <- theta; tp$sigma_log_L1 <- tp$sigma_log_L1 + h
#' tm <- theta; tm$sigma_log_L1 <- tm$sigma_log_L1 - h
#' max(abs(a[[3]] - (mv_sigma(d, tp) - mv_sigma(d, tm)) / (2 * h)))
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
#'
#' @description
#' Returns the standard deviations and correlations of the response, whichever
#' side the parametrization carries, with the closed-form Jacobian
#' [mv_sd_cor()] supplies. A PRECISION parametrization reports two further
#' blocks, which are what it describes directly: the conditional variances
#' \eqn{1/\Omega_{jj} = \operatorname{Var}(Y_j \mid Y_{-j})}, and above two
#' dimensions the partial correlations
#' \eqn{-\Omega_{jk}/\sqrt{\Omega_{jj}\Omega_{kk}}}, the correlation of two
#' coordinates given all the others, which is zero exactly where the precision
#' has a zero.
#'
#' @details
#' # What a precision's diagonal means
#'
#' The quantity with a reading is the conditional VARIANCE, so the diagonal
#' quantities [mv_sd_cor()] produces from \eqn{\Omega} are square roots of the
#' wrong thing; they are dropped and \eqn{1/\Omega_{jj}} is reported instead.
#' Its ratio to the marginal variance is \eqn{1 - R_j^2} for the regression of
#' that coordinate on all the others.
#'
#' At \eqn{p = 2} there is nothing to condition on, so the partial correlation
#' IS the correlation and is not printed twice.
#'
#' # The parametrization's own quantities
#'
#' Whatever the matrix parametrization declares through
#' `parameters7::param_readable()` is appended as a further block. An AR(1)
#' covariance is about a scale and a correlation; a log-Cholesky one declares
#' nothing and the summary stops at the standard deviations.
#'
#' @param distrib An [MvGaussianDistrib] object, from [mvgaussian1_distrib()].
#' @param theta A named list of parameters, already aligned by the generic.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list with `value`, `jacobian`, `transform` and `block`, as
#'   [mv_derived()] documents. A covariance parametrization gives
#'   \eqn{p(p+1)/2} quantities; a precision one adds \eqn{p} conditional
#'   variances and, above \eqn{p = 2}, \eqn{p(p-1)/2} partial correlations.
#'
#' @section Notation:
#' \eqn{\Sigma} is the covariance, \eqn{\Omega = \Sigma^{-1}} the precision,
#' \eqn{p} the dimension, \eqn{R_j^2} the coefficient of determination of the
#' regression of coordinate \eqn{j} on the others, and \eqn{Y_{-j}} the
#' response with that coordinate removed.
#'
#' @seealso [mv_sd_cor()] for the closed-form Jacobian,
#'   [mv_param_block()] for the appended block, [mv_summary()] for the printed
#'   result, and [mv_derived()] for the generic.
#'
#' @examples
#' d <- mvgaussian1_distrib(2)
#' theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
#'               sigma_L2.1 = 0.5)
#' mv_derived(d, theta)$value
#'
#' # The precision side reports the same law's standard deviations and
#' # correlation, and adds the conditional variances.
#' o <- mvgaussian2_distrib(2, parameters7::log_cholesky(2))
#' th_o <- list(mu1 = 0, mu2 = 0, omega_log_L1 = 0, omega_log_L2 = 0,
#'              omega_L2.1 = 0.5)
#' od <- mv_derived(o, th_o)
#' od$value
#' od$block
#'
#' # A conditional variance is 1 / Omega_jj, and is at most the marginal one.
#' Om <- parameters7::param_value(o@param, unlist(th_o)[3:5])
#' c(conditional = 1 / Om[1, 1], marginal = mv_sigma(o, th_o)[1, 1])
#'
#' # At three dimensions the partial correlations appear as a block of their
#' # own, the partial and the marginal no longer coinciding.
#' o3 <- mvgaussian2_distrib(3, parameters7::log_cholesky(3))
#' th3 <- as.list(stats::setNames(
#'   c(0, 0, 0, 0, 0, 0, 0.5, -0.4, 0.3), o3@params))
#' unique(mv_derived(o3, th3)$block)
#'
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
#'
#' @description
#' Returns the square roots of the diagonal of the SCALE matrix and the
#' correlations it implies, with the closed-form Jacobian [mv_sd_cor()]
#' supplies. The correlations are the response's as well: the covariance is
#' \eqn{\nu\Sigma/(\nu-2)} where it exists, and a positive multiple of a matrix
#' leaves its correlations alone. The diagonal quantities are NOT standard
#' deviations of the response, and are named `scale_sd_v1`, ..., `scale_sd_vp`
#' and blocked as `"Scale standard deviations"` to say so.
#'
#' @details
#' The degrees of freedom appear in no quantity here. \eqn{\nu} is already an
#' interpretable parameter on its own scale, so [confint.distrib_fit()]
#' reports it and this generic does not repeat it. That is also why the
#' \eqn{\nu} column of the Jacobian is zero.
#'
#' @param distrib An [MvStudentTDistrib] object, from
#'   [mvstudent_t1_distrib()].
#' @param theta A named list of parameters, already aligned by the generic.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list with `value`, `jacobian`, `transform` and `block`, as
#'   [mv_derived()] documents: \eqn{p} scale standard deviations on the log
#'   scale and \eqn{p(p-1)/2} correlations on Fisher's \eqn{z}, plus whatever
#'   the matrix parametrization declares.
#'
#' @section Notation:
#' \eqn{\Sigma} is the scale matrix, \eqn{\nu} the degrees of freedom and
#' \eqn{p} the dimension.
#'
#' @seealso [mv_sd_cor()] for the closed-form Jacobian,
#'   [mv_derived.MvGaussianDistrib()], whose diagonal quantities are genuine
#'   standard deviations, [variance.MvStudentTDistrib()] for the covariance,
#'   and [mv_derived()] for the generic.
#'
#' @examples
#' d <- mvstudent_t1_distrib(2)
#' theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
#'               sigma_L2.1 = 0.5, nu = 6)
#' der <- mv_derived(d, theta)
#' der$value
#' der$block
#'
#' # The diagonal quantities are the scale matrix's, not the response's: the
#' # response is more spread than they say by sqrt(nu / (nu - 2)).
#' c(scale_sd = der$value[["scale_sd_v1"]],
#'   response_sd = sqrt(variance(d, theta)[1, 1]))
#'
#' # The correlation is the same either way, the two matrices being positive
#' # multiples of each other.
#' cv <- function(m) m[1, 2] / sqrt(m[1, 1] * m[2, 2])
#' c(reported = der$value[["cor_v1_v2"]],
#'   from_covariance = cv(variance(d, theta)))
#'
#' # nu enters no quantity here, so its Jacobian column is zero.
#' der$jacobian[, "nu"]
#'
#' @keywords internal
S7::method(mv_derived, MvStudentTDistrib) <- function(distrib, theta, ...) {
  p <- distrib@n_dim
  sigma <- unname(mv_sigma(distrib, theta))
  a <- mv_sigma_derivs(distrib, theta, n_before = p)
  out <- mv_sd_cor(sigma, a, distrib@params,
                   sd_label = "scale_sd",
                   sd_block = "Scale standard deviations")

  if (!isTRUE(distrib@inverted)) {
    return(mv_append_block(out, mv_param_block(distrib, theta)))
  }

  # The inverse parametrization's own reading. Sigma^-1 and its derivatives are
  # the matrix parameter's matrix directly, so they are taken from it.
  #
  # A partial correlation is minus the correlation of the inverse, and here it
  # IS the response's: the response's precision is (nu - 2) / nu times this
  # matrix, and a positive multiple cancels out of -M_jk / sqrt(M_jj M_kk).
  # The diagonal is the case that does NOT carry over. For a gaussian
  # 1 / Omega_jj is the conditional variance; here it is the Schur complement
  # of the SCALE matrix, so it is a conditional SCALE, and the conditional law
  # of a t is a t whose scale carries a further factor depending on the
  # conditioning values. It is named cscale to say so, as the diagonal
  # quantities above are named scale_sd.
  s <- distrib@param
  v <- mv_flat_theta(distrib, theta)
  eta <- v[p + seq_len(s@n_free)]
  om <- unname(parameters7::param_value(s, eta))
  aw <- vector("list", distrib@n_params)
  aw[p + seq_len(s@n_free)] <- lapply(parameters7::param_d1(s, eta), unname)

  pc <- mv_sd_cor(om, aw, distrib@params)
  # In two dimensions there is nothing to condition on, so the partial
  # correlation IS the correlation and printing it again would be noise.
  keep <- grepl("^cor_", names(pc$value)) & p >= 3L
  cscale <- 1 / diag(om)
  nm_cs <- sprintf("cscale_v%d", seq_len(p))
  jac_cs <- matrix(0, p, distrib@n_params,
                   dimnames = list(nm_cs, distrib@params))
  for (k in seq_along(aw)) {
    if (is.null(aw[[k]])) next
    jac_cs[, k] <- -diag(aw[[k]]) / diag(om)^2
  }

  # The partial correlations are named by renaming the correlations they come
  # from, and every one of the four fields is renamed with them: the consumer
  # indexes them positionally, so a stale name there is inert today and is a
  # trap for whoever indexes by name next.
  nm_pc <- sub("^cor_", "pcor_", names(pc$value)[keep])
  extra <- list(
    value = c(stats::setNames(cscale, nm_cs),
              stats::setNames(-pc$value[keep], nm_pc)),
    jacobian = rbind(jac_cs,
                     `rownames<-`(-pc$jacobian[keep, , drop = FALSE], nm_pc)),
    transform = c(stats::setNames(rep("log", p), nm_cs),
                  stats::setNames(unname(pc$transform[keep]), nm_pc)),
    block = c(stats::setNames(rep("Conditional scales", p), nm_cs),
              stats::setNames(rep("Partial correlations", sum(keep)), nm_pc))
  )

  mv_append_block(mv_append_block(out, extra), mv_param_block(distrib, theta))
}


#' @title Interpretable Estimates of a Multivariate Fit
#'
#' @description
#' Reports the standard deviations, correlations and whatever else
#' [mv_derived()] declares for the fitted distribution, each with its standard
#' error and its confidence interval. A multivariate fit estimates the free
#' values of a \pkg{parameters7} parametrization, and those are coordinates:
#' the estimate and standard error of `sigma_log_L2` answer no question anybody
#' asked. This function carries the fit's variance matrix onto the quantities
#' that do.
#'
#' @details
#' # The delta method
#'
#' \deqn{\widehat{\mathrm{Var}}\{g(\hat\theta)\}
#'   = J \, \widehat{\mathrm{Var}}(\hat\theta) \, J^\top,
#'   \qquad J = \partial g/\partial\theta,}
#' with \eqn{J} taken from [mv_derived()], in closed form for the two families
#' that ship and by one central difference for any other.
#'
#' # Where each interval is built
#'
#' On the scale the quantity declares, then mapped back, which is the
#' discipline [fit_distrib()] applies to a univariate parameter. A standard
#' deviation is intervalled on the log scale and a correlation on Fisher's
#' \eqn{z}, so neither interval can leave the set its quantity lives in. An
#' interval on the raw scale routinely puts a correlation above one.
#'
#' The standard error is carried onto that scale by the transform's own
#' derivative, so the reported `Std. Error` column stays on the quantity's
#' natural scale while the limits are computed on the transformed one.
#'
#' @param object A [distrib_fit()] of a MULTIVARIATE distribution. A univariate
#'   fit is rejected: there the parameters are already the interpretable
#'   quantities and [confint.distrib_fit()] is the function to call.
#' @param level The confidence level, a single number strictly inside
#'   \eqn{(0, 1)}. Defaults to the fit's own `level`. Anything else is an error.
#'
#' @return A data frame with one row per quantity and four columns:
#'   `Estimate`, `Std. Error`, and the two confidence limits, named after the
#'   percentiles they are (`2.5%` and `97.5%` at the default level). Its row
#'   names are the quantity names, and it carries the attribute `"block"`, a
#'   character vector naming the group each row belongs to, by which `print()`
#'   lays the summary out.
#'
#' @section Notation:
#' \eqn{\theta} is the parameter vector, \eqn{g} the map to the reported
#' quantities, \eqn{J} its Jacobian and \eqn{z} Fisher's transform of a
#' correlation.
#'
#' @seealso [mv_derived()] for the quantities and the Jacobian,
#'   [confint.distrib_fit()] for the coordinates themselves, and
#'   [fit_distrib()] for the fit.
#'
#' @examples
#' set.seed(1)
#' d <- mvgaussian1_distrib(2)
#' truth <- list(mu1 = 0, mu2 = 1, sigma_log_L1 = 0,
#'               sigma_log_L2 = 0, sigma_L2.1 = 0.7)
#' y <- distrib_rng(d, 500, truth)
#' fit <- fit_distrib(d, y)
#'
#' # The standard deviations and the correlation, which is what one reads.
#' mv_summary(fit)
#'
#' # The coordinates the fit was estimated on, beside them.
#' confint(fit)
#'
#' # Every interval stays inside the set its quantity lives in: the standard
#' # deviations are positive and the correlation is inside (-1, 1).
#' s <- mv_summary(fit)
#' c(sd_lower_positive = all(s[1:2, 3] > 0),
#'   cor_inside = s[3, 3] > -1 && s[3, 4] < 1)
#'
#' # A wider level widens the interval and leaves the estimate alone.
#' mv_summary(fit, level = 0.99)
#'
#' # A univariate fit is rejected by name.
#' try(mv_summary(fit_distrib(gaussian1_distrib(), rnorm(50))))
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
