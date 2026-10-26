#' @include multivariate.R mvgaussian_distrib.R mvstudent_t_distrib.R
NULL

#' @title Interpretable Quantities of a Multivariate Distribution
#'
#' @description
#' Returns the quantities that a reader of a multivariate fit reads: standard
#' deviations, correlations, and the quantities that the matrix parametrization
#' of the family declares. It also returns the Jacobian that carries standard
#' errors onto those quantities, and the scale on which the confidence interval
#' of each one is built.
#'
#' @details
#' # Coordinates and quantities
#'
#' The free values of a \pkg{parameters7} parametrization are coordinates chosen
#' so that an optimizer can move freely. The logarithm of the second diagonal
#' entry of a Cholesky factor has an estimate and a standard error, but these do
#' not describe a quantity of interest. This generic names the quantities that
#' do and supplies \eqn{\partial g/\partial\theta}, so that [mv_summary()] can
#' apply the delta method to them.
#'
#' # The interval scale that each quantity declares
#'
#' Each quantity carries the scale on which its interval is built, as
#' [fit_distrib()] builds a univariate interval on the link scale and maps it
#' back. The interval of a standard deviation is built on the log scale, so it
#' cannot reach zero; that of a correlation on Fisher's
#' \eqn{z = \mathrm{artanh}(\rho)}, so it cannot leave \eqn{(-1, 1)}; that of an
#' unconstrained quantity on its own scale.
#' On the raw scale the interval of a correlation can exceed one.
#'
#' # What the default method reports
#'
#' The method registered on [multivariate_distrib()] returns the distinct
#' entries of the matrix that [mv_sigma()] produces, named `sigma_v1_v1` and so
#' on, with a Jacobian from one central difference per parameter, followed by
#' the quantities that the matrix parametrization declares through
#' `param_readable()`. A family whose matrix is not a covariance reports that
#' matrix on its own scale. The multivariate gaussian and Student t families
#' override this method with closed-form Jacobians; the Dirichlet and the
#' multinomial use it, and report the covariance of the response followed by the
#' probabilities of the simplex.
#'
#' @param distrib An object inheriting from [multivariate_distrib()].
#' @param theta A named list or vector of parameters, where each component is a
#'   single number. The generic aligns it before dispatch, so any order is
#'   accepted.
#' @param ... Passed to methods. The methods of the package do not use it.
#'
#' @return A named list with
#'   \describe{
#'     \item{`value`}{a named numeric vector of the quantities;}
#'     \item{`jacobian`}{a numeric matrix with one row per quantity and one
#'       column per parameter of `distrib`, in `distrib@params` order;}
#'     \item{`transform`}{a character vector, one of `"identity"`, `"log"`,
#'       `"atanh"` or `"logit"` per quantity, naming the scale on which its
#'       interval is built;}
#'     \item{`block`}{a character vector labeling the group to which each
#'       quantity belongs and by which the printed summary is laid out.}
#'   }
#'   All four are named after the quantities and are the same length.
#'
#' @section Notation:
#' \eqn{\theta} is the parameter vector, \eqn{g} the map to the reported
#' quantities, \eqn{\Sigma} the matrix that the family carries and \eqn{\rho} a
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
#' # Each quantity declares the scale on which its interval is built.
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


#' @title Quantities Declared by the Matrix Parametrization
#'
#' @description
#' Returns the block that a \pkg{parameters7} parametrization declares through
#' `param_readable()`, widened to the parameter vector of the distribution: the
#' Jacobian columns of the parametrization are placed at the positions of its
#' free values in `distrib@params`, and every other column is zero. An AR(1)
#' covariance, for instance, declares a scale and a correlation, which a
#' summary reports beside the standard deviations.
#'
#' @details
#' A parametrization that declares nothing gives `NULL`, and the summary is
#' then the result of [mv_derived()] without it. The label ends in
#' `(precision)` when the declared quantities describe a precision, because the
#' parametrization itself does not record whether its matrix is a covariance or
#' a precision. They describe a precision when the distribution is inverted and
#' the parametrization is not an inverse one, or the other way round: an inverse
#' parametrization such as `parameters7::ar1_inv()` declares, with
#' `inverse = TRUE`, the quantities of the matrix that it inverts, so under an
#' inverted distribution they describe the covariance.
#'
#' The positions of the free values are found by name: the family names each
#' free value with a prefix followed by the free name, after the location for
#' the gaussian and the Student t and first for the Dirichlet and the
#' multinomial. Where a name does not identify exactly one parameter, the
#' positions after the location are used. A multivariate family written outside
#' this package may have no `param` property, so the function tests for it with
#' `S7::prop_names()` and returns `NULL` when it is absent.
#'
#' @param distrib A [multivariate_distrib()] object.
#' @param theta A named list of parameters, already aligned.
#'
#' @return `NULL` where the parametrization declares nothing, or a named list
#'   in the shape of the return value of [mv_derived()]: `value`, `jacobian`,
#'   `transform` and `block`, with the Jacobian as wide as `distrib@params`,
#'   and a logical attribute `precision`, `TRUE` when the quantities describe
#'   a precision.
#'
#' @seealso [parameters7::param_readable()] for what a parametrization
#'   declares, [mv_append_block()] for how the result is joined on, and
#'   [mv_derived.multivariate_distrib()] and [mv_own_block()] for the
#'   consumers.
#'
#' @examples
#' # A log-Cholesky covariance declares no block.
#' d <- mvgaussian1_distrib(2)
#' theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
#'               sigma_L2.1 = 0.5)
#' th <- distributions7:::align_theta(d, theta)
#' is.null(distributions7:::mv_param_block(d, th))
#'
#' # An AR(1) covariance declares a scale and a correlation.
#' a <- mvgaussian1_distrib(3, parameters7::ar1(3))
#' th <- as.list(stats::setNames(c(0, 0, 0, 0.1, 0.3), a@params))
#' pb <- distributions7:::mv_param_block(a, distributions7:::align_theta(a, th))
#' pb$value
#' pb$transform
#' pb$block
#'
#' # The Jacobian is as wide as the parameter vector of the distribution, with
#' # the mean columns zero.
#' dim(pb$jacobian)
#' pb$jacobian[, 1:3]
#'
#' # The simplex of a Dirichlet declares its probabilities, whose free values
#' # come first among the parameters.
#' dd <- dirichlet_distrib(3)
#' thd <- list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8)
#' pbd <- distributions7:::mv_param_block(dd, distributions7:::align_theta(dd, thd))
#' pbd$value
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
  # Where the structure's free values sit among the parameters: after the
  # location for the elliptical families, first for the simplex-valued ones.
  # Each is found by its name, which the family writes as a prefix and the
  # free name; a name that does not identify one parameter falls back to the
  # position after the location.
  pos <- vapply(s@free_names, function(f) {
    hit <- which(endsWith(distrib@params, paste0("_", f)))
    if (length(hit) == 1L) hit else NA_integer_
  }, integer(1))
  if (anyNA(pos) || anyDuplicated(pos)) pos <- p + seq_len(s@n_free)
  r <- parameters7::param_readable(s, v[pos])
  if (is.null(r)) return(NULL)

  nm <- names(r$value)
  jac <- matrix(0, length(nm), distrib@n_params,
                dimnames = list(nm, distrib@params))
  jac[, pos] <- r$jacobian
  # Which matrix the family carries is the distribution's business, not the
  # structure's: the same structure carries a covariance on one side of a
  # model and a precision on the other. An inverse structure (ar1_inv() and
  # the like) declares the quantities of the matrix it inverts, so the two
  # inversions cancel: ar1_inv() under mvgaussian2_distrib() describes the
  # AR(1) covariance.
  inv <- "inverted" %in% S7::prop_names(distrib) && isTRUE(distrib@inverted)
  inv <- xor(inv, isTRUE(r$inverse))
  label <- if (inv) {
    paste(r$label, "(precision)")
  } else {
    r$label
  }
  structure(list(value = r$value, jacobian = jac, transform = r$transform,
                 block = stats::setNames(rep(label, length(nm)), nm)),
            precision = inv)
}


#' @title Quantities Declared by a Structured Parametrization
#'
#' @description
#' Returns the block that [mv_param_block()] returns, with the common variance
#' `scale` of the structure replaced by its square root where the matrix is a
#' covariance or a scale matrix, or `NULL` where the structure declares no
#' block. A structured matrix such as `compound_symmetry()` or `ar1()` is fixed
#' by these few quantities, so a summary reports them in place of the standard
#' deviation of every coordinate and the correlation of every pair, which would
#' repeat them.
#'
#' @details
#' The standard deviation is \eqn{\sqrt{s}} with \eqn{s} the `scale`, so its
#' Jacobian row is the `scale` row times \eqn{1/(2\sqrt{s})}; its interval is
#' built on the log scale, as that of `scale` was. Where the declared quantities
#' describe a precision (the attribute `precision` of the result of
#' [mv_param_block()]) the `scale` is a precision's common diagonal and is
#' reported as it is.
#'
#' @param distrib A multivariate distribution object.
#' @param theta Its parameters, aligned.
#' @param sd_label The name given to the square root of `scale`.
#'
#' @return The result of [mv_param_block()], possibly with one row renamed and
#'   transformed, or `NULL`.
#'
#' @seealso [mv_param_block()] for the declared block, [mv_derived()] for the
#'   consumer.
#'
#' @examples
#' d <- mvgaussian1_distrib(3, parameters7::compound_symmetry(3))
#' th <- distributions7:::align_theta(d, list(mu1 = 0, mu2 = 0, mu3 = 0,
#'   sigma_log_scale = log(4), sigma_logit_rho = 0))
#' distributions7:::mv_own_block(d, th)$value
#'
#' @keywords internal
mv_own_block <- function(distrib, theta, sd_label = "sd") {
  pb <- mv_param_block(distrib, theta)
  if (is.null(pb) || isTRUE(attr(pb, "precision"))) return(pb)
  k <- which(names(pb$value) == "scale")
  if (length(k) != 1L) return(pb)
  s <- pb$value[[k]]
  pb$jacobian[k, ] <- pb$jacobian[k, ] / (2 * sqrt(s))
  pb$value[[k]] <- sqrt(s)
  pb$transform[[k]] <- "log"
  nm <- names(pb$value)
  nm[k] <- sd_label
  names(pb$value) <- names(pb$transform) <- names(pb$block) <- nm
  rownames(pb$jacobian) <- nm
  pb
}


#' @title Append One Block of Derived Quantities to Another
#'
#' @description
#' Concatenates two results in the shape of [mv_derived()], stacking the
#' `value`, `transform` and `block` vectors and row-binding the two Jacobians.
#' It joins the quantities of a parametrization, from [mv_param_block()], onto
#' the standard deviations and correlations.
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
#' Returns the distinct entries of the matrix that [mv_sigma()] produces, named
#' `sigma_vi_vj` after their coordinates, with a Jacobian from one central
#' difference in each parameter, followed by the block that [mv_param_block()]
#' returns when the parametrization declares one. It is the method of a family
#' that registers none of its own, among them the Dirichlet and the
#' multinomial, whose matrix is the covariance of the response. The entries
#' are variances and covariances, where a closed-form method reports standard
#' deviations and correlations.
#'
#' @details
#' The diagonal entries are variances and their intervals are built on the log
#' scale; the off-diagonal ones are unconstrained given the diagonal and their
#' intervals are built on their own scale. The closed-form methods report
#' correlations on Fisher's \eqn{z}, whose intervals stay inside \eqn{(-1, 1)};
#' the interval of a covariance has no such bound.
#'
#' @param distrib An object inheriting from [multivariate_distrib()] that has
#'   no method of its own.
#' @param theta A named list of parameters, already aligned by the generic.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list with `value`, `jacobian`, `transform` and `block`, as
#'   [mv_derived()] documents. There are `p * (p + 1) / 2` matrix entries,
#'   blocked as `"Variances"` and `"Covariances"`, followed by the quantities
#'   that the matrix parametrization declares through `param_readable()` when
#'   it declares any.
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
#' # The entries are the lower triangle of the matrix.
#' mv_sigma(d, theta)
#'
#' # The closed-form method reports standard deviations and a correlation in
#' # place of variances and a covariance.
#' mv_derived(d, theta)$value
#'
#' # The Dirichlet uses this method: its covariance, then its probabilities.
#' dd <- dirichlet_distrib(3)
#' mv_derived(dd, list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8))$block
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
#' entries into standard deviations and correlations and computes the Jacobian
#' explicitly. Writing \eqn{\Sigma = DRD} with \eqn{D} the
#' diagonal of standard deviations,
#' \deqn{\frac{\partial s_j}{\partial\theta_k} = \frac{A_k[j,j]}{2 s_j},
#'   \qquad
#'   \frac{\partial \rho_{jl}}{\partial\theta_k}
#'     = \frac{A_k[j,l]}{s_j s_l}
#'       - \frac{\rho_{jl}}{2}\left(\frac{A_k[j,j]}{\Sigma_{jj}}
#'         + \frac{A_k[l,l]}{\Sigma_{ll}}\right),}
#' with \eqn{A_k = \partial\Sigma/\partial\theta_k}. The function computes no
#' further derivative: the arrays are those that the caller obtained from the
#' parametrization.
#'
#' @param sigma A \eqn{p \times p} symmetric positive definite numeric matrix.
#' @param a A list of \eqn{p \times p} numeric matrices, one per parameter of
#'   the distribution and in `params` order, holding
#'   \eqn{\partial\Sigma/\partial\theta_k}. For a parameter that does not
#'   affect the matrix, the entry is `NULL` or a matrix of zeros; both are
#'   accepted.
#' @param params The distribution's parameter names, which become the column
#'   names of the Jacobian.
#' @param sd_label The prefix for the diagonal quantities. Defaults to `"sd"`;
#'   the Student t passes `"scale_sd"`, because its diagonal quantities are not
#'   standard deviations of the response.
#' @param cor_label The prefix for the off-diagonal quantities. Defaults to
#'   `"cor"`.
#' @param sd_block,cor_block The block labels by which the printed summary
#'   groups the quantities.
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
#' @seealso [mv_sigma_derivs()] for the arrays that it takes,
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
#' # The Jacobian is compared with a numerical one taken through mv_sigma().
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
#' Returns \eqn{\partial\Sigma/\partial\theta_k} for every parameter of a
#' multivariate distribution, in `distrib@params` order, so that the list can
#' be indexed by parameter position without any bookkeeping at the call site.
#' The parametrization supplies the derivatives in its own free values; the
#' location components and, for a Student \eqn{t}, the degrees of freedom leave
#' the matrix alone, and their entries of the list are `NULL`.
#'
#' @details
#' Where the parametrization carries the precision the arrays are carried onto
#' the covariance by \eqn{\partial\Sigma/\partial\eta_k = -\Sigma A_k \Sigma},
#' so the result is always in the covariance, on whichever side the model is
#' written.
#'
#' @param distrib A [multivariate_distrib()] object with a `param` property.
#' @param theta A named list of parameters, already aligned.
#' @param n_before The number of parameters before the free values of the
#'   matrix parametrization, which is \eqn{p} for the gaussian and the Student
#'   \eqn{t}.
#'
#' @return A list as long as `distrib@params`. The entries at the positions of
#'   the free values of the parametrization are \eqn{p \times p} numeric
#'   matrices; every other entry is `NULL`.
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
#' # One entry per parameter, and the two mean entries are NULL.
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


#' @title Interpretable Quantities of a Multivariate Gaussian
#' @name mv_derived.MvGaussianDistrib
#'
#' @description
#' Returns the quantities that the chosen parametrization describes, with their
#' Jacobian. Where the matrix parametrization declares its own block through
#' `parameters7::param_readable()` (a compound symmetry, an AR(1), an AR(p)),
#' only that block is returned, and the common variance `scale` is reported as
#' its square root `sd` on the covariance side. Otherwise a covariance
#' parametrization reports the standard deviations and the correlations, and a
#' precision parametrization reports the conditional variances
#' \eqn{1/\Omega_{jj} = \operatorname{Var}(Y_j \mid Y_{-j})} and the partial
#' correlations \eqn{-\Omega_{jk}/\sqrt{\Omega_{jj}\Omega_{kk}}}, the
#' correlation of two coordinates given all the others.
#'
#' @details
#' # One parametrization, one set of quantities
#'
#' A model written on the precision is reported on the precision, so the
#' marginal standard deviations and correlations are not reported beside its own
#' quantities; they remain available from `mv_sigma()` and `variance()`. In the
#' same way a structured matrix is fixed by its few parameters, and the standard
#' deviation of every coordinate and the correlation of every pair would repeat
#' them.
#'
#' # What a precision's diagonal means
#'
#' The diagonal of a precision is read through the conditional variance, so
#' \eqn{1/\Omega_{jj}} is reported in place of a square root of
#' \eqn{\Omega_{jj}}. Its ratio to the marginal variance is \eqn{1 - R_j^2} for
#' the regression of that coordinate on all the others. At \eqn{p = 2} the
#' partial correlation equals the correlation.
#'
#' @param distrib An [MvGaussianDistrib] object, from [mvgaussian1_distrib()]
#'   or [mvgaussian2_distrib()].
#' @param theta A named list of parameters, already aligned by the generic.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list with `value`, `jacobian`, `transform` and `block`, as
#'   [mv_derived()] documents: the structure's own block where it declares one;
#'   otherwise \eqn{p} standard deviations and \eqn{p(p-1)/2} correlations for
#'   a covariance, or \eqn{p} conditional variances and \eqn{p(p-1)/2} partial
#'   correlations for a precision.
#'
#' @section Notation:
#' \eqn{\Sigma} is the covariance, \eqn{\Omega = \Sigma^{-1}} the precision,
#' \eqn{p} the dimension, \eqn{R_j^2} the coefficient of determination of the
#' regression of coordinate \eqn{j} on the others, and \eqn{Y_{-j}} the
#' response without that coordinate.
#'
#' @seealso [mv_sd_cor()] for the closed-form Jacobian,
#'   [mv_own_block()] for the structure's own block, [mv_summary()] for the
#'   printed result, and [mv_derived()] for the generic.
#'
#' @examples
#' d <- mvgaussian1_distrib(2)
#' theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
#'               sigma_L2.1 = 0.5)
#' mv_derived(d, theta)$value
#'
#' # The precision side reports the conditional variances and the partial
#' # correlation, and nothing of the covariance.
#' o <- mvgaussian2_distrib(2, parameters7::log_cholesky(2))
#' th_o <- list(mu1 = 0, mu2 = 0, omega_log_L1 = 0, omega_log_L2 = 0,
#'              omega_L2.1 = 0.5)
#' od <- mv_derived(o, th_o)
#' od$value
#'
#' # A conditional variance is 1 / Omega_jj, and is at most the marginal one.
#' Om <- parameters7::param_value(o@param, unlist(th_o)[3:5])
#' c(conditional = 1 / Om[1, 1], marginal = mv_sigma(o, th_o)[1, 1])
#'
#' # A compound symmetry reports its standard deviation and its correlation.
#' cs <- mvgaussian1_distrib(3, parameters7::compound_symmetry(3))
#' mv_derived(cs, list(mu1 = 0, mu2 = 0, mu3 = 0, sigma_log_scale = log(4),
#'                     sigma_logit_rho = 0))$value
#'
#' @keywords internal
S7::method(mv_derived, MvGaussianDistrib) <- function(distrib, theta, ...) {
  own <- mv_own_block(distrib, theta, sd_label = "sd")
  if (!is.null(own)) return(own)

  p <- distrib@n_dim
  if (!isTRUE(distrib@inverted)) {
    sigma <- unname(mv_sigma(distrib, theta))
    a <- mv_sigma_derivs(distrib, theta, n_before = p)
    return(mv_sd_cor(sigma, a, distrib@params))
  }

  # The precision's own reading. Omega and its derivatives are the matrix
  # parameter's matrix directly, with no inversion.
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
  # precision's diagonal means is the conditional variance 1/Omega_jj.
  pc <- mv_sd_cor(omega, aw, distrib@params)
  keep <- grepl("^cor_", names(pc$value))
  nm_pc <- sub("^cor_", "pcor_", names(pc$value)[keep])
  cvar <- 1 / diag(omega)
  nm_cvar <- sprintf("cvar_v%d", seq_len(p))
  jac_cvar <- matrix(0, p, distrib@n_params,
                     dimnames = list(nm_cvar, distrib@params))
  for (l in seq_along(aw)) {
    if (is.null(aw[[l]])) next
    # d(1/w_jj)/d eta = -A[j,j] / w_jj^2
    jac_cvar[, l] <- -diag(aw[[l]]) / diag(omega)^2
  }

  list(
    value = c(stats::setNames(cvar, nm_cvar),
              stats::setNames(-pc$value[keep], nm_pc)),
    jacobian = rbind(jac_cvar,
                     `rownames<-`(-pc$jacobian[keep, , drop = FALSE], nm_pc)),
    transform = c(stats::setNames(rep("log", p), nm_cvar),
                  stats::setNames(rep("atanh", sum(keep)), nm_pc)),
    block = c(stats::setNames(rep("Conditional variances", p), nm_cvar),
              stats::setNames(rep("Partial correlations", sum(keep)), nm_pc))
  )
}


#' @title Scale Standard Deviations and Correlations of a Multivariate t
#' @name mv_derived.MvStudentTDistrib
#'
#' @description
#' Returns the square roots of the diagonal of the scale matrix and the
#' correlations that it implies, with the closed-form Jacobian that
#' [mv_sd_cor()] supplies. The correlations are also those of the response: the
#' covariance is \eqn{\nu\Sigma/(\nu-2)} where it exists, and a positive
#' multiple of a matrix has the same correlations. The diagonal quantities are
#' not standard deviations of the response, and for that reason they are named
#' `scale_sd_v1`, ..., `scale_sd_vp` and blocked as `"Scale standard
#' deviations"`.
#'
#' The rule of [mv_derived.MvGaussianDistrib()] holds here as well: a
#' structured parametrization reports only its own block, its `scale` as the
#' square root `scale_sd` on the scale-matrix side, and the inverse
#' parametrization reports the conditional scales \eqn{1/M_{jj}} and the
#' partial correlations and nothing of the scale matrix.
#'
#' @details
#' The quantities do not include the degrees of freedom: \eqn{\nu} is
#' interpretable on its own scale, so [confint.distrib_fit()] reports it and
#' this generic does not repeat it. For the same reason the \eqn{\nu} column of
#' the Jacobian is zero.
#'
#' @param distrib An [MvStudentTDistrib] object, from
#'   [mvstudent_t1_distrib()] or [mvstudent_t2_distrib()].
#' @param theta A named list of parameters, already aligned by the generic.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list with `value`, `jacobian`, `transform` and `block`, as
#'   [mv_derived()] documents: the structure's own block where it declares
#'   one; otherwise \eqn{p} scale standard deviations on the log scale and
#'   \eqn{p(p-1)/2} correlations on Fisher's \eqn{z}, or for the inverse
#'   parametrization \eqn{p} conditional scales and \eqn{p(p-1)/2} partial
#'   correlations.
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
#' # The correlation is the same either way, because the two matrices are
#' # positive multiples of each other.
#' cv <- function(m) m[1, 2] / sqrt(m[1, 1] * m[2, 2])
#' c(reported = der$value[["cor_v1_v2"]],
#'   from_covariance = cv(variance(d, theta)))
#'
#' # The Jacobian column of nu is zero.
#' der$jacobian[, "nu"]
#'
#' @keywords internal
S7::method(mv_derived, MvStudentTDistrib) <- function(distrib, theta, ...) {
  own <- mv_own_block(distrib, theta, sd_label = "scale_sd")
  if (!is.null(own)) return(own)

  p <- distrib@n_dim
  if (!isTRUE(distrib@inverted)) {
    sigma <- unname(mv_sigma(distrib, theta))
    a <- mv_sigma_derivs(distrib, theta, n_before = p)
    return(mv_sd_cor(sigma, a, distrib@params,
                     sd_label = "scale_sd",
                     sd_block = "Scale standard deviations"))
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
  keep <- grepl("^cor_", names(pc$value))
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

  extra
}


#' @title Interpretable Estimates of a Multivariate Fit
#'
#' @description
#' Reports the standard deviations, correlations and the other quantities that
#' [mv_derived()] declares for the fitted distribution, each with its standard
#' error and its confidence interval. A multivariate fit estimates the free
#' values of a \pkg{parameters7} parametrization, and those are coordinates: the
#' estimate and standard error of `sigma_log_L2` do not describe a quantity of
#' interest. This function carries the variance matrix of the fit onto
#' quantities that do.
#'
#' @details
#' # The delta method
#'
#' \deqn{\widehat{\mathrm{Var}}\{g(\hat\theta)\}
#'   = J \, \widehat{\mathrm{Var}}(\hat\theta) \, J^\top,
#'   \qquad J = \partial g/\partial\theta,}
#' where \eqn{J} is taken from [mv_derived()], in closed form for the
#' multivariate gaussian and Student \eqn{t} and by one central difference for
#' the Dirichlet, the multinomial and any other family.
#'
#' # Where each interval is built
#'
#' Each interval is built on the scale that its quantity declares and then
#' mapped back, as [fit_distrib()] does for a univariate parameter. The interval
#' of a standard deviation is built on the log scale and that of a correlation
#' on Fisher's \eqn{z}, so each interval stays inside the set in which its
#' quantity lies. An interval on the raw scale can put a correlation above one.
#'
#' The standard error is carried onto that scale by the transform's own
#' derivative, so the reported `Std. Error` column stays on the quantity's
#' natural scale while the limits are computed on the transformed one.
#'
#' @param object A [distrib_fit()] of a multivariate distribution. A univariate
#'   fit is rejected: there the parameters are already the interpretable
#'   quantities and [confint.distrib_fit()] is the function to call.
#' @param level The confidence level, a single number strictly inside
#'   \eqn{(0, 1)}. Defaults to the fit's own `level`. Anything else is an error.
#' @param information An [observed_info()] object (the default) or an
#'   [expected_info()] object, which selects the information matrix whose
#'   inverse gives the standard errors; see [vcov.distrib_fit()].
#'
#' @return A data frame with one row per quantity and four columns:
#'   `Estimate`, `Std. Error`, and the two confidence limits, named after the
#'   percentiles that they represent (`2.5%` and `97.5%` at the default level).
#'   Its row names are the quantity names, and it carries the attribute
#'   `"block"`, a character vector naming the group to which each row belongs
#'   and by which `print()` lays the summary out.
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
#' # The summary reports the standard deviations and the correlation.
#' mv_summary(fit)
#'
#' # The coordinates on which the fit was estimated are shown beside them.
#' confint(fit)
#'
#' # Every interval stays inside the set of its quantity: the standard
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
mv_summary <- function(object, level = object@level,
                       information = observed_info()) {
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
  v <- fit_report(object, information)$vcov
  # a negative variance gives NaN, not a standard error of zero
  var_d <- diag(der$jacobian %*% v %*% t(der$jacobian))
  var_d[!is.na(var_d) & var_d < 0] <- NaN
  se <- sqrt(var_d)

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
