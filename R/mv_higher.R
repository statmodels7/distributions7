#' @include partition_sums.R
#' @include dirichlet_distrib.R
#' @include multinomial_distrib.R
#' @include mvstudent_t_distrib.R
NULL

# Third and fourth derivatives of the two simplex-valued families.
#
# Both log-densities are a sum of terms each depending on ONE coordinate of a
# simplex -- log Gamma(alpha_j) for the Dirichlet, log p_j for the multinomial
# -- so the multivariate Faa di Bruno collapses to one univariate partition
# sum per coordinate, and nothing needs the mixed derivative tensors a general
# composition would. The inner derivatives are written out; the map's come from
# parameters7, whose contract carries orders three and four.

#' Derivative Tensors of a Simplex-Valued Map
#'
#' @description
#' Collects \eqn{\partial^{S}\mu_j} for every multiset \eqn{S} of free indices
#' up to the requested order, keyed by the sorted tuple, as a list of numeric
#' vectors over the coordinates \eqn{j}.
#'
#' @param s A \pkg{parameters7} `simplex` parameter.
#' @param eta The free vector.
#' @param order The highest order required, 1 to 4.
#'
#' @return A named list of numeric vectors, keyed as `"1"`, `"1,2"`
#'   and so on.
#'
#' @keywords internal
simplex_map_tensors <- function(s, eta, order) {
  out <- list()
  for (k in seq_len(order)) {
    d <- switch(k,
      parameters7::param_d1(s, eta), parameters7::param_d2(s, eta),
      parameters7::param_d3(s, eta), parameters7::param_d4(s, eta)
    )
    tup <- parameters7::param_tuple_indices(s, k)
    for (i in seq_along(tup)) {
      out[[paste(sort(tup[[i]]), collapse = ",")]] <- as.numeric(d[[i]])
    }
  }
  out
}

#' One Univariate Chain Rule Over a Multivariate Map
#'
#' @description
#' Evaluates \eqn{\partial^{S} f(u(v))} for a scalar function \eqn{f} of one
#' coordinate, by the partition form of Faa di Bruno
#' \deqn{\partial^{S} f(u) = \sum_{\pi} f^{(\lvert\pi\rvert)}(u)
#'       \prod_{B \in \pi} \partial^{B} u,}
#' the sum running over the set partitions of the multiset \eqn{S}. Everything
#' is vectorized over the coordinates, so one call serves every \eqn{j} at
#' once.
#'
#' @param tuple An integer vector of free-value indices, with repeats.
#' @param fd A list whose \eqn{m}-th element is \eqn{f^{(m)}} evaluated at
#'   every coordinate.
#' @param ud A named list of the map's derivative tensors, keyed by sorted
#'   tuple as [simplex_map_tensors()] returns them.
#'
#' @return A numeric vector over the coordinates.
#'
#' @keywords internal
chain_univariate <- function(tuple, fd, ud) {
  acc <- 0
  for (part in index_partitions(tuple)) {
    term <- fd[[length(part)]]
    for (b in part) term <- term * ud[[paste(sort(b), collapse = ",")]]
    acc <- acc + term
  }
  acc
}

#' Derivative Tensors of a Dirichlet's Shape Vector
#'
#' @description
#' \eqn{\alpha = \phi\,\mu(\eta)} is bilinear in the concentration and the
#' mean, so \eqn{\partial^{S}\alpha} is \eqn{\phi\,\partial^{S}\mu} when
#' \eqn{S} names no \eqn{\phi}, is \eqn{\partial^{S'}\mu} when it names one,
#' and is zero when it names two or more.
#'
#' @param s A \pkg{parameters7} `simplex` parameter.
#' @param eta The mean's free vector.
#' @param phi The concentration.
#' @param order The highest order required, 1 to 4.
#'
#' @return A named list of numeric vectors over the coordinates, keyed by
#'   sorted tuple over the composite index set, \eqn{\phi} last.
#'
#' @keywords internal
dirichlet_map_tensors <- function(s, eta, phi, order) {
  k <- s@n_free
  mu <- simplex_map_tensors(s, eta, min(order, 4L))
  mu0 <- as.numeric(parameters7::param_value(s, eta))
  zero <- rep(0, length(mu0))
  out <- list()
  for (o in seq_len(order)) {
    for (t in tuple_indices_upto(k + 1L, o)) {
      nphi <- sum(t == k + 1L)
      rest <- t[t <= k]
      key <- paste(sort(t), collapse = ",")
      out[[key]] <- if (nphi >= 2L) {
        zero
      } else if (nphi == 1L) {
        if (!length(rest)) mu0 else mu[[paste(sort(rest), collapse = ",")]]
      } else {
        phi * mu[[paste(sort(rest), collapse = ",")]]
      }
    }
  }
  out
}

#' Index Tuples of a Given Width Over a Number of Variables
#'
#' @description
#' The same enumeration `parameters7::param_tuple_indices()` uses, taken
#' over a count rather than over a parameter, so that a composite index set
#' formed by appending a coordinate can be enumerated without building an
#' object for it.
#'
#' @param d The number of variables.
#' @param order The tuple width, 1 to 4.
#'
#' @return A list of integer vectors.
#'
#' @keywords internal
tuple_indices_upto <- function(d, order) {
  numericals7::tuple_indices(d, order)
}

#' Expected Higher Derivatives of a Multivariate Family by Sampling
#'
#' @description
#' Averages the closed-form observed derivatives over a sample from the family,
#' which is the route the multivariate branch takes throughout: a quadrature
#' over the simplex has no counterpart to the one-dimensional split at
#' quantiles, and the Bartlett route would need the score's own higher
#' derivatives.
#'
#' @param distrib A multivariate distribution object.
#' @param y The observed response, read only for its number of rows.
#' @param theta A named list of parameters.
#' @param order The derivative order, 3 or 4.
#' @param nsim The Monte Carlo sample size.
#'
#' @return A named list of numeric vectors, each constant across observations.
#'
#' @keywords internal
mv_expected_higher <- function(distrib, y, theta, order, nsim) {
  big <- distrib_rng(distrib, nsim, theta)
  ob <- if (S7::S7_inherits(distrib, DirichletDistrib)) {
    dirichlet_higher(distrib, big, theta, order)
  } else {
    multinomial_higher(distrib, big, theta, order)
  }
  n <- n_obs(distrib, y)
  lapply(ob, function(v) rep(mean(v), n))
}


#' @title Dirichlet Analytical Third and Fourth Derivatives
#' @name distrib_deriv3.DirichletDistrib
#' @description
#' Closed form. The log-density is
#' \eqn{\log\Gamma(\phi) - \sum_j \log\Gamma(\alpha_j) +
#' \sum_j(\alpha_j - 1)\log y_j} with \eqn{\alpha = \phi\mu(\eta)}, so every
#' term depends on one coordinate of the simplex and the chain rule is one
#' univariate partition sum per coordinate. The concentration also enters
#' directly, through \eqn{\log\Gamma(\phi)}, which contributes only to the
#' component all of whose indices name it.
#' @param distrib A `DirichletDistrib` object.
#' @param y A matrix with one row per observation.
#' @param theta A named list of parameters.
#' @param expected Whether to return expected derivatives.
#' @param scale Either `"parameter"` or `"link"`; handled by the generic.
#' @param approx The approximation used when `expected` is `TRUE`.
#' @param nsim The Monte Carlo sample size.
#' @param ... Unused.
#' @return A named list of third-derivative components.
#' @seealso [dirichlet_distrib()]
S7::method(distrib_deriv3, DirichletDistrib) <- function(distrib, y, theta,
                                                         expected = FALSE,
                                                         scale = c("parameter", "link"),
                                                         approx = c("integrate", "bartlett", "mc", "opg"),
                                                         nsim = 10000, ...) {
  if (expected) return(mv_expected_higher(distrib, y, theta, 3L, nsim))
  dirichlet_higher(distrib, y, theta, 3L)
}

#' @rdname distrib_deriv3.DirichletDistrib
#' @name distrib_deriv4.DirichletDistrib
S7::method(distrib_deriv4, DirichletDistrib) <- function(distrib, y, theta,
                                                         expected = FALSE,
                                                         scale = c("parameter", "link"),
                                                         approx = c("integrate", "bartlett", "mc", "opg"),
                                                         nsim = 10000, ...) {
  if (expected) return(mv_expected_higher(distrib, y, theta, 4L, nsim))
  dirichlet_higher(distrib, y, theta, 4L)
}

#' Assemble a Dirichlet's Higher Derivatives
#'
#' @description
#' Runs the univariate chain rule over the shape vector for every component of
#' the requested order and adds the two terms that do not pass through it: the
#' response, which enters \eqn{\alpha_j} linearly, and \eqn{\log\Gamma(\phi)}.
#'
#' @param distrib A `DirichletDistrib` object.
#' @param y A matrix with one row per observation.
#' @param theta A named list of parameters.
#' @param order The derivative order, 3 or 4.
#'
#' @return A named list of numeric vectors, one per component.
#'
#' @keywords internal
dirichlet_higher <- function(distrib, y, theta, order) {
  s <- distrib@param
  v <- mv_flat_theta(distrib, theta)
  k <- s@n_free
  eta <- v[seq_len(k)]
  phi <- v[[k + 1L]]
  alpha <- phi * as.numeric(parameters7::param_value(s, eta))

  y <- if (is.matrix(y)) y else matrix(y, nrow = 1L)
  ly <- log(y)

  ud <- dirichlet_map_tensors(s, eta, phi, order)
  # f(a) = -log Gamma(a), so f^(m)(a) = -psi^(m-1)(a)
  fd <- lapply(seq_len(order), function(m) -psigamma(alpha, m - 1L))

  idx <- deriv_indices(distrib@params, order)
  out <- lapply(idx, function(t) {
    key <- paste(sort(t), collapse = ",")
    # the log-Gamma of the concentration enters only the component all of
    # whose indices name it, and neither term below carries the response
    direct <- if (all(t == k + 1L)) psigamma(phi, order - 1L) else 0
    base <- sum(chain_univariate(t, fd, ud)) + direct
    # the response enters alpha linearly, so it multiplies the map's own
    # tensor and is the only term that varies by observation
    base + as.numeric(ly %*% ud[[key]])
  })
  stats::setNames(out, deriv_names(distrib@params, order))
}


#' @title Multinomial Analytical Third and Fourth Derivatives
#' @name distrib_deriv3.MultinomialDistrib
#' @description
#' Closed form. The log-mass is \eqn{\sum_j y_j \log p_j} up to a constant, so
#' each term depends on one coordinate of the simplex and the chain rule is one
#' univariate partition sum per coordinate, with
#' \eqn{f^{(m)}(p) = (-1)^{m-1}(m-1)!\,p^{-m}}.
#' @param distrib A `MultinomialDistrib` object.
#' @param y A matrix with one row per observation.
#' @param theta A named list of parameters.
#' @param expected Whether to return expected derivatives.
#' @param scale Either `"parameter"` or `"link"`; handled by the generic.
#' @param approx The approximation used when `expected` is `TRUE`.
#' @param nsim The Monte Carlo sample size.
#' @param ... Unused.
#' @return A named list of third-derivative components.
#' @seealso [multinomial_distrib()]
S7::method(distrib_deriv3, MultinomialDistrib) <- function(distrib, y, theta,
                                                            expected = FALSE,
                                                            scale = c("parameter", "link"),
                                                            approx = c("integrate", "bartlett", "mc", "opg"),
                                                            nsim = 10000, ...) {
  if (expected) return(mv_expected_higher(distrib, y, theta, 3L, nsim))
  multinomial_higher(distrib, y, theta, 3L)
}

#' @rdname distrib_deriv3.MultinomialDistrib
#' @name distrib_deriv4.MultinomialDistrib
S7::method(distrib_deriv4, MultinomialDistrib) <- function(distrib, y, theta,
                                                            expected = FALSE,
                                                            scale = c("parameter", "link"),
                                                            approx = c("integrate", "bartlett", "mc", "opg"),
                                                            nsim = 10000, ...) {
  if (expected) return(mv_expected_higher(distrib, y, theta, 4L, nsim))
  multinomial_higher(distrib, y, theta, 4L)
}

#' Assemble a Multinomial's Higher Derivatives
#'
#' @description
#' Runs the univariate chain rule over the probability vector for every
#' component of the requested order, weighted by the summed counts.
#'
#' @param distrib A `MultinomialDistrib` object.
#' @param y A matrix with one row per observation.
#' @param theta A named list of parameters.
#' @param order The derivative order, 3 or 4.
#'
#' @return A named list of numeric vectors, one per component.
#'
#' @keywords internal
multinomial_higher <- function(distrib, y, theta, order) {
  s <- distrib@param
  eta <- mv_flat_theta(distrib, theta)
  prob <- as.numeric(parameters7::param_value(s, eta))

  y <- if (is.matrix(y)) y else matrix(y, nrow = 1L)

  ud <- simplex_map_tensors(s, eta, order)
  # f(p) = log p, so f^(m)(p) = (-1)^(m-1) (m-1)! p^-m
  fd <- lapply(seq_len(order), function(m) {
    (-1)^(m - 1L) * factorial(m - 1L) / prob^m
  })

  idx <- deriv_indices(distrib@params, order)
  out <- lapply(idx, function(t) as.numeric(y %*% chain_univariate(t, fd, ud)))
  stats::setNames(out, deriv_names(distrib@params, order))
}


# --- the multivariate Student t ---------------------------------------------
#
# The log-density is
#   c(nu) - (1/2) log|Sigma| - ((nu+p)/2) log(1 + q/nu),   q = r' Sigma^-1 r,
# and every term is elementary in nu: unlike the univariate skew t, whose nu
# derivatives carry the derivative of a Student t distribution function in its
# degrees of freedom, nothing here is a distribution function. So all four
# orders close.

#' Derivatives of the Quadratic Form of a Multivariate Student t
#'
#' @description
#' Evaluates \eqn{\partial^{B} q} for a multiset \eqn{B} of mean and matrix
#' indices, with \eqn{q = r'\Sigma^{-1}r}. The form is quadratic in the mean,
#' so a block naming three or more mean coordinates is zero.
#'
#' @param b An integer vector of composite indices, mean coordinates first.
#' @param r The \eqn{n \times p} matrix of residuals.
#' @param pget The accessor for \eqn{\partial^{t}\Sigma^{-1}}, as returned by
#'   [mvg_ptensors()].
#' @param p The dimension.
#'
#' @return A numeric vector of length `nrow(r)`, or a scalar recycled by
#'   the caller.
#'
#' @keywords internal
mvt_q_deriv <- function(b, r, pget, p) {
  is_mu <- b <= p
  n_mu <- sum(is_mu)
  n <- nrow(r)
  if (n_mu >= 3L) return(rep(0, n))
  et <- b[!is_mu] - p
  P <- pget(et)
  if (n_mu == 0L) return(rowSums((r %*% P) * r))
  if (n_mu == 1L) return(-2 * (r %*% P)[, b[is_mu]])
  ij <- b[is_mu]
  rep(2 * P[ij[1L], ij[2L]], n)
}

#' Assemble a Multivariate Student t's Higher Derivatives
#'
#' @description
#' Splits a component into its mean-and-matrix part and its \eqn{\nu} part.
#' The log-determinant and the \eqn{\nu} constant each contribute to one kind
#' of component only; the remaining term is a univariate chain rule in \eqn{q}
#' whose outer derivatives, \eqn{(-1)^{m-1}(m-1)!/(\nu+q)^{m}}, are then
#' differentiated in \eqn{\nu} by the Leibniz rule against the prefactor
#' \eqn{(\nu+p)/2}, which is linear.
#'
#' @param distrib A `MvStudentTDistrib` object.
#' @param y An \eqn{n \times p} matrix of observations.
#' @param theta A named list of parameters.
#' @param order The derivative order, 1 to 4.
#'
#' @return A named list of numeric vectors, one per component.
#'
#' @keywords internal
mvt_higher <- function(distrib, y, theta, order) {
  y <- as_mv_matrix(distrib, y)
  n <- nrow(y)
  p <- distrib@n_dim
  pc <- mvt_pieces(distrib, theta)
  nu <- pc$nu
  s <- pc$s
  nfree <- s@n_free
  inu <- p + nfree + 1L
  r <- sweep(y, 2L, pc$mu, "-")
  q <- rowSums((r %*% pc$sigma_inv) * r)

  pt <- mvg_ptensors(pc, order)
  ldfun <- switch(order,
    parameters7::param_dlogdet, parameters7::param_d2logdet,
    parameters7::param_d3logdet, parameters7::param_d4logdet
  )
  ld <- ldfun(s, pc$eta)
  ldkey <- vapply(parameters7::param_tuple_indices(s, order),
                  function(t) paste(sort(t), collapse = ","), character(1))

  # d^a/dnu^a of c(nu) = lgamma((nu+p)/2) - lgamma(nu/2) - (p/2) log(nu pi)
  cnu <- function(a) {
    0.5^a * (psigamma((nu + p) / 2, a - 1L) - psigamma(nu / 2, a - 1L)) -
      (p / 2) * (-1)^(a - 1L) * factorial(a - 1L) / nu^a
  }
  # d^a/dnu^a of (nu + q)^-m
  wa <- function(m, a) {
    if (a < 0L) return(rep(0, n))
    (-1)^a * (factorial(m + a - 1L) / factorial(m - 1L)) / (nu + q)^(m + a)
  }

  idx <- deriv_indices(distrib@params, order)
  out <- lapply(idx, function(t) {
    a <- sum(t == inu)
    rest <- t[t != inu]
    acc <- rep(0, n)

    # the constant in nu, present only when every index names nu
    if (a == order) acc <- acc + cnu(a)

    # the log-determinant, present only when every index names the matrix
    if (a == 0L && length(rest) == order && all(rest > p)) {
      acc <- acc - 0.5 * ld[[match(paste(sort(rest - p), collapse = ","), ldkey)]]
    }

    # -((nu+p)/2) log(1 + q/nu), differentiated in the rest and then in nu
    if (!length(rest)) {
      # v(nu) = log(nu + q) - log(nu)
      v <- if (a == 0L) base::log1p(q / nu) else
        (-1)^(a - 1L) * factorial(a - 1L) * ((nu + q)^-a - nu^-a)
      vprev <- if (a == 0L) rep(0, n) else if (a == 1L) base::log1p(q / nu) else
        (-1)^(a - 2L) * factorial(a - 2L) * ((nu + q)^-(a - 1L) - nu^-(a - 1L))
      acc <- acc - ((nu + p) / 2 * v + a / 2 * vprev)
    } else {
      for (part in index_partitions(rest)) {
        m <- length(part)
        term <- rep(1, n)
        for (b in part) term <- term * mvt_q_deriv(b, r, pt$get, p)
        coef <- (-1)^(m - 1L) * factorial(m - 1L) *
          ((nu + p) / 2 * wa(m, a) + a / 2 * wa(m, a - 1L))
        acc <- acc - coef * term
      }
    }
    acc
  })
  stats::setNames(out, deriv_names(distrib@params, order))
}


#' @title Multivariate Student t Third and Fourth Derivatives
#' @name distrib_deriv3.MvStudentTDistrib
#' @description
#' Closed form. Every term of the log-density is elementary in \eqn{\nu}, so
#' the obstruction the univariate skew t meets -- the derivative of a Student t
#' distribution function in its degrees of freedom -- does not arise here, and
#' all four orders close.
#' @param distrib A `MvStudentTDistrib` object.
#' @param y An \eqn{n \times p} matrix of observations.
#' @param theta A named list of parameters.
#' @param expected Logical; the expectation is approximated by sampling.
#' @param scale Either `"parameter"` or `"link"`; handled by the generic.
#' @param approx Strategy label; sampling is the only multivariate route.
#' @param nsim Monte Carlo sample size.
#' @param ... Unused.
#' @return A named list of third-derivative component vectors.
#' @seealso [mvstudent_t_distrib()]
S7::method(distrib_deriv3, MvStudentTDistrib) <- function(distrib, y, theta,
                                                           expected = FALSE,
                                                           scale = c("parameter", "link"),
                                                           approx = c("integrate", "bartlett", "mc", "opg"),
                                                           nsim = 10000, ...) {
  if (expected) {
    big <- distrib_rng(distrib, nsim, theta)
    ob <- mvt_higher(distrib, big, theta, 3L)
    n <- n_obs(distrib, y)
    return(lapply(ob, function(v) rep(mean(v), n)))
  }
  mvt_higher(distrib, y, theta, 3L)
}

#' @rdname distrib_deriv3.MvStudentTDistrib
#' @name distrib_deriv4.MvStudentTDistrib
S7::method(distrib_deriv4, MvStudentTDistrib) <- function(distrib, y, theta,
                                                           expected = FALSE,
                                                           scale = c("parameter", "link"),
                                                           approx = c("integrate", "bartlett", "mc", "opg"),
                                                           nsim = 10000, ...) {
  if (expected) {
    big <- distrib_rng(distrib, nsim, theta)
    ob <- mvt_higher(distrib, big, theta, 4L)
    n <- n_obs(distrib, y)
    return(lapply(ob, function(v) rep(mean(v), n)))
  }
  mvt_higher(distrib, y, theta, 4L)
}
