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

#' @title Derivative Arrays of a Simplex-Valued Map
#'
#' @description
#' Collects \eqn{\partial^{S}\mu_j} for every multiset \eqn{S} of free indices
#' up to the requested order, keyed by the sorted tuple, as a list of numeric
#' vectors over the coordinates \eqn{j}. The arrays come from
#' \pkg{parameters7}, which supplies orders three and four for a `simplex`. The
#' function derives nothing itself and only adds the keying, which lets the
#' chain rule in [chain_univariate()] look a block up by its index multiset.
#'
#' @param s A `parameters7::simplex()` parametrization.
#' @param eta The free vector, a numeric vector of length `s@n_free`.
#' @param order The highest order required, a single whole number from 1 to 4.
#'   A larger value signals an error.
#'
#' @return A named list of numeric vectors, each as long as the simplex has
#'   coordinates, keyed `"1"`, `"2"`, `"1,1"`, `"1,2"` and so on, with the
#'   indices sorted and comma-separated. Its length is the number of index
#'   multisets of width 1 to `order`.
#'
#' @section Notation:
#' \eqn{\mu} is the point of the simplex that the parametrization produces,
#' \eqn{\eta} its free vector, and \eqn{S} a multiset of free-value indices.
#'
#' @seealso [chain_univariate()], which consumes these,
#'   [dirichlet_map_tensors()] for the version that carries a concentration as
#'   well, and [parameters7::simplex()] for the parametrization.
#'
#' @examples
#' s <- parameters7::simplex(3)
#' eta <- c(0.3, -0.2)
#' mt <- distributions7:::simplex_map_tensors(s, eta, 2)
#' names(mt)
#'
#' # The first-order entry is the parametrization's own first derivative.
#' mt[["1"]]
#' as.numeric(parameters7::param_d1(s, eta)[[1]])
#'
#' # Every derivative of a point on the simplex sums to zero, the value
#' # summing to one at every eta.
#' vapply(mt, sum, numeric(1))
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

#' @title One Univariate Chain Rule Over a Multivariate Map
#'
#' @description
#' Evaluates \eqn{\partial^{S} f(u(v))} for a scalar function \eqn{f} of one
#' coordinate, by the partition form of Faa di Bruno
#' \deqn{\partial^{S} f(u) = \sum_{\pi} f^{(\lvert\pi\rvert)}(u)
#'       \prod_{B \in \pi} \partial^{B} u,}
#' the sum running over the set partitions of the multiset \eqn{S}. Everything
#' is vectorized over the coordinates, so one call serves every \eqn{j} at once.
#'
#' @details
#' The two simplex-valued families use this at orders three and four. Their
#' log-densities are sums of terms that each depend on one coordinate,
#' \eqn{\log\Gamma(\alpha_j)} for the Dirichlet and \eqn{\log p_j} for the
#' multinomial, so the multivariate Faa di Bruno formula reduces to one
#' univariate partition sum per coordinate, and the mixed derivative arrays are
#' never formed. The multinomial uses it only on a chart other than the additive
#' log-ratio one, where its derivatives are cumulants (see [mn_cumulant()]).
#'
#' @param tuple An integer vector of free-value indices, with repeats allowed.
#'   Its length is the derivative order.
#' @param fd A list whose \eqn{m}-th element is \eqn{f^{(m)}} evaluated at every
#'   coordinate, a numeric vector each. It must run to `length(tuple)`.
#' @param ud A named list of the map's derivative arrays, keyed by sorted tuple
#'   as [simplex_map_tensors()] returns them. Every block of every partition of
#'   `tuple` must have a key in it.
#'
#' @return A numeric vector over the coordinates, as long as the elements of
#'   `fd`.
#'
#' @section Notation:
#' \eqn{S} is a multiset of free-value indices, \eqn{\pi} a set partition of it,
#' \eqn{B} a block of that partition, and \eqn{u} the map being composed with.
#'
#' @seealso [index_partitions()] for the enumeration it sums over,
#'   [simplex_map_tensors()] for the argument, and [dirichlet_higher()] and
#'   [multinomial_higher()] for the consumers.
#'
#' @examples
#' s <- parameters7::simplex(3)
#' eta <- c(0.3, -0.2)
#' ud <- distributions7:::simplex_map_tensors(s, eta, 2)
#'
#' # Compose f(u) = log(u) with the map: f'(u) = 1/u, f''(u) = -1/u^2.
#' u <- as.numeric(parameters7::param_value(s, eta))
#' fd <- list(1 / u, -1 / u^2)
#'
#' # The mixed second derivative of log(mu_j) in the two free values.
#' got <- distributions7:::chain_univariate(c(1L, 2L), fd, ud)
#' got
#'
#' # Against a difference of the first-order chain rule, taken directly.
#' h <- 1e-5
#' d1 <- function(e) {
#'   uu <- as.numeric(parameters7::param_value(s, e))
#'   as.numeric(parameters7::param_d1(s, e)[[1]]) / uu
#' }
#' max(abs(got - (d1(eta + c(0, h)) - d1(eta - c(0, h))) / (2 * h)))
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

#' @title Derivative Arrays of a Dirichlet's Shape Vector
#'
#' @description
#' Collects \eqn{\partial^{S}\alpha} for every multiset \eqn{S} over the
#' Dirichlet's composite index set, the simplex's free values followed by the
#' concentration. The shape vector \eqn{\alpha = \phi\,\mu(\eta)} is bilinear
#' in the concentration and the mean, so the whole array follows from
#' [simplex_map_tensors()] by three cases: \eqn{\phi\,\partial^{S}\mu} when
#' \eqn{S} names no \eqn{\phi}, \eqn{\partial^{S'}\mu} when it names one, and
#' exactly zero when it names two or more.
#'
#' @param s A `parameters7::simplex()` parametrization of the mean.
#' @param eta The mean's free vector, a numeric vector of length `s@n_free`.
#' @param phi The concentration, a single positive number.
#' @param order The highest order required, a single whole number from 1 to 4.
#'
#' @return A named list of numeric vectors over the coordinates, keyed by sorted
#'   tuple over the composite index set with \eqn{\phi} last, so at
#'   `s@n_free = 2` the concentration is index 3 and the key `"3,3"` holds a
#'   vector of zeros.
#'
#' @section Notation:
#' \eqn{\alpha} is the Dirichlet's shape vector, \eqn{\phi} its concentration,
#' \eqn{\mu} the mean on the simplex, \eqn{\eta} the mean's free vector, and
#' \eqn{S'} the multiset \eqn{S} with its one \eqn{\phi} index removed.
#'
#' @seealso [simplex_map_tensors()] for the arrays it lifts,
#'   [tuple_indices_upto()] for the enumeration over the composite index set,
#'   and [dirichlet_higher()] for the consumer.
#'
#' @examples
#' s <- parameters7::simplex(3)
#' eta <- c(0.3, -0.2)
#' dt <- distributions7:::dirichlet_map_tensors(s, eta, 8, 2)
#' names(dt)
#'
#' # The key naming phi alone is the mean itself: alpha is phi times mu.
#' dt[["3"]]
#' as.numeric(parameters7::param_value(s, eta))
#'
#' # Naming phi twice gives exactly zero, alpha being linear in it.
#' dt[["3,3"]]
#'
#' # And a key naming no phi is phi times the mean's own derivative.
#' all.equal(dt[["1"]], 8 * as.numeric(parameters7::param_d1(s, eta)[[1]]))
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

#' @title Index Tuples of a Given Width Over a Number of Variables
#'
#' @description
#' Enumerates the index multisets of exactly `order` indices drawn from `1:d`,
#' in the order that `parameters7::param_tuple_indices()` uses. It takes a count
#' of variables in place of a parametrization, so a composite index set formed
#' by appending a coordinate can be enumerated without building an object of
#' that shape. The Dirichlet needs this for its simplex free values followed by
#' its concentration. The work is done by `numericals7::tuple_indices()`, which
#' this function wraps.
#'
#' @param d The number of variables, a single positive whole number.
#' @param order The tuple width, a single whole number from 1 to 4.
#'
#' @return A list of integer vectors, each of length `order`, sorted within
#'   itself. Its length is `choose(d + order - 1, order)`.
#'
#' @seealso [numericals7::tuple_indices()] for the enumeration itself,
#'   [dirichlet_map_tensors()] for the consumer, and [deriv_names()] for the
#'   same enumeration rendered as component names.
#'
#' @examples
#' # Two indices over three variables: the three repeats then the three pairs.
#' distributions7:::tuple_indices_upto(3, 2)
#'
#' # The count is choose(d + order - 1, order).
#' c(got = length(distributions7:::tuple_indices_upto(4, 3)),
#'   expected = choose(4 + 3 - 1, 3))
#'
#' # Every tuple is sorted, which is what makes a comma-joined key canonical.
#' all(vapply(distributions7:::tuple_indices_upto(3, 3),
#'            function(t) !is.unsorted(t), TRUE))
#'
#' @keywords internal
tuple_indices_upto <- function(d, order) {
  numericals7::tuple_indices(d, order)
}

# Commented out on 2026-10-12 (T4): the expected third and fourth orders of the Dirichlet and the multinomial
# are closed form (dirichlet_higher(), multinomial_higher() with
# expected = TRUE), so the sampling route has no caller.
# Kept rather than deleted.
# @title Expected Higher Derivatives of a Simplex-Valued Family by Sampling
#
# @description
# Averages the closed-form observed third or fourth derivatives over a sample
# drawn from the family. This is the route the multivariate branch takes
# throughout: a quadrature over the simplex has no counterpart to the
# one-dimensional split at quantiles, and the Bartlett route would need the
# score's own higher derivatives. The result therefore carries Monte Carlo
# error of order `nsim^(-1/2)`, and a caller who needs it reproducible sets a
# seed first.
#
# @details
# Only the two simplex-valued families reach here. The branch is taken on
# whether `distrib` inherits `DirichletDistrib`, everything else being routed
# to [multinomial_higher()], so this function must not be called for any third
# family.
#
# @param distrib A `DirichletDistrib` or `MultinomialDistrib` object.
# @param y The observed response, read only for its number of rows.
# @param theta A named list of parameters, each component a single number.
# @param order The derivative order, `3L` or `4L`.
# @param nsim The Monte Carlo sample size, a single positive whole number.
#
# @return A named list of numeric vectors of length `n_obs(distrib, y)`, each
#   constant across observations, keyed as `deriv_names(distrib@params, order)`.
#
# @seealso [dirichlet_higher()] and [multinomial_higher()] for the observed
#   derivatives it averages, and [distrib_deriv3.DirichletDistrib()] for the
#   method that calls it.
#
# @examples
# d <- dirichlet_distrib(3)
# theta <- list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8)
# set.seed(1)
# y <- distrib_rng(d, 5, theta)
#
# set.seed(5)
# e3 <- distributions7:::mv_expected_higher(d, y, theta, 3L, 3000)
#
# # One value per observation, and every vector constant.
# lengths(e3)[1:3]
# all(vapply(e3, function(z) diff(range(z)) == 0, TRUE))
#
# # Two seeds give two answers, this being a sample rather than an integral.
# set.seed(6)
# c(first = e3[[1]][1],
#   second = distributions7:::mv_expected_higher(d, y, theta, 3L, 3000)[[1]][1])
#
# @keywords internal
# mv_expected_higher <- function(distrib, y, theta, order, nsim) {
#   big <- distrib_rng(distrib, nsim, theta)
#   ob <- if (S7::S7_inherits(distrib, DirichletDistrib)) {
#     dirichlet_higher(distrib, big, theta, order)
#   } else {
#     multinomial_higher(distrib, big, theta, order)
#   }
#   n <- n_obs(distrib, y)
#   lapply(ob, function(v) rep(mean(v), n))
# }


#' @title Dirichlet Third Derivatives
#' @name distrib_deriv3.DirichletDistrib
#'
#' @description
#' Computes every third derivative of the log-density in the parameters, in
#' closed form. The log-density is
#' \deqn{\ell = \log\Gamma(\phi) - \sum_j \log\Gamma(\alpha_j)
#'   + \sum_j(\alpha_j - 1)\log y_j, \qquad \alpha = \phi\,\mu(\eta),}
#' so every term depends on one coordinate of the simplex and the chain rule
#' reduces to a univariate partition sum per coordinate, run by
#' [chain_univariate()] over the arrays that [dirichlet_map_tensors()] supplies.
#' The concentration also enters directly, through \eqn{\log\Gamma(\phi)}, which
#' contributes only to the component whose indices all name it.
#'
#' @details
#' The shape vector is bilinear in \eqn{\phi} and \eqn{\mu}, so the arrays of
#' the map are zero for a component that names \eqn{\phi} twice or more.
#'
#' On the additive log-ratio chart of [parameters7::simplex()] the term
#' \eqn{-\log\Gamma(\alpha_j)} is split as
#' \eqn{\log\alpha_j - \log\Gamma(1 + \alpha_j)}. The first part equals
#' \eqn{\log\phi + \eta_j - \mathrm{LSE}(\eta)}, with \eqn{\eta_p = 0}, and its
#' derivatives are added in closed form: \eqn{p(-1)^{m-1}(m-1)!/\phi^m} on the
#' component of order \eqn{m} in \eqn{\phi} alone, and \eqn{-p} times a cumulant
#' of the categorical law (see [mn_cumulant()]) on a component that names only
#' the free values of the mean. The second part is regular at a small shape, so
#' the partition sum runs over \eqn{\psi^{(m-1)}(1 + \alpha_j)}. Applied to
#' \eqn{\log\Gamma(\alpha_j)} directly, the chain rule would cancel at a small
#' shape, terms of order one leaving a result of the order of \eqn{\alpha_j}. On
#' another chart a coordinate with \eqn{\alpha_j < 1} is scaled by
#' \eqn{\alpha_j} through
#' \eqn{\alpha^m\psi^{(m-1)}(\alpha) = (-1)^m(m-1)! + \alpha^m\psi^{(m-1)}(1 + \alpha)},
#' so that \eqn{\psi^{(m-1)}(\alpha_j)} does not overflow.
#'
#' With `expected = TRUE` the expectation is exact. The response enters only
#' through \eqn{\sum_j(\alpha_j - 1)\log y_j}, linearly, and
#' \eqn{\mathrm{E}\log y_j = \psi(\alpha_j) - \psi(\phi)}.
#'
#' @param distrib A `DirichletDistrib` object, from [dirichlet_distrib()].
#' @param y A numeric matrix with one row per observation and one column per
#'   coordinate, each row a point of the open simplex. A row with a zero
#'   coordinate is outside the support and gives a non-finite component through
#'   \eqn{\log y_j}.
#' @param theta A named list of parameters, each component a single number:
#'   `mean_alr1`, ..., `mean_alr(p-1)` and `phi`.
#' @param expected Logical of length 1. When `TRUE` the expectation is returned,
#'   in closed form (see Details). Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`. The method
#'   returns the parameter scale and the generic applies the link scale to its
#'   result.
#' @param approx Not read: the expectation is in closed form. It is present so
#'   that the signature matches the generic's.
#' @param nsim Not read. It is present so that the signature matches the
#'   generic's.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @return A named list of numeric vectors of length \eqn{n}, keyed and ordered
#'   as `deriv_names(distrib@params, 3)`. At \eqn{p = 3} there are ten
#'   components. With `expected = TRUE` every vector is constant.
#'
#' @section Notation:
#' \eqn{\alpha} is the shape vector, \eqn{\phi} the concentration, \eqn{\mu}
#' the mean on the simplex, \eqn{\eta} its free vector, \eqn{p} the number of
#' coordinates, \eqn{\Gamma} the gamma function, \eqn{\psi} the digamma
#' function and \eqn{\ell} the log-density of one observation.
#'
#' @seealso [distrib_deriv4.DirichletDistrib()] for the next order,
#'   [distrib_hessian.DirichletDistrib()] for the second,
#'   [dirichlet_higher()] for the shared engine, and [distrib_deriv3()] for the
#'   generic.
#'
#' @examples
#' d <- dirichlet_distrib(3)
#' theta <- list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8)
#' set.seed(1)
#' y <- distrib_rng(d, 5, theta)
#'
#' d3 <- distrib_deriv3(d, y, theta)
#' names(d3)
#'
#' # Against a central difference of the Hessian in phi.
#' h <- 1e-4
#' tp <- theta; tp$phi <- tp$phi + h
#' tm <- theta; tm$phi <- tm$phi - h
#' c(exact = sum(d3[["mean_alr1_mean_alr1_phi"]]),
#'   stencil = (sum(distrib_hessian(d, y, tp)[["mean_alr1_mean_alr1"]]) -
#'              sum(distrib_hessian(d, y, tm)[["mean_alr1_mean_alr1"]])) / (2 * h))
#'
#' # The pure-phi component carries no response: it is psi''(phi) minus the
#' # mu_j^3-weighted psi''(phi mu_j), the same at every observation.
#' mu <- as.numeric(parameters7::param_value(d@param, c(0.3, -0.2)))
#' c(component = d3[["phi_phi_phi"]][1],
#'   direct = psigamma(8, 2) - sum(mu^3 * psigamma(8 * mu, 2)))
#'
#' # The expectation is exact, so it does not depend on the seed.
#' distrib_deriv3(d, y, theta, expected = TRUE)[["phi_phi_phi"]][1]
S7::method(distrib_deriv3, DirichletDistrib) <- function(distrib, y, theta,
                                                         expected = FALSE,
                                                         scale = c("parameter", "link"),
                                                         approx = c("integrate", "bartlett", "mc", "opg"),
                                                         nsim = 10000, ...) {
  dirichlet_higher(distrib, y, theta, 3L, expected = expected)
}

#' @title Dirichlet Fourth Derivatives
#' @name distrib_deriv4.DirichletDistrib
#'
#' @description
#' Computes every fourth derivative of the log-density in the parameters, in
#' closed form, by the same construction as
#' [distrib_deriv3.DirichletDistrib()] one order up: a univariate partition sum
#' per coordinate of the simplex, run by [chain_univariate()] over the arrays
#' that [dirichlet_map_tensors()] supplies, plus the direct
#' \eqn{\log\Gamma(\phi)} term on the component whose indices all name the
#' concentration.
#'
#' @details
#' Run at order one, the same assembly reproduces
#' [distrib_gradient.DirichletDistrib()], which is written separately; on the
#' additive log-ratio chart [distrib_hessian.DirichletDistrib()] calls it at
#' order two. The split of \eqn{\log\Gamma(\alpha_j)} that keeps a small shape
#' accurate is described in [distrib_deriv3.DirichletDistrib()].
#'
#' With `expected = TRUE` the expectation is exact, as for the third order.
#'
#' @param distrib A `DirichletDistrib` object, from [dirichlet_distrib()].
#' @param y A numeric matrix with one row per observation and one column per
#'   coordinate, each row a point of the open simplex. A row with a zero
#'   coordinate is outside the support and gives a non-finite component through
#'   \eqn{\log y_j}.
#' @param theta A named list of parameters, each component a single number:
#'   `mean_alr1`, ..., `mean_alr(p-1)` and `phi`.
#' @param expected Logical of length 1. When `TRUE` the expectation is returned,
#'   in closed form (see Details). Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`. The method
#'   returns the parameter scale and the generic applies the link scale to its
#'   result.
#' @param approx Not read: the expectation is in closed form. It is present so
#'   that the signature matches the generic's.
#' @param nsim Not read. It is present so that the signature matches the
#'   generic's.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @return A named list of numeric vectors of length \eqn{n}, keyed and ordered
#'   as `deriv_names(distrib@params, 4)`. At \eqn{p = 3} there are fifteen
#'   components. With `expected = TRUE` every vector is constant.
#'
#' @section Notation:
#' \eqn{\alpha} is the shape vector, \eqn{\phi} the concentration, \eqn{\mu}
#' the mean on the simplex, \eqn{\eta} its free vector and \eqn{\ell} the
#' log-density of one observation.
#'
#' @seealso [distrib_deriv3.DirichletDistrib()] for the order below,
#'   [dirichlet_higher()] for the shared engine, and [distrib_deriv4()] for the
#'   generic.
#'
#' @examples
#' d <- dirichlet_distrib(3)
#' theta <- list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8)
#' set.seed(1)
#' y <- distrib_rng(d, 5, theta)
#'
#' d4 <- distrib_deriv4(d, y, theta)
#' length(d4)
#'
#' # Against a central difference of the third order in phi.
#' h <- 1e-4
#' tp <- theta; tp$phi <- tp$phi + h
#' tm <- theta; tm$phi <- tm$phi - h
#' c(exact = sum(d4[["mean_alr1_mean_alr1_phi_phi"]]),
#'   stencil = (sum(distrib_deriv3(d, y, tp)[["mean_alr1_mean_alr1_phi"]]) -
#'              sum(distrib_deriv3(d, y, tm)[["mean_alr1_mean_alr1_phi"]])) /
#'             (2 * h))
S7::method(distrib_deriv4, DirichletDistrib) <- function(distrib, y, theta,
                                                         expected = FALSE,
                                                         scale = c("parameter", "link"),
                                                         approx = c("integrate", "bartlett", "mc", "opg"),
                                                         nsim = 10000, ...) {
  dirichlet_higher(distrib, y, theta, 4L, expected = expected)
}

#' Assemble a Dirichlet's Higher Derivatives
#'
#' @description
#' Runs the univariate chain rule over the shape vector for every component of
#' the requested order and adds the terms that do not pass through it: the
#' response, which enters \eqn{\alpha_j} linearly, \eqn{\log\Gamma(\phi)}, and
#' on the additive log-ratio chart the derivatives of
#' \eqn{\sum_j \log\alpha_j} (see [distrib_deriv3.DirichletDistrib()]).
#'
#' @param distrib A `DirichletDistrib` object.
#' @param y A matrix with one row per observation.
#' @param theta A named list of parameters.
#' @param order The derivative order, 1 to 4.
#' @param expected Logical of length 1. When `TRUE` the response term is
#'   replaced by its expectation, through
#'   \eqn{\mathrm{E}\log y_j = \psi(\alpha_j) - \psi(\phi)}.
#'
#' @return A named list of numeric vectors of length `nrow(y)`, keyed as
#'   `deriv_names(distrib@params, order)`.
#'
#' @keywords internal
dirichlet_higher <- function(distrib, y, theta, order, expected = FALSE) {
  s <- distrib@param
  v <- mv_flat_theta(distrib, theta)
  k <- s@n_free
  eta <- v[seq_len(k)]
  phi <- v[[k + 1L]]
  alpha <- phi * as.numeric(parameters7::param_value(s, eta))

  y <- if (is.matrix(y)) y else matrix(y, nrow = 1L)
  n <- nrow(y)

  ud <- dirichlet_map_tensors(s, eta, phi, order)
  alr <- mn_is_alr(s)
  if (alr) {
    # -log Gamma(a) = log a - log Gamma(1 + a). On the additive log-ratio
    # chart sum_j log a_j = K log phi + sum_(j<K) eta_j - K LSE(eta), whose
    # derivatives are added below in closed form, and log Gamma(1 + a) is
    # regular at a tiny shape: the chain rule on -log Gamma(a) itself cancels
    # there, terms of order one leaving a result of the order of a.
    prob <- exp(mn_log_probs(s, eta))
    nk <- length(alpha)
    fd <- lapply(seq_len(order), function(m) -psigamma(1 + alpha, m - 1L))
    uds <- ud
  } else {
    # f(a) = -log Gamma(a), so f^(m)(a) = -psi^(m-1)(a). Each coordinate is
    # scaled by sc = min(a, 1): the chain rule multiplies f^(m) by m tensors
    # of the map, each proportional to a, so f^(m) sc^m and the tensors over
    # sc give the same product without the overflow of psi^(m-1)(a) ~ a^-m at
    # a tiny shape. a^m psi^(m-1)(a) = (-1)^m (m-1)! + a^m psi^(m-1)(1 + a).
    sc <- pmin(alpha, 1)
    small <- alpha < 1
    fd <- lapply(seq_len(order), function(m) {
      out <- numeric(length(alpha))
      out[!small] <- -psigamma(alpha[!small], m - 1L)
      out[small] <- -((-1)^m * factorial(m - 1L) +
                        alpha[small]^m * psigamma(1 + alpha[small], m - 1L))
      out
    })
    uds <- lapply(ud, function(u) u / sc)
  }

  # the response enters only through sum_j (alpha_j - 1) log y_j, linearly;
  # its expectation is E log y_j = psi(alpha_j) - psi(phi)
  ely <- if (expected) digamma(alpha) - digamma(phi) else NULL
  ly <- if (expected) NULL else log(y)

  idx <- deriv_indices(distrib@params, order)
  out <- lapply(idx, function(t) {
    key <- paste(sort(t), collapse = ",")
    # the log-Gamma of the concentration enters only the component all of
    # whose indices name it, and neither term below carries the response
    nphi <- sum(t == k + 1L)
    direct <- if (nphi == order) psigamma(phi, order - 1L) else 0
    if (alr && nphi == order) {
      direct <- direct + nk * (-1)^(order - 1L) * factorial(order - 1L) / phi^order
    }
    if (alr && nphi == 0L) {
      direct <- direct + if (order == 1L) 1 - nk * prob[t] else
        -nk * mn_cumulant(prob, t)
    }
    base <- sum(chain_univariate(t, fd, uds)) + direct
    # the response enters alpha linearly, so it multiplies the map's own
    # tensor and is the only term that varies by observation
    if (expected) {
      rep(base + sum(ely * ud[[key]]), n)
    } else {
      base + as.numeric(ly %*% ud[[key]])
    }
  })
  stats::setNames(out, deriv_names(distrib@params, order))
}


#' @title Multinomial Third Derivatives
#' @name distrib_deriv3.MultinomialDistrib
#'
#' @description
#' Computes every third derivative of the log-mass in the parameters, in closed
#' form. Up to a constant in the counts the log-mass is
#' \eqn{\sum_j y_j \log p_j}. On the additive log-ratio chart of
#' [parameters7::simplex()] this equals
#' \eqn{\sum_{j<p} y_j\eta_j - N\,\mathrm{LSE}(\eta)}, with \eqn{N} the total of
#' the row, so every derivative of order two or more is \eqn{-N} times a joint
#' cumulant of the indicators of the categorical law (see [mn_cumulant()]). The
#' cumulant is a sum of products of probabilities, computed without dividing by
#' one, so it stays finite when a probability underflows. On another chart
#' each term depends on one coordinate of the simplex and the chain rule
#' reduces to a univariate partition sum per coordinate, with
#' \deqn{f^{(m)}(p) = (-1)^{m-1}(m-1)!\,p^{-m},}
#' weighted by the counts.
#'
#' @details
#' With `expected = TRUE` the expectation is exact. The derivatives are linear
#' in the counts, so the method replaces each count by its expectation
#' \eqn{n p_j}. On the additive log-ratio chart this replaces \eqn{N} by
#' \eqn{n}, and the result equals the observed value at every row on the
#' support.
#'
#' @param distrib A `MultinomialDistrib` object, from [multinomial_distrib()].
#' @param y A numeric matrix with one row per observation and one column per
#'   category, each row a vector of counts summing to `distrib@size`.
#' @param theta A named list of parameters, each component a single number:
#'   `probs_alr1`, ..., `probs_alr(p-1)`.
#' @param expected Logical of length 1. When `TRUE` the expectation is returned,
#'   in closed form (see Details). Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`. The method
#'   returns the parameter scale and the generic applies the link scale to its
#'   result.
#' @param approx Not read: the expectation is in closed form. It is present so
#'   that the signature matches the generic's.
#' @param nsim Not read. It is present so that the signature matches the
#'   generic's.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @return A named list of numeric vectors of length \eqn{n}, keyed and ordered
#'   as `deriv_names(distrib@params, 3)`. At \eqn{p = 3} there are four
#'   components, the simplex spending \eqn{p - 1} free values. With
#'   `expected = TRUE` every vector is constant.
#'
#' @section Notation:
#' \eqn{p_j} is the probability of category \eqn{j}, \eqn{y_j} its count,
#' \eqn{n} the number of trials, \eqn{\eta} the simplex's free vector and
#' \eqn{\ell} the log-mass of one observation.
#'
#' @seealso [distrib_deriv4.MultinomialDistrib()] for the next order,
#'   [multinomial_higher()] for the shared engine, [mn_cumulant()] for the
#'   cumulants, and [distrib_deriv3()] for the generic.
#'
#' @examples
#' d <- multinomial_distrib(3, size = 10)
#' theta <- list(probs_alr1 = 0.4, probs_alr2 = -0.3)
#' set.seed(2)
#' y <- distrib_rng(d, 5, theta)
#'
#' d3 <- distrib_deriv3(d, y, theta)
#' vapply(d3, sum, numeric(1))
#'
#' # Against a central difference of the Hessian.
#' h <- 1e-4
#' tp <- theta; tp$probs_alr2 <- tp$probs_alr2 + h
#' tm <- theta; tm$probs_alr2 <- tm$probs_alr2 - h
#' c(exact = sum(d3[["probs_alr1_probs_alr1_probs_alr2"]]),
#'   stencil = (sum(distrib_hessian(d, y, tp)[["probs_alr1_probs_alr1"]]) -
#'              sum(distrib_hessian(d, y, tm)[["probs_alr1_probs_alr1"]])) /
#'             (2 * h))
#'
#' # The expectation equals the sum over the finite support.
#' sup <- mv_support(d, theta)
#' w <- distrib_pdf(d, sup, theta)
#' c(expected = distrib_deriv3(d, y, theta, expected = TRUE)[[1]][1],
#'   sum = sum(w * distrib_deriv3(d, sup, theta)[[1]]))
S7::method(distrib_deriv3, MultinomialDistrib) <- function(distrib, y, theta,
                                                            expected = FALSE,
                                                            scale = c("parameter", "link"),
                                                            approx = c("integrate", "bartlett", "mc", "opg"),
                                                            nsim = 10000, ...) {
  multinomial_higher(distrib, y, theta, 3L, expected = expected)
}

#' @title Multinomial Fourth Derivatives
#' @name distrib_deriv4.MultinomialDistrib
#'
#' @description
#' Computes every fourth derivative of the log-mass in the parameters, in
#' closed form, by the same construction as
#' [distrib_deriv3.MultinomialDistrib()] one order up: \eqn{-N} times a cumulant
#' of the categorical law on the additive log-ratio chart, and a univariate
#' partition sum per coordinate of the simplex, with
#' \eqn{f^{(m)}(p) = (-1)^{m-1}(m-1)!\,p^{-m}} and the counts as weights, on
#' another chart.
#'
#' @details
#' Run at orders one and two, the same assembly reproduces
#' [distrib_gradient.MultinomialDistrib()] and
#' [distrib_hessian.MultinomialDistrib()], which are written separately.
#'
#' With `expected = TRUE` the expectation is exact, as for the third order.
#'
#' @param distrib A `MultinomialDistrib` object, from [multinomial_distrib()].
#' @param y A numeric matrix with one row per observation and one column per
#'   category, each row a vector of counts summing to `distrib@size`.
#' @param theta A named list of parameters, each component a single number:
#'   `probs_alr1`, ..., `probs_alr(p-1)`.
#' @param expected Logical of length 1. When `TRUE` the expectation is returned,
#'   in closed form (see Details). Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`. The method
#'   returns the parameter scale and the generic applies the link scale to its
#'   result.
#' @param approx Not read: the expectation is in closed form. It is present so
#'   that the signature matches the generic's.
#' @param nsim Not read. It is present so that the signature matches the
#'   generic's.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @return A named list of numeric vectors of length \eqn{n}, keyed and ordered
#'   as `deriv_names(distrib@params, 4)`. At \eqn{p = 3} there are five
#'   components. With `expected = TRUE` every vector is constant.
#'
#' @section Notation:
#' \eqn{p_j} is the probability of category \eqn{j}, \eqn{y_j} its count,
#' \eqn{N} the total of a row, \eqn{\eta} the simplex's free vector and
#' \eqn{\ell} the log-mass of one observation.
#'
#' @seealso [distrib_deriv3.MultinomialDistrib()] for the order below,
#'   [multinomial_higher()] for the shared engine, and [distrib_deriv4()] for
#'   the generic.
#'
#' @examples
#' d <- multinomial_distrib(3, size = 10)
#' theta <- list(probs_alr1 = 0.4, probs_alr2 = -0.3)
#' set.seed(2)
#' y <- distrib_rng(d, 5, theta)
#'
#' d4 <- distrib_deriv4(d, y, theta)
#' names(d4)
#'
#' # Against a central difference of the third order.
#' h <- 1e-4
#' tp <- theta; tp$probs_alr2 <- tp$probs_alr2 + h
#' tm <- theta; tm$probs_alr2 <- tm$probs_alr2 - h
#' c(exact = sum(d4[["probs_alr1_probs_alr1_probs_alr2_probs_alr2"]]),
#'   stencil = (sum(distrib_deriv3(d, y, tp)[["probs_alr1_probs_alr1_probs_alr2"]]) -
#'              sum(distrib_deriv3(d, y, tm)[["probs_alr1_probs_alr1_probs_alr2"]])) /
#'             (2 * h))
S7::method(distrib_deriv4, MultinomialDistrib) <- function(distrib, y, theta,
                                                            expected = FALSE,
                                                            scale = c("parameter", "link"),
                                                            approx = c("integrate", "bartlett", "mc", "opg"),
                                                            nsim = 10000, ...) {
  multinomial_higher(distrib, y, theta, 4L, expected = expected)
}

#' @title Whether a Simplex Parameter Is the Additive Log-Ratio Chart
#'
#' @description
#' Returns `TRUE` when the probability parameter is a [parameters7::simplex()],
#' whose free values are the additive log-ratios against the last category.
#' On that chart the log-probabilities are \eqn{\eta_j - \mathrm{LSE}(\eta)}
#' and their derivatives are cumulants of the categorical law, which
#' [mn_log_probs()] and [mn_cumulant()] compute without dividing by a
#' probability.
#'
#' @param s The probability parameter of a [multinomial_distrib()] or a
#'   [dirichlet_distrib()].
#'
#' @return A single logical.
#'
#' @seealso [multinomial_higher()] and [distrib_hessian.MultinomialDistrib()].
#'
#' @examples
#' distributions7:::mn_is_alr(parameters7::simplex(3))
#'
#' @keywords internal
mn_is_alr <- function(s) {
  S7::S7_inherits(s, parameters7::SimplexParam)
}


#' @title The Log-Probabilities on the Additive Log-Ratio Chart
#'
#' @description
#' Computes \eqn{\log p_j = \eta_j - \mathrm{LSE}(\eta_1, \ldots,
#' \eta_{K-1}, 0)}, with \eqn{\eta_K = 0}, the log-sum-exp shifted by its
#' largest argument. A probability that underflows keeps its logarithm.
#'
#' @param s A [parameters7::simplex()] parameter.
#' @param eta Its free vector.
#'
#' @return A numeric vector of length `s@n_free + 1`.
#'
#' @seealso [mn_is_alr()] and [mn_cumulant()].
#'
#' @examples
#' distributions7:::mn_log_probs(parameters7::simplex(3), c(-800, 0))
#'
#' @keywords internal
mn_log_probs <- function(s, eta) {
  e <- c(as.numeric(eta), 0)
  m <- max(e)
  e - (m + log(sum(exp(e - m))))
}


#' @title A Joint Cumulant of the Indicators of a Categorical Law
#'
#' @description
#' Computes the joint cumulant of the indicators \eqn{1\{Z = t_1\}, \ldots,
#' 1\{Z = t_r\}} of a categorical variable with probabilities `prob`, which is
#' the derivative of the log-sum-exp \eqn{\mathrm{LSE}(\eta)} in
#' \eqn{\eta_{t_1}, \ldots, \eta_{t_r}}. It is the sum over the set partitions
#' \eqn{\pi} of the positions of
#' \eqn{(-1)^{\lvert\pi\rvert-1}(\lvert\pi\rvert-1)!\prod_{B\in\pi}\mu_B},
#' where \eqn{\mu_B = p_k} if every index of \eqn{B} is \eqn{k} and zero
#' otherwise.
#'
#' @param prob The probability vector.
#' @param t An integer vector of category indices, of length 1 to 4.
#'
#' @return A single number.
#'
#' @seealso [multinomial_higher()], the consumer.
#'
#' @examples
#' p <- c(0.2, 0.3, 0.5)
#' # The second cumulant is p_k (delta_kl - p_l).
#' c(distributions7:::mn_cumulant(p, c(1, 1)), 0.2 * 0.8)
#' c(distributions7:::mn_cumulant(p, c(1, 2)), -0.2 * 0.3)
#'
#' @keywords internal
mn_cumulant <- function(prob, t) {
  acc <- 0
  for (pt in numericals7::set_partitions(length(t))) {
    q <- length(pt)
    term <- (-1)^(q - 1L) * factorial(q - 1L)
    for (b in pt) {
      v <- t[b]
      term <- if (all(v == v[1L])) term * prob[v[1L]] else 0
      if (term == 0) break
    }
    acc <- acc + term
  }
  acc
}


#' Assemble a Multinomial's Higher Derivatives
#'
#' @description
#' Computes every component of the requested order. On the additive log-ratio
#' chart a component of order one is \eqn{y_k - N p_k} and a component of
#' higher order is \eqn{-N} times a cumulant from [mn_cumulant()], with \eqn{N}
#' the total of the row. On another chart the method runs the univariate chain
#' rule over the probability vector for every component, weighted by the
#' counts of each category.
#'
#' @param distrib A `MultinomialDistrib` object.
#' @param y A matrix with one row per observation.
#' @param theta A named list of parameters.
#' @param order The derivative order, 1 to 4.
#' @param expected Logical of length 1. When `TRUE` each count is replaced by
#'   its expectation \eqn{n p_j}, which gives the expected derivatives exactly.
#'
#' @return A named list of numeric vectors of length `nrow(y)`, keyed as
#'   `deriv_names(distrib@params, order)`.
#'
#' @keywords internal
multinomial_higher <- function(distrib, y, theta, order, expected = FALSE) {
  s <- distrib@param
  eta <- mv_flat_theta(distrib, theta)
  y <- if (is.matrix(y)) y else matrix(y, nrow = 1L)
  n <- nrow(y)
  idx <- deriv_indices(distrib@params, order)

  if (mn_is_alr(s)) {
    # sum_j y_j log p_j = sum_(j<K) y_j eta_j - N LSE(eta), so every
    # derivative of order two or more is -N times a cumulant of the
    # categorical law, a sum of products of probabilities with no division
    prob <- exp(mn_log_probs(s, eta))
    tot <- if (expected) rep(distrib@size, n) else rowSums(y)
    out <- lapply(idx, function(t) {
      if (length(t) == 1L) {
        if (expected) rep(0, n) else y[, t] - tot * prob[t]
      } else {
        -tot * mn_cumulant(prob, t)
      }
    })
    return(stats::setNames(out, deriv_names(distrib@params, order)))
  }

  prob <- as.numeric(parameters7::param_value(s, eta))
  ud <- simplex_map_tensors(s, eta, order)
  # f(p) = log p, so f^(m)(p) = (-1)^(m-1) (m-1)! p^-m
  fd <- lapply(seq_len(order), function(m) {
    (-1)^(m - 1L) * factorial(m - 1L) / prob^m
  })
  # linear in the counts, whose expectation is size * p
  yy <- if (expected) matrix(distrib@size * prob, n, length(prob), byrow = TRUE) else y
  out <- lapply(idx, function(t) as.numeric(yy %*% chain_univariate(t, fd, ud)))
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
#' @param pget The accessor for \eqn{\partial^{t}\Sigma^{-1}}, the `get`
#'   element of the list that [mvg_ptensors()] returns.
#' @param p The dimension.
#'
#' @return A numeric vector of length `nrow(r)`.
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
#' Splits a component into its mean-and-matrix part and its \eqn{\nu} part. The
#' log-determinant contributes only to a component whose indices all name the
#' matrix, with the sign set by the orientation of the chart (see
#' [mvg_ptensors()]). A component whose indices all name \eqn{\nu} is computed
#' by [mvt_cp_d()] and [mvt_h_d()], which avoid the cancellation of the direct
#' forms at large \eqn{\nu}. Every other component is a univariate chain rule
#' in \eqn{q} over the set partitions of its mean and matrix indices, whose
#' outer derivatives \eqn{(-1)^{m-1}(m-1)!/(\nu+q)^{m}} are differentiated in
#' \eqn{\nu} by the Leibniz rule against the linear prefactor \eqn{(\nu+p)/2}.
#' The two terms of the Leibniz rule are combined into one fraction, since
#' written separately they agree to leading order as \eqn{\nu} grows.
#'
#' @param distrib A `MvStudentTDistrib` object.
#' @param y An \eqn{n \times p} matrix of observations.
#' @param theta A named list of parameters.
#' @param order The derivative order, 1 to 4.
#'
#' @return A named list of numeric vectors of length \eqn{n}, keyed as
#'   `deriv_names(distrib@params, order)`.
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
  q <- mv_quad(r, s, pc$eta, distrib@inverted)

  pt <- mvg_ptensors(pc, order, inverted = distrib@inverted)
  ldfun <- switch(order,
    parameters7::param_dlogdet, parameters7::param_d2logdet,
    parameters7::param_d3logdet, parameters7::param_d4logdet
  )
  ld <- ldfun(s, pc$eta)
  ldkey <- vapply(parameters7::param_tuple_indices(s, order),
                  function(t) paste(sort(t), collapse = ","), character(1))

  # (nu + p)/2 d^a/dnu^a (nu + q)^-m + a/2 d^(a-1)/dnu^(a-1) (nu + q)^-m,
  # combined into one fraction: written as two terms, they agree to leading
  # order as nu grows. With s = nu + q it is
  # (-1)^a (m+a-2)!/(m-1)! [nu (m-1) + p (m+a-1) - a q] / (2 s^(m+a)).
  is <- 1 / (nu + q)
  kpart <- function(m, a) {
    if (a == 0L) return(0.5 * ((nu + p) * is) * is^(m - 1L))
    br <- nu * (m - 1L) + p * (m + a - 1L) - a * q
    0.5 * (-1)^a * (factorial(m + a - 2L) / factorial(m - 1L)) *
      (br * is) * is^(m + a - 1L)
  }

  idx <- deriv_indices(distrib@params, order)
  out <- lapply(idx, function(t) {
    a <- sum(t == inu)
    rest <- t[t != inu]
    acc <- rep(0, n)

    # the log-determinant, present only when every index names the matrix
    if (a == 0L && length(rest) == order && all(rest > p)) {
      acc <- acc + pt$sign_ld *
        ld[[match(paste(sort(rest - p), collapse = ","), ldkey)]]
    }

    if (!length(rest)) {
      # every index names nu: c_p^(a) - h^(a), each without its cancellation
      acc <- acc + mvt_cp_d(nu, p, a) - mvt_h_d(q / nu, nu, p, a)
    } else {
      # -((nu+p)/2) log(1 + q/nu), differentiated in the rest and then in nu
      for (part in index_partitions(rest)) {
        m <- length(part)
        term <- rep(1, n)
        for (b in part) term <- term * mvt_q_deriv(b, r, pt$get, p)
        coef <- (-1)^(m - 1L) * factorial(m - 1L) * kpart(m, a)
        acc <- acc - coef * term
      }
    }
    acc
  })
  stats::setNames(out, deriv_names(distrib@params, order))
}


#' @title Multivariate Student t Third Derivatives
#' @name distrib_deriv3.MvStudentTDistrib
#'
#' @description
#' Computes every third derivative of the log-density in the parameters, in
#' closed form. The log-density carries `lgamma`, a logarithm and a quadratic
#' form, each elementary in \eqn{\nu}. The location-and-matrix part of a
#' component is handled by the array expansion that the gaussian uses through
#' [mvg_ptensors()], and the \eqn{\nu} part by differentiating
#' \eqn{(-1)^{m-1}(m-1)!/(\nu+q)^m} against the linear prefactor
#' \eqn{(\nu+p)/2}. The components in \eqn{\nu} alone come from [mvt_cp_d()]
#' and [mvt_h_d()].
#'
#' @details
#' Run at orders one and two, the same assembly reproduces
#' [distrib_gradient.MvStudentTDistrib()] and
#' [distrib_hessian.MvStudentTDistrib()], which are written separately, for
#' the scale matrix and for its inverse.
#'
#' For the gaussian, every component with three location indices is zero,
#' because the quadratic form is quadratic in \eqn{\mu}. Here the weight
#' depends on \eqn{q}, and \eqn{\ell^{(\mu_1\mu_1\mu_1)}} is an ordinary
#' non-zero number.
#'
#' With `expected = TRUE` the expectation is taken by sampling and carries
#' Monte Carlo error of order `nsim^(-1/2)`.
#'
#' @param distrib An [MvStudentTDistrib] object, from [mvstudent_t1_distrib()]
#'   or [mvstudent_t2_distrib()].
#' @param y An \eqn{n \times p} numeric matrix of observations. A vector of
#'   length \eqn{p} is read as a single observation.
#' @param theta A named list of parameters, each component a single number.
#' @param expected Logical of length 1. When `TRUE` the expectation is returned,
#'   by sampling: `nsim` draws are made and the observed components averaged.
#'   Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`. The method
#'   returns the parameter scale and the generic applies the link scale to its
#'   result. The two differ here, because `nu` carries a log link by default.
#' @param approx Not read: `expected = TRUE` always takes the expectation by
#'   sampling. It is present so that the signature matches the generic's.
#' @param nsim The number of draws used when `expected = TRUE`. Defaults to
#'   `10000`. Not read otherwise.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @return A named list of numeric vectors of length \eqn{n}, keyed and ordered
#'   as `deriv_names(distrib@params, 3)`. At \eqn{p = 2} on an unstructured
#'   scale matrix there are 56 components. With `expected = TRUE` every vector
#'   is constant.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\Sigma} the scale matrix, \eqn{\nu} the
#' degrees of freedom, \eqn{q} the squared Mahalanobis distance, \eqn{p} the
#' dimension and \eqn{\ell} the log-density of one observation.
#'
#' @seealso [distrib_deriv4.MvStudentTDistrib()] for the next order,
#'   [distrib_hessian.MvStudentTDistrib()] for the second, [mvt_higher()] for
#'   the shared engine, [distrib_deriv3.MvGaussianDistrib()] for the limiting
#'   family, and [distrib_deriv3()] for the generic.
#'
#' @examples
#' d <- mvstudent_t1_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)
#' set.seed(3)
#' y <- distrib_rng(d, 4, theta)
#'
#' d3 <- distrib_deriv3(d, y, theta)
#' length(d3)
#'
#' # A gaussian's three-location component is exactly zero and this one is not.
#' c(t = sum(d3[["mu1_mu1_mu1"]]),
#'   gaussian = sum(distrib_deriv3(mvgaussian1_distrib(2), y,
#'                                 theta[1:5])[["mu1_mu1_mu1"]]))
#'
#' # Against a central difference of the Hessian in nu.
#' h <- 1e-4
#' tp <- theta; tp$nu <- tp$nu + h
#' tm <- theta; tm$nu <- tm$nu - h
#' c(exact = sum(d3[["mu1_mu2_nu"]]),
#'   stencil = (sum(distrib_hessian(d, y, tp)[["mu1_mu2"]]) -
#'              sum(distrib_hessian(d, y, tm)[["mu1_mu2"]])) / (2 * h))
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

#' @title Multivariate Student t Fourth Derivatives
#' @name distrib_deriv4.MvStudentTDistrib
#'
#' @description
#' Computes every fourth derivative of the log-density in the parameters, in
#' closed form, by the same construction as
#' [distrib_deriv3.MvStudentTDistrib()] one order up: the location and matrix
#' part through the array expansion of [mvg_ptensors()], the \eqn{\nu} part by
#' differentiating \eqn{(-1)^{m-1}(m-1)!/(\nu+q)^m} against the linear
#' prefactor \eqn{(\nu+p)/2}, and the components in \eqn{\nu} alone through
#' [mvt_cp_d()] and [mvt_h_d()].
#'
#' @details
#' Run at orders one and two, the same assembly reproduces the closed-form
#' gradient and Hessian, which are written separately, for the scale matrix
#' and for its inverse.
#'
#' With `expected = TRUE` the expectation is taken by sampling and carries
#' Monte Carlo error of order `nsim^(-1/2)`.
#'
#' @param distrib An [MvStudentTDistrib] object, from [mvstudent_t1_distrib()]
#'   or [mvstudent_t2_distrib()].
#' @param y An \eqn{n \times p} numeric matrix of observations. A vector of
#'   length \eqn{p} is read as a single observation.
#' @param theta A named list of parameters, each component a single number.
#' @param expected Logical of length 1. When `TRUE` the expectation is returned,
#'   by sampling: `nsim` draws are made and the observed components averaged.
#'   Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`. The method
#'   returns the parameter scale and the generic applies the link scale to its
#'   result. The two differ here, because `nu` carries a log link by default.
#' @param approx Not read: `expected = TRUE` always takes the expectation by
#'   sampling. It is present so that the signature matches the generic's.
#' @param nsim The number of draws used when `expected = TRUE`. Defaults to
#'   `10000`. Not read otherwise.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @return A named list of numeric vectors of length \eqn{n}, keyed and ordered
#'   as `deriv_names(distrib@params, 4)`. At \eqn{p = 2} on an unstructured
#'   scale matrix there are 126 components. With `expected = TRUE` every vector
#'   is constant.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\Sigma} the scale matrix, \eqn{\nu} the
#' degrees of freedom, \eqn{q} the squared Mahalanobis distance, \eqn{p} the
#' dimension and \eqn{\ell} the log-density of one observation.
#'
#' @seealso [distrib_deriv3.MvStudentTDistrib()] for the order below,
#'   [mvt_higher()] for the shared engine, and [distrib_deriv4()] for the
#'   generic.
#'
#' @examples
#' d <- mvstudent_t1_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)
#' set.seed(3)
#' y <- distrib_rng(d, 4, theta)
#'
#' d4 <- distrib_deriv4(d, y, theta)
#' length(d4)
#'
#' # Against a central difference of the third order in nu.
#' h <- 1e-4
#' tp <- theta; tp$nu <- tp$nu + h
#' tm <- theta; tm$nu <- tm$nu - h
#' c(exact = sum(d4[["mu1_mu2_nu_nu"]]),
#'   stencil = (sum(distrib_deriv3(d, y, tp)[["mu1_mu2_nu"]]) -
#'              sum(distrib_deriv3(d, y, tm)[["mu1_mu2_nu"]])) / (2 * h))
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
