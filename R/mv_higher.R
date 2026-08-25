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
#' \pkg{parameters7}, whose contract carries orders three and four for a
#' `simplex`, so nothing is derived here; what the function adds is the keying,
#' which lets the chain rule in [chain_univariate()] look a block up by its
#' index multiset.
#'
#' @param s A `parameters7::simplex()` parametrization.
#' @param eta The free vector, a numeric vector of length `s@n_free`.
#' @param order The highest order required, a single whole number from 1 to 4.
#'   Above 4 the `switch` returns `NULL` and the call fails.
#'
#' @return A named list of numeric vectors, each as long as the simplex has
#'   coordinates, keyed `"1"`, `"2"`, `"1,1"`, `"1,2"` and so on, with the
#'   indices sorted and comma-separated. Its length is the number of index
#'   multisets of width 1 to `order`.
#'
#' @section Notation:
#' \eqn{\mu} is the point of the simplex the parametrization produces,
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
#' This is why the two simplex-valued families are cheap at orders three and
#' four. Their log-densities are sums of terms each depending on ONE coordinate,
#' \eqn{\log\Gamma(\alpha_j)} for the Dirichlet and \eqn{\log p_j} for the
#' multinomial, so the multivariate Faa di Bruno collapses to one univariate
#' partition sum per coordinate and no mixed derivative array is ever formed.
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
#' concentration. The shape vector \eqn{\alpha = \phi\,\mu(\eta)} is BILINEAR
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
#'   tuple over the composite index set with \eqn{\phi} LAST, so at
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
#' Enumerates the index multisets of exactly `order` indices drawn from
#' `1:d`, in the order `parameters7::param_tuple_indices()` uses. It takes a
#' COUNT rather than a parametrization, so a composite index set formed by
#' appending a coordinate can be enumerated without building an object of that
#' shape. The Dirichlet's set is the one that needs this: its simplex free
#' values followed by its concentration. The work is done by
#' `numericals7::tuple_indices()`, the toolkit's one copy of the enumeration.
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

#' @title Expected Higher Derivatives of a Simplex-Valued Family by Sampling
#'
#' @description
#' Averages the closed-form observed third or fourth derivatives over a sample
#' drawn from the family. This is the route the multivariate branch takes
#' throughout: a quadrature over the simplex has no counterpart to the
#' one-dimensional split at quantiles, and the Bartlett route would need the
#' score's own higher derivatives. The result therefore carries Monte Carlo
#' error of order `nsim^(-1/2)`, and a caller who needs it reproducible sets a
#' seed first.
#'
#' @details
#' Only the two simplex-valued families reach here. The branch is taken on
#' whether `distrib` inherits `DirichletDistrib`, everything else being routed
#' to [multinomial_higher()], so this function must not be called for any third
#' family.
#'
#' @param distrib A `DirichletDistrib` or `MultinomialDistrib` object.
#' @param y The observed response, read only for its number of rows.
#' @param theta A named list of parameters, each component a single number.
#' @param order The derivative order, `3L` or `4L`.
#' @param nsim The Monte Carlo sample size, a single positive whole number.
#'
#' @return A named list of numeric vectors of length `n_obs(distrib, y)`, each
#'   constant across observations, keyed as `deriv_names(distrib@params, order)`.
#'
#' @seealso [dirichlet_higher()] and [multinomial_higher()] for the observed
#'   derivatives it averages, and [distrib_deriv3.DirichletDistrib()] for the
#'   method that calls it.
#'
#' @examples
#' d <- dirichlet_distrib(3)
#' theta <- list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8)
#' set.seed(1)
#' y <- distrib_rng(d, 5, theta)
#'
#' set.seed(5)
#' e3 <- distributions7:::mv_expected_higher(d, y, theta, 3L, 3000)
#'
#' # One value per observation, and every vector constant.
#' lengths(e3)[1:3]
#' all(vapply(e3, function(z) diff(range(z)) == 0, TRUE))
#'
#' # Two seeds give two answers, this being a sample rather than an integral.
#' set.seed(6)
#' c(first = e3[[1]][1],
#'   second = distributions7:::mv_expected_higher(d, y, theta, 3L, 3000)[[1]][1])
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


#' @title Dirichlet Third Derivatives
#' @name distrib_deriv3.DirichletDistrib
#'
#' @description
#' Computes every third derivative of the log-density in the parameters, in
#' closed form. The log-density is
#' \deqn{\ell = \log\Gamma(\phi) - \sum_j \log\Gamma(\alpha_j)
#'   + \sum_j(\alpha_j - 1)\log y_j, \qquad \alpha = \phi\,\mu(\eta),}
#' so every term depends on ONE coordinate of the simplex and the chain rule
#' collapses to a univariate partition sum per coordinate, run by
#' [chain_univariate()] over the arrays [dirichlet_map_tensors()] supplies. The
#' concentration also enters directly, through \eqn{\log\Gamma(\phi)}, which
#' contributes only to the component all of whose indices name it.
#'
#' @details
#' The shape vector is bilinear in \eqn{\phi} and \eqn{\mu}, so a component
#' naming \eqn{\phi} twice or more gets nothing from the map, and only the
#' direct \eqn{\log\Gamma(\phi)} term survives there. That is what keeps the
#' order-3 and order-4 assemblies as cheap as the order-2 one.
#'
#' With `expected = TRUE` the expectation is taken by sampling and carries
#' Monte Carlo error of order `nsim^(-1/2)`.
#' @param distrib A `DirichletDistrib` object, from [dirichlet_distrib()].
#' @param y A numeric matrix with one row per observation and one column per
#'   coordinate, each row a point of the open simplex. A row with a zero
#'   coordinate is outside the support and gives a non-finite component through
#'   \eqn{\log y_j}.
#' @param theta A named list of parameters, each component a single number:
#'   `mean_alr1`, ..., `mean_alr(p-1)` and `phi`.
#' @param expected Logical of length 1. When `TRUE` the expectation is returned,
#'   by SAMPLING through [mv_expected_higher()]. Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, handled by the
#'   generic before dispatch.
#' @param approx Ignored: sampling is the only multivariate route to an
#'   expectation. Present so that the signature matches the generic's.
#' @param nsim The number of draws used when `expected = TRUE`. Defaults to
#'   `10000`. Ignored otherwise.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @return A named list of numeric vectors of length \eqn{n}, keyed and ordered
#'   as `deriv_names(distrib@params, 3)`. At \eqn{p = 3} there are ten
#'   components. With `expected = TRUE` every vector is constant.
#'
#' @section Notation:
#' \eqn{\alpha} is the shape vector, \eqn{\phi} the concentration, \eqn{\mu}
#' the mean on the simplex, \eqn{\eta} its free vector, \eqn{\Gamma} the gamma
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
#' # Against one stencil on the analytic Hessian, which shares no algebra with
#' # the partition sum.
#' h <- 1e-4
#' tp <- theta; tp$phi <- tp$phi + h
#' tm <- theta; tm$phi <- tm$phi - h
#' c(exact = sum(d3[["mean_alr1_mean_alr1_phi"]]),
#'   stencil = (sum(distrib_hessian(d, y, tp)[["mean_alr1_mean_alr1"]]) -
#'              sum(distrib_hessian(d, y, tm)[["mean_alr1_mean_alr1"]])) / (2 * h))
#'
#' # The pure-phi component gets nothing from the map beyond the shapes
#' # themselves: it is psi''(phi) minus the mu_j^3-weighted psi''(phi mu_j),
#' # and so is the same at every observation.
#' mu <- as.numeric(parameters7::param_value(d@param, c(0.3, -0.2)))
#' c(component = d3[["phi_phi_phi"]][1],
#'   direct = psigamma(8, 2) - sum(mu^3 * psigamma(8 * mu, 2)))
#'
S7::method(distrib_deriv3, DirichletDistrib) <- function(distrib, y, theta,
                                                         expected = FALSE,
                                                         scale = c("parameter", "link"),
                                                         approx = c("integrate", "bartlett", "mc", "opg"),
                                                         nsim = 10000, ...) {
  if (expected) return(mv_expected_higher(distrib, y, theta, 3L, nsim))
  dirichlet_higher(distrib, y, theta, 3L)
}

#' @title Dirichlet Fourth Derivatives
#' @name distrib_deriv4.DirichletDistrib
#'
#' @description
#' Computes every fourth derivative of the log-density in the parameters, in
#' closed form, by the same construction as
#' [distrib_deriv3.DirichletDistrib()] one order up: a univariate partition sum
#' per coordinate of the simplex, run by [chain_univariate()] over the arrays
#' [dirichlet_map_tensors()] supplies, plus the direct \eqn{\log\Gamma(\phi)}
#' term on the component all of whose indices name the concentration.
#'
#' @details
#' The license for this order is that the SAME assembly run at orders one and
#' two reproduces the hand-written score and information, which are derived
#' separately and are already under [check_distrib()]. Agreement at the orders
#' that can be checked is what authorizes the order that cannot.
#'
#' With `expected = TRUE` the expectation is taken by sampling and carries
#' Monte Carlo error of order `nsim^(-1/2)`.
#' @param distrib A `DirichletDistrib` object, from [dirichlet_distrib()].
#' @param y A numeric matrix with one row per observation and one column per
#'   coordinate, each row a point of the open simplex. A row with a zero
#'   coordinate is outside the support and gives a non-finite component through
#'   \eqn{\log y_j}.
#' @param theta A named list of parameters, each component a single number:
#'   `mean_alr1`, ..., `mean_alr(p-1)` and `phi`.
#' @param expected Logical of length 1. When `TRUE` the expectation is returned,
#'   by SAMPLING through [mv_expected_higher()]. Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, handled by the
#'   generic before dispatch.
#' @param approx Ignored: sampling is the only multivariate route to an
#'   expectation. Present so that the signature matches the generic's.
#' @param nsim The number of draws used when `expected = TRUE`. Defaults to
#'   `10000`. Ignored otherwise.
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
#' # Against one stencil on the analytic third order.
#' h <- 1e-4
#' tp <- theta; tp$phi <- tp$phi + h
#' tm <- theta; tm$phi <- tm$phi - h
#' c(exact = sum(d4[["mean_alr1_mean_alr1_phi_phi"]]),
#'   stencil = (sum(distrib_deriv3(d, y, tp)[["mean_alr1_mean_alr1_phi"]]) -
#'              sum(distrib_deriv3(d, y, tm)[["mean_alr1_mean_alr1_phi"]])) /
#'             (2 * h))
#'
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


#' @title Multinomial Third Derivatives
#' @name distrib_deriv3.MultinomialDistrib
#'
#' @description
#' Computes every third derivative of the log-mass in the parameters, in closed
#' form. Up to a constant in the counts the log-mass is
#' \eqn{\sum_j y_j \log p_j}, so each term depends on ONE coordinate of the
#' simplex and the chain rule collapses to a univariate partition sum per
#' coordinate, with
#' \deqn{f^{(m)}(p) = (-1)^{m-1}(m-1)!\,p^{-m}.}
#' The combinatorial factor carries no parameter, so it drops out of every
#' derivative and the counts enter only as weights.
#'
#' @details
#' The support is finite, so the expectation of any component is an exact sum
#' over `mv_support(distrib, theta)` and needs no sampling. `expected = TRUE`
#' nonetheless routes to [mv_expected_higher()], which draws; a caller who
#' wants the exact expectation forms the sum against
#' [distrib_pdf.MultinomialDistrib()] directly.
#' @param distrib A `MultinomialDistrib` object, from [multinomial_distrib()].
#' @param y A numeric matrix with one row per observation and one column per
#'   category, each row a vector of counts summing to `distrib@size`.
#' @param theta A named list of parameters, each component a single number:
#'   `probs_alr1`, ..., `probs_alr(p-1)`.
#' @param expected Logical of length 1. When `TRUE` the expectation is returned,
#'   by SAMPLING through [mv_expected_higher()]. Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, handled by the
#'   generic before dispatch.
#' @param approx Ignored: sampling is the only multivariate route to an
#'   expectation. Present so that the signature matches the generic's.
#' @param nsim The number of draws used when `expected = TRUE`. Defaults to
#'   `10000`. Ignored otherwise.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @return A named list of numeric vectors of length \eqn{n}, keyed and ordered
#'   as `deriv_names(distrib@params, 3)`. At \eqn{p = 3} there are four
#'   components, the simplex spending \eqn{p - 1} free values. With
#'   `expected = TRUE` every vector is constant.
#'
#' @section Notation:
#' \eqn{p_j} is the probability of category \eqn{j}, \eqn{y_j} its count,
#' \eqn{\eta} the simplex's free vector and \eqn{\ell} the log-mass of one
#' observation.
#'
#' @seealso [distrib_deriv4.MultinomialDistrib()] for the next order,
#'   [multinomial_higher()] for the shared engine, [mv_support()] for the exact
#'   route to an expectation, and [distrib_deriv3()] for the generic.
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
#' # Against one stencil on the analytic Hessian, component by component.
#' h <- 1e-4
#' tp <- theta; tp$probs_alr2 <- tp$probs_alr2 + h
#' tm <- theta; tm$probs_alr2 <- tm$probs_alr2 - h
#' c(exact = sum(d3[["probs_alr1_probs_alr1_probs_alr2"]]),
#'   stencil = (sum(distrib_hessian(d, y, tp)[["probs_alr1_probs_alr1"]]) -
#'              sum(distrib_hessian(d, y, tm)[["probs_alr1_probs_alr1"]])) /
#'             (2 * h))
#'
#' # The support is finite, so the mass over it sums to one exactly and an
#' # expectation can be taken as a sum rather than a sample.
#' sup <- mv_support(d, theta)
#' c(states = nrow(sup), mass = sum(distrib_pdf(d, sup, theta)))
#'
S7::method(distrib_deriv3, MultinomialDistrib) <- function(distrib, y, theta,
                                                            expected = FALSE,
                                                            scale = c("parameter", "link"),
                                                            approx = c("integrate", "bartlett", "mc", "opg"),
                                                            nsim = 10000, ...) {
  if (expected) return(mv_expected_higher(distrib, y, theta, 3L, nsim))
  multinomial_higher(distrib, y, theta, 3L)
}

#' @title Multinomial Fourth Derivatives
#' @name distrib_deriv4.MultinomialDistrib
#'
#' @description
#' Computes every fourth derivative of the log-mass in the parameters, in
#' closed form, by the same construction as
#' [distrib_deriv3.MultinomialDistrib()] one order up: a univariate partition
#' sum per coordinate of the simplex, with
#' \eqn{f^{(m)}(p) = (-1)^{m-1}(m-1)!\,p^{-m}} and the counts as weights.
#'
#' @details
#' The license for this order is that the SAME assembly run at orders one and
#' two reproduces the hand-written score and information, which are derived
#' separately and are already under [check_distrib()].
#' @param distrib A `MultinomialDistrib` object, from [multinomial_distrib()].
#' @param y A numeric matrix with one row per observation and one column per
#'   category, each row a vector of counts summing to `distrib@size`.
#' @param theta A named list of parameters, each component a single number:
#'   `probs_alr1`, ..., `probs_alr(p-1)`.
#' @param expected Logical of length 1. When `TRUE` the expectation is returned,
#'   by SAMPLING through [mv_expected_higher()]. Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, handled by the
#'   generic before dispatch.
#' @param approx Ignored: sampling is the only multivariate route to an
#'   expectation. Present so that the signature matches the generic's.
#' @param nsim The number of draws used when `expected = TRUE`. Defaults to
#'   `10000`. Ignored otherwise.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @return A named list of numeric vectors of length \eqn{n}, keyed and ordered
#'   as `deriv_names(distrib@params, 4)`. At \eqn{p = 3} there are five
#'   components. With `expected = TRUE` every vector is constant.
#'
#' @section Notation:
#' \eqn{p_j} is the probability of category \eqn{j}, \eqn{y_j} its count,
#' \eqn{\eta} the simplex's free vector and \eqn{\ell} the log-mass of one
#' observation.
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
#' # Against one stencil on the analytic third order.
#' h <- 1e-4
#' tp <- theta; tp$probs_alr2 <- tp$probs_alr2 + h
#' tm <- theta; tm$probs_alr2 <- tm$probs_alr2 - h
#' c(exact = sum(d4[["probs_alr1_probs_alr1_probs_alr2_probs_alr2"]]),
#'   stencil = (sum(distrib_deriv3(d, y, tp)[["probs_alr1_probs_alr1_probs_alr2"]]) -
#'              sum(distrib_deriv3(d, y, tm)[["probs_alr1_probs_alr1_probs_alr2"]])) /
#'             (2 * h))
#'
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


#' @title Multivariate Student t Third Derivatives
#' @name distrib_deriv3.MvStudentTDistrib
#'
#' @description
#' Computes every third derivative of the log-density in the parameters, in
#' closed form. Nothing is obstructed here: the log-density carries
#' `lgamma`, a logarithm and a quadratic form, each elementary in \eqn{\nu}, so
#' there is no distribution function to differentiate in its degrees of
#' freedom. The component splits into a location-and-matrix part, handled by
#' the same array expansion the gaussian uses through [mvg_ptensors()], and a
#' \eqn{\nu} part obtained by differentiating
#' \eqn{(-1)^{m-1}(m-1)!/(\nu+q)^m} against the linear prefactor
#' \eqn{(\nu+p)/2}.
#'
#' @details
#' The license for this order is that the SAME assembly run at orders one and
#' two reproduces the hand-written score and information, which are derived
#' separately and are already under [check_distrib()]. Against one stencil on
#' the analytic Hessian, every component here agrees to a relative
#' \eqn{1.2\times10^{-8}} or better.
#'
#' Unlike the gaussian's, no component vanishes. There a tuple with three
#' location indices is exactly zero, the quadratic form being quadratic in
#' \eqn{\mu}; here the weight depends on \eqn{q}, so
#' \eqn{\ell^{(\mu_1\mu_1\mu_1)}} is an ordinary non-zero number.
#'
#' With `expected = TRUE` the expectation is taken by sampling and carries
#' Monte Carlo error of order `nsim^(-1/2)`.
#' @param distrib An [MvStudentTDistrib] object, from [mvstudent_t_distrib()].
#' @param y An \eqn{n \times p} numeric matrix of observations. A vector of
#'   length \eqn{p} is read as a single observation.
#' @param theta A named list of parameters, each component a single number.
#' @param expected Logical of length 1. When `TRUE` the expectation is returned,
#'   by SAMPLING: `nsim` draws are made and the observed components averaged.
#'   Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, handled by the
#'   generic before dispatch. The two differ here, `nu` carrying a log link by
#'   default.
#' @param approx Ignored: sampling is the only multivariate route to an
#'   expectation. Present so that the signature matches the generic's.
#' @param nsim The number of draws used when `expected = TRUE`. Defaults to
#'   `10000`. Ignored otherwise.
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
#' d <- mvstudent_t_distrib(2)
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
#'   gaussian = sum(distrib_deriv3(mvgaussian_distrib(2), y,
#'                                 theta[1:5])[["mu1_mu1_mu1"]]))
#'
#' # Against one stencil on the analytic Hessian.
#' h <- 1e-4
#' tp <- theta; tp$nu <- tp$nu + h
#' tm <- theta; tm$nu <- tm$nu - h
#' c(exact = sum(d3[["mu1_mu2_nu"]]),
#'   stencil = (sum(distrib_hessian(d, y, tp)[["mu1_mu2"]]) -
#'              sum(distrib_hessian(d, y, tm)[["mu1_mu2"]])) / (2 * h))
#'
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
#' part through [mvg_ptensors()]'s array expansion, and the \eqn{\nu} part by
#' differentiating \eqn{(-1)^{m-1}(m-1)!/(\nu+q)^m} against the linear
#' prefactor \eqn{(\nu+p)/2}.
#'
#' @details
#' The license for this order is that the SAME assembly run at orders one and
#' two reproduces the hand-written score and information, which are derived
#' separately and are already under [check_distrib()].
#'
#' With `expected = TRUE` the expectation is taken by sampling and carries
#' Monte Carlo error of order `nsim^(-1/2)`.
#' @param distrib An [MvStudentTDistrib] object, from [mvstudent_t_distrib()].
#' @param y An \eqn{n \times p} numeric matrix of observations. A vector of
#'   length \eqn{p} is read as a single observation.
#' @param theta A named list of parameters, each component a single number.
#' @param expected Logical of length 1. When `TRUE` the expectation is returned,
#'   by SAMPLING: `nsim` draws are made and the observed components averaged.
#'   Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, handled by the
#'   generic before dispatch. The two differ here, `nu` carrying a log link by
#'   default.
#' @param approx Ignored: sampling is the only multivariate route to an
#'   expectation. Present so that the signature matches the generic's.
#' @param nsim The number of draws used when `expected = TRUE`. Defaults to
#'   `10000`. Ignored otherwise.
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
#' d <- mvstudent_t_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)
#' set.seed(3)
#' y <- distrib_rng(d, 4, theta)
#'
#' d4 <- distrib_deriv4(d, y, theta)
#' length(d4)
#'
#' # Against one stencil on the analytic third order.
#' h <- 1e-4
#' tp <- theta; tp$nu <- tp$nu + h
#' tm <- theta; tm$nu <- tm$nu - h
#' c(exact = sum(d4[["mu1_mu2_nu_nu"]]),
#'   stencil = (sum(distrib_deriv3(d, y, tp)[["mu1_mu2_nu"]]) -
#'              sum(distrib_deriv3(d, y, tm)[["mu1_mu2_nu"]])) / (2 * h))
#'
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
