#' @include distrib.R generics.R moments.R
NULL

#' @title S7 Classes for a Reparametrized Distribution
#' @name ReparamContinuousDistrib
#'
#' @description
#' What [reparametrize()] returns: the same law as its parent, written in
#' different coordinates. There is one class per kind of parent, so that a
#' continuous parent keeps the defaults registered on [continuous_distrib()]
#' and a discrete one those of [discrete_distrib()]. Nothing else distinguishes
#' the two: every method in this file is registered on both, from the same
#' body.
#'
#' @details
#' The three properties beyond a distribution's own are what the methods read.
#' `parent_distrib` is the law being rewritten and is delegated to for the
#' density, the distribution function, the quantile function, the generator and
#' the response derivatives, all of which a change of parametrization leaves
#' alone. `reparam_map` carries the new parameters to the parent's, and
#' `reparam_derivs` carries the map's partial derivatives, which is where the
#' exactness of the parameter derivatives comes from.
#'
#' @param distrib_name The name of the family, a single string.
#' @param dimension `"univariate"` or `"multivariate"`. Always the first for a
#'   reparametrization, [reparametrize()] refusing a multivariate parent.
#' @param params A character vector naming the new parameters.
#' @param params_bounds A named list of length-two numeric vectors, the open
#'   interval each new parameter lives in.
#' @param link_params A named list of \pkg{linkfunctions7} links, one per new
#'   parameter. Note the name: it is `link_params` on the object and `links` in
#'   [reparametrize()]'s signature.
#' @param params_interpretation A named character vector describing each new
#'   parameter.
#' @param params_smooth A named logical vector saying which parameters the
#'   log-density is differentiable in.
#' @param bounds A length-two numeric vector, the support.
#' @param n_params The number of new parameters.
#' @param parent_distrib The distribution being rewritten.
#' @param reparam_map The map from the new parameters to the parent's, a
#'   function of one named list returning another.
#' @param reparam_derivs The function returning the map's keyed partial tables,
#'   as [reparam_tables()] consumes them.
#'
#' @return An object of class `ReparamContinuousDistrib` or
#'   `ReparamDiscreteDistrib`, carrying every property of a `distrib` plus the
#'   three above.
#'
#' @seealso [reparametrize()], the constructor that fills these in;
#'   [is_reparam()] to test for either class;
#'   [reparam_theta()] and [reparam_tables()], the two readers.
#'
#' @examples
#' # The class a continuous parent gives.
#' d <- reparametrize(
#'   gaussian1_distrib(),
#'   map = function(psi) list(mu = psi$mu, sigma = sqrt(psi$sigma2)),
#'   params = c("mu", "sigma2"),
#'   bounds = list(mu = c(-Inf, Inf), sigma2 = c(0, Inf)),
#'   links = list(mu = linkfunctions7::identity_link(),
#'                sigma2 = linkfunctions7::log_link())
#' )
#' class(d)
#' d@parent_distrib@distrib_name
ReparamContinuousDistrib <- S7::new_class("ReparamContinuousDistrib",
  parent = continuous_distrib,
  properties = list(
    parent_distrib = S7::class_any,
    reparam_map = S7::class_function,
    reparam_derivs = S7::class_function
  )
)

#' @rdname ReparamContinuousDistrib
#' @name ReparamDiscreteDistrib
ReparamDiscreteDistrib <- S7::new_class("ReparamDiscreteDistrib",
  parent = discrete_distrib,
  properties = list(
    parent_distrib = S7::class_any,
    reparam_map = S7::class_function,
    reparam_derivs = S7::class_function
  )
)

#' Is This a Reparametrized Distribution?
#'
#' @description
#' `TRUE` for either of the two wrapper classes, [ReparamContinuousDistrib()]
#' and [ReparamDiscreteDistrib()], and `FALSE` for anything else. The test
#' exists because the two classes have different parents, so no single
#' `S7_inherits()` call answers it.
#'
#' @param distrib An object inheriting from `distrib`, or anything else.
#'
#' @return A single logical.
#'
#' @seealso [reparametrize()], which builds them;
#'   [ReparamContinuousDistrib()] for the classes.
#'
#' @examples
#' d <- reparametrize(
#'   gaussian1_distrib(),
#'   map = function(psi) list(mu = psi$mu, sigma = sqrt(psi$sigma2)),
#'   params = c("mu", "sigma2"),
#'   bounds = list(mu = c(-Inf, Inf), sigma2 = c(0, Inf)),
#'   links = list(mu = linkfunctions7::identity_link(),
#'                sigma2 = linkfunctions7::log_link())
#' )
#' c(distributions7:::is_reparam(d),
#'   distributions7:::is_reparam(gaussian2_distrib()))
#'
#' @keywords internal
is_reparam <- function(distrib) {
  S7::S7_inherits(distrib, ReparamContinuousDistrib) ||
    S7::S7_inherits(distrib, ReparamDiscreteDistrib)
}

#' The Parent's Parameters at the New Ones
#'
#' @description
#' Runs the map on plain numbers and returns the parent's parameters in the
#' parent's own order, which is the form every probability function needs
#' before delegating. The new parameters are aligned and validated first, so a missing
#' or out-of-bounds component throws here and not several frames down.
#'
#' @details
#' A map that fails to return some parent parameter is caught with a message
#' naming the parameter and the parent's full list. The check is cheap and is
#' repeated on every call: a map is ordinary R and may return different names
#' at different points, so the construction cannot settle it once.
#'
#' @param distrib A reparametrized distribution.
#' @param theta A named list of the new parameters, on the new parameter scale.
#'   Components may be vectors.
#'
#' @return A named list of the parent's parameters, in the parent's order.
#'
#' @seealso [reparam_tables()] for the map's derivatives;
#'   [reparametrize()]; [ReparamContinuousDistrib()].
#'
#' @keywords internal
reparam_theta <- function(distrib, theta) {
  theta <- align_theta(distrib, theta)
  out <- distrib@reparam_map(theta[seq_len(distrib@n_params)])
  parent <- distrib@parent_distrib
  missing <- setdiff(parent@params, names(out))
  if (length(missing)) {
    stop(sprintf(
      "The map did not return '%s'. It must return every parameter of '%s': %s.",
      missing[1L], parent@distrib_name, paste(parent@params, collapse = ", ")
    ), call. = FALSE)
  }
  out[parent@params]
}

#' The Map Derivatives as Keyed Tables
#'
#' @description
#' Returns, for every parent parameter, the partial derivatives of that
#' component of the map with respect to the new parameters, keyed by the sorted
#' tuple of new-parameter positions: `"1"`, `"1,2"`, `"2,2,3,3"` and so on. A
#' missing key is an exact zero, so a map with many vanishing partials costs
#' nothing for them.
#'
#' @details
#' # Where the tables come from
#'
#' When the family supplies `map_derivs`, they are its hand-written closed
#' forms; the shipped second parametrizations all do, and the formulas live in
#' [reparam_map_derivs()]. Otherwise each needed partial comes from one
#' finite-difference stencil of [numericals7::fd_derivative()] applied to the
#' analytic map, a single stencil per order and never a chain of differences,
#' at the accuracy that construction carries: about \eqn{10^{-8}} at first
#' order, fading with the order. Exact tables are therefore the recommendation
#' for any family fitted in earnest.
#'
#' @param distrib A reparametrized distribution.
#' @param theta A named list of the new parameters, already aligned by
#'   [reparam_theta()]'s caller. Only the first `n_params` components are read.
#'
#' @return A list over the parent's parameters, each element a keyed list of
#'   that component's partial derivatives.
#'
#' @seealso [reparam_stencil_derivs()] for the numerical route;
#'   [reparam_map_derivs()] for the hand-written tables;
#'   [chain_derivatives()], the consumer.
#'
#' @keywords internal
reparam_tables <- function(distrib, theta) {
  q <- distrib@n_params
  psi <- theta[seq_len(q)]
  distrib@reparam_derivs(psi)
}

#' One Stencil Per Map Partial
#'
#' @description
#' The numerical fallback behind [reparam_tables()]: every partial of the map
#' is one central stencil of [numericals7::fd_derivative()] in each direction,
#' applied to the analytic map, iterated across **distinct** directions. A
#' cross-variable composition of two first differences is one mixed stencil and
#' not the same-variable nesting the toolkit forbids.
#'
#' @details
#' The accuracy is about \eqn{10^{-8}} at first order and fades with the order,
#' which is why a family fitted in earnest supplies `map_derivs` instead. What
#' the fallback buys is that a user's own reparametrization works immediately,
#' with no derivatives written: the same bargain
#' [distrib_grad_y.continuous_distrib()] offers a density-only family.
#'
#' @param map The map function, of a named list of the new parameters returning
#'   a named list of the parent's.
#' @param params A character vector naming the new parameters.
#' @param parent_params A character vector naming the parent's, which fixes the
#'   order of the returned list.
#'
#' @return A function of the new parameters, usable as an object's
#'   `reparam_derivs`, returning the keyed tables [reparam_tables()] describes.
#'
#' @seealso [reparam_tables()], which calls the result;
#'   [numericals7::fd_derivative()] for the stencil;
#'   [reparam_map_derivs()] for the exact alternative.
#'
#' @keywords internal
reparam_stencil_derivs <- function(map, params, parent_params) {
  q <- length(params)
  pp <- length(parent_params)
  function(psi) {
    n <- max(lengths(psi))
    rows <- if (n > 1L) seq_len(n) else 1L
    tabs <- lapply(seq_len(pp), function(i) list())
    tuples <- unlist(lapply(1:4, function(r) {
      g <- expand.grid(rep(list(seq_len(q)), r))
      keys <- apply(g, 1L, function(z) paste(sort(z), collapse = ","))
      unique(keys)
    }))
    for (key in tuples) {
      tup <- as.integer(strsplit(key, ",")[[1L]])
      counts <- tabulate(tup, nbins = q)
      vals <- vapply(rows, function(r) {
        base <- lapply(psi, function(v) v[[min(r, length(v))]])
        f <- function(x) {
          th <- base
          for (j in which(counts > 0L)) th[[j]] <- x[[match(j, which(counts > 0L))]]
          unlist(map(stats::setNames(th, params))[parent_params])
        }
        dirs <- which(counts > 0L)
        # iterate single stencils across distinct directions
        gcur <- f
        for (j in dirs) {
          ord <- counts[j]
          gprev <- gcur
          xj <- base[[j]]
          gcur <- local({
            jj <- j; oo <- ord; gp <- gprev; dd <- dirs
            function(x) {
              vapply(seq_len(pp), function(i) {
                numericals7::fd_derivative(function(z) {
                  xx <- x
                  xx[[match(jj, dd)]] <- z
                  gp(xx)[i]
                }, xj, order = oo)
              }, numeric(1))
            }
          })
        }
        gcur(lapply(dirs, function(j) base[[j]]))
      }, numeric(pp))
      vals <- matrix(vals, nrow = pp)
      for (i in seq_len(pp)) {
        if (any(abs(vals[i, ]) > 1e-12)) tabs[[i]][[key]] <- vals[i, ]
      }
    }
    tabs
  }
}

#' Parameter Derivatives of a Reparametrized Distribution
#'
#' @description
#' The body the four derivative methods share: aligns the new parameters, maps
#' them to the parent's, fetches the map's tables and hands all three to
#' [chain_derivatives()]. The `order` and `expected` arguments are what
#' separate the four registrations.
#'
#' @param distrib A reparametrized distribution.
#' @param y The response, a numeric vector.
#' @param theta A named list of the new parameters, on the new parameter scale.
#' @param order The derivative order, 1 to 4.
#' @param expected Should the expected derivatives be carried? A single
#'   logical, `FALSE` by default. `TRUE` is meaningful from order 2 up, the
#'   score having mean zero.
#'
#' @return A named list of numeric vectors, keyed as
#'   [`deriv_names(distrib@params, order)`][deriv_names].
#'
#' @seealso [chain_derivatives()], which does the work;
#'   [distrib_gradient.ReparamContinuousDistrib()] and its three siblings, the
#'   registrations; [reparametrize()].
#'
#' @keywords internal
reparam_chain <- function(distrib, y, theta, order, expected = FALSE) {
  theta <- align_theta(distrib, theta)
  n <- max(n_obs(distrib@parent_distrib, y), 1L)
  chain_derivatives(
    parent = distrib@parent_distrib,
    y = y,
    th_par = reparam_theta(distrib, theta),
    maps = reparam_tables(distrib, theta),
    new_params = distrib@params,
    order = order,
    expected = expected
  )
}

#' The Partition Sum Itself
#'
#' @description
#' Carries a parent's derivatives of the log-density into new coordinates,
#' given the keyed partial tables of the map. It is separated from
#' [reparam_chain()] so that a family written in its own right, not obtained
#' through [reparametrize()], can use the same machinery instead of a second
#' copy of it: the hand-written second parametrizations all do.
#'
#' @details
#' # The identity
#'
#' \deqn{\ell^{(I)}(\psi) = \sum_{\pi} \sum_{i_1 \dots i_{|\pi|}}
#'       \ell^{(i_1 \dots i_{|\pi|})}(\theta)
#'       \prod_{B \in \pi} \frac{\partial^{|B|}\theta_{i_B}}{\partial\psi_B},}
#' the sum over set partitions of the multi-index. Every order the parent has
#' in closed form therefore survives in closed form.
#'
#' # Under expectation
#'
#' Expectation is linear and the map deterministic, so the expected derivatives
#' obey the same formula, and the term carrying the parent's score drops
#' because the score has mean zero. That is why `D[[1]]` is left `NULL` when
#' `expected` is `TRUE`: a parent with an exact expected information gives an
#' exact one here.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density, \eqn{\theta} the parent's parameters,
#' \eqn{\psi} the new ones, \eqn{\pi} a set partition of the multi-index and
#' \eqn{B} one of its blocks.
#'
#' @param parent The distribution whose derivatives are being carried.
#' @param y The response, a numeric vector.
#' @param th_par The parent's parameters as plain numbers, from
#'   [reparam_theta()].
#' @param maps The keyed partial tables of the map, as [reparam_tables()]
#'   returns them.
#' @param new_params A character vector naming the new parameters, which names
#'   and orders the result.
#' @param order The derivative order, 1 to 4.
#' @param expected Should the expected derivatives be carried? A single
#'   logical, `FALSE` by default.
#'
#' @return A named list of numeric vectors, keyed as
#'   [`deriv_names(new_params, order)`][deriv_names].
#'
#' @seealso [chain_assemble()], the sum it delegates to;
#'   [reparam_chain()], its caller for a [reparametrize()] object;
#'   [chain_cdf_deriv()], the same construction on the distribution function.
#'
#' @keywords internal
chain_derivatives <- function(parent, y, th_par, maps, new_params, order,
                              expected = FALSE) {
  n <- max(n_obs(parent, y), 1L)

  # The parent's derivatives of every order up to the one asked for. Under
  # expectation the score contributes nothing, by the first Bartlett identity.
  D <- vector("list", order)
  for (k in seq_len(order)) {
    D[[k]] <- if (expected) {
      switch(k,
        NULL,
        distrib_expected_hessian(parent, y, th_par),
        distrib_deriv3(parent, y, th_par, expected = TRUE),
        distrib_deriv4(parent, y, th_par, expected = TRUE)
      )
    } else {
      switch(k,
        distrib_gradient(parent, y, th_par),
        distrib_hessian(parent, y, th_par),
        distrib_deriv3(parent, y, th_par),
        distrib_deriv4(parent, y, th_par)
      )
    }
  }

  chain_assemble(D, parent@params, maps, new_params, order, n)
}

#' Faa di Bruno Over Set Partitions, on Tables
#'
#' @description
#' The partition sum of [chain_derivatives()], taking the inner derivatives as
#' tables instead of fetching them from a distribution. With the two split, a
#' family whose derivatives are easiest in coordinates it does not expose can
#' carry them into the ones it does, without a second copy of the enumeration.
#'
#' @details
#' [enet_distrib()] is the case the split was made for: its derivatives are
#' compact in the two rates \eqn{(a, c)} and the family is parametrized by
#' \eqn{(\lambda, \alpha)}, which is a bilinear map away. The same function
#' also serves the location-scale cdf derivatives, whose inner coordinate is
#' the standardized quantile.
#'
#' The key into an inner table is **built** from the parameter names and never
#' parsed out of one, so a parameter whose own name contains an underscore is
#' safe. That is the shape of mistake the package records having made once, in
#' three separate fallbacks at once.
#'
#' @param D A list of length `order`. `D[[k]]` is the inner derivative table of
#'   order \eqn{k}, keyed as
#'   [`deriv_names(inner_params, k)`][deriv_names]. `D[[1]]` may be `NULL`
#'   where the first-order term is known to vanish.
#' @param inner_params A character vector naming the inner coordinates.
#' @param maps Per inner coordinate, a keyed table of the map's partials in the
#'   outer coordinates, as [reparam_tables()] returns them.
#' @param new_params A character vector naming the outer parameters.
#' @param order The derivative order, 1 to 4.
#' @param n The length to recycle a constant component to, so that every
#'   component comes back one value per observation.
#'
#' @return A named list of numeric vectors, keyed as
#'   [`deriv_names(new_params, order)`][deriv_names].
#'
#' @seealso [chain_derivatives()], its caller;
#'   [loc_scale_cdf_deriv_k()] and [enet_distrib()], the two other consumers;
#'   [numericals7::set_partitions()] for the enumeration.
#'
#' @keywords internal
chain_assemble <- function(D, inner_params, maps, new_params, order, n) {
  p <- length(inner_params)

  # h^i_B for one block, from the keyed table; a missing key is a zero
  hvec <- function(i, tup) {
    v <- maps[[i]][[paste(sort(tup), collapse = ",")]]
    if (is.null(v)) 0 else v
  }

  # The inner derivative at a multiset of its own indices, keyed by name.
  # The key is BUILT from the names, never parsed out of one, so a parameter
  # whose own name contains an underscore is safe.
  pvec <- function(k, ids) {
    D[[k]][[paste(inner_params[sort(ids)], collapse = "_")]]
  }

  parts <- numericals7::set_partitions(order)
  assign_grid <- lapply(seq_len(order), function(nb) {
    as.matrix(expand.grid(rep(list(seq_len(p)), nb)))
  })

  idx <- deriv_indices(new_params, order)
  nm <- deriv_names(new_params, order)

  out <- lapply(seq_along(nm), function(m) {
    I <- idx[[m]]
    acc <- 0
    for (pi in parts) {
      # The number of blocks IS the order of the parent derivative the term
      # carries: one block means the score and a high derivative of the map,
      # while `order` blocks mean the parent's own top derivative and the
      # Jacobian. Under expectation the one-block term drops, which is what
      # D[[1]] being NULL above arranges.
      nb <- length(pi)
      g <- assign_grid[[nb]]
      for (r in seq_len(nrow(g))) {
        ids <- as.integer(g[r, ])
        hprod <- 1
        for (b in seq_len(nb)) {
          hb <- hvec(ids[b], I[pi[[b]]])
          if (all(hb == 0)) {
            hprod <- 0
            break
          }
          hprod <- hprod * hb
        }
        if (all(hprod == 0)) next
        d <- pvec(nb, ids)
        if (is.null(d)) next
        acc <- acc + d * hprod
      }
    }
    if (length(acc) == 1L) rep(acc, length.out = n) else acc
  })
  stats::setNames(out, nm)
}


# --- S7 METHODS IMPLEMENTATION ---

#' @title Density of a Reparametrized Distribution
#' @name distrib_pdf.ReparamContinuousDistrib
#' @aliases distrib_pdf.ReparamDiscreteDistrib
#'
#' @description
#' The parent's density, evaluated at the mapped parameters. A change of
#' parametrization does not change the law, so nothing is recomputed here: the
#' new parameters go through [reparam_theta()] and the parent answers.
#'
#' @param distrib A reparametrized distribution, from [reparametrize()].
#' @param y A numeric vector of observations.
#' @param theta A named list of the new parameters, on the new parameter scale.
#' @param log Should the log-density be returned? A single logical, `FALSE` by
#'   default.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of densities, the length of `y` recycled against
#'   `theta`.
#'
#' @seealso [reparametrize()]; [reparam_theta()] for the map;
#'   [distrib_gradient.ReparamContinuousDistrib()], where the parametrization
#'   does enter.
#'
#' @examples
#' d <- reparametrize(
#'   gaussian1_distrib(),
#'   map = function(psi) list(mu = psi$mu, sigma = sqrt(psi$sigma2)),
#'   params = c("mu", "sigma2"),
#'   bounds = list(mu = c(-Inf, Inf), sigma2 = c(0, Inf)),
#'   links = list(mu = linkfunctions7::identity_link(),
#'                sigma2 = linkfunctions7::log_link())
#' )
#' th <- list(mu = 1, sigma2 = 4)
#'
#' # The same numbers as the family written out by hand.
#' all.equal(distrib_pdf(d, c(0, 1, 2), th),
#'           distrib_pdf(gaussian2_distrib(), c(0, 1, 2), th))
#'
#' @keywords internal
reparam_pdf <- function(distrib, y, theta, log = FALSE, ...) {
  distrib_pdf(distrib@parent_distrib, y, reparam_theta(distrib, theta), log = log)
}
S7::method(distrib_pdf, ReparamContinuousDistrib) <- reparam_pdf
S7::method(distrib_pdf, ReparamDiscreteDistrib) <- reparam_pdf

#' @title Distribution Function of a Reparametrized Distribution
#' @name distrib_cdf.ReparamContinuousDistrib
#' @aliases distrib_cdf.ReparamDiscreteDistrib
#'
#' @description
#' The parent's distribution function, evaluated at the mapped parameters. The
#' `lower.tail` and `log.p` arguments are passed through and applied by the
#' parent, so they are not applied twice.
#'
#' @param distrib A reparametrized distribution, from [reparametrize()].
#' @param q A numeric vector of quantiles.
#' @param theta A named list of the new parameters, on the new parameter scale.
#' @param lower.tail Should probabilities be \eqn{P(Y \le q)}? A single
#'   logical, `TRUE` by default.
#' @param log.p Should log-probabilities be returned? A single logical, `FALSE`
#'   by default.
#' @param ... Passed to the parent's method.
#'
#' @return A numeric vector of probabilities, the length of `q` recycled
#'   against `theta`.
#'
#' @seealso [reparametrize()];
#'   [distrib_grad_cdf.ReparamContinuousDistrib()] for its derivatives.
#'
#' @examples
#' d <- reparametrize(
#'   gaussian1_distrib(),
#'   map = function(psi) list(mu = psi$mu, sigma = sqrt(psi$sigma2)),
#'   params = c("mu", "sigma2"),
#'   bounds = list(mu = c(-Inf, Inf), sigma2 = c(0, Inf)),
#'   links = list(mu = linkfunctions7::identity_link(),
#'                sigma2 = linkfunctions7::log_link())
#' )
#' distrib_cdf(d, c(0, 1, 2), list(mu = 1, sigma2 = 4))
#'
#' @keywords internal
reparam_cdf <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE, ...) {
  distrib_cdf(distrib@parent_distrib, q, reparam_theta(distrib, theta),
              lower.tail = lower.tail, log.p = log.p, ...)
}
S7::method(distrib_cdf, ReparamContinuousDistrib) <- reparam_cdf
S7::method(distrib_cdf, ReparamDiscreteDistrib) <- reparam_cdf

#' @title Quantile Function of a Reparametrized Distribution
#' @name distrib_quantile.ReparamContinuousDistrib
#' @aliases distrib_quantile.ReparamDiscreteDistrib
#'
#' @description
#' The parent's quantile function, evaluated at the mapped parameters. Where
#' the parent inverts its distribution function numerically, so does this; the
#' map costs one evaluation and does not enter the search.
#'
#' @param distrib A reparametrized distribution, from [reparametrize()].
#' @param p A numeric vector of probabilities in \eqn{[0,1]}, or of
#'   log-probabilities when `log.p` is `TRUE`.
#' @param theta A named list of the new parameters, on the new parameter scale.
#' @param lower.tail Are the probabilities \eqn{P(Y \le q)}? A single logical,
#'   `TRUE` by default.
#' @param log.p Is `p` a log-probability? A single logical, `FALSE` by default.
#' @param ... Passed to the parent's method.
#'
#' @return A numeric vector of quantiles, the length of `p` recycled against
#'   `theta`.
#'
#' @seealso [reparametrize()];
#'   [distrib_cdf.ReparamContinuousDistrib()], which this inverts.
#'
#' @examples
#' d <- reparametrize(
#'   gaussian1_distrib(),
#'   map = function(psi) list(mu = psi$mu, sigma = sqrt(psi$sigma2)),
#'   params = c("mu", "sigma2"),
#'   bounds = list(mu = c(-Inf, Inf), sigma2 = c(0, Inf)),
#'   links = list(mu = linkfunctions7::identity_link(),
#'                sigma2 = linkfunctions7::log_link())
#' )
#' distrib_quantile(d, c(0.25, 0.5, 0.75), list(mu = 1, sigma2 = 4))
#'
#' @keywords internal
reparam_quantile <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE, ...) {
  distrib_quantile(distrib@parent_distrib, p, reparam_theta(distrib, theta),
                   lower.tail = lower.tail, log.p = log.p, ...)
}
S7::method(distrib_quantile, ReparamContinuousDistrib) <- reparam_quantile
S7::method(distrib_quantile, ReparamDiscreteDistrib) <- reparam_quantile

#' @title Random Generation from a Reparametrized Distribution
#' @name distrib_rng.ReparamContinuousDistrib
#' @aliases distrib_rng.ReparamDiscreteDistrib
#'
#' @description
#' The parent's generator, at the mapped parameters. The draws come from the
#' parent's own method and consume its random numbers, so a seed set before the
#' call gives the same sample the parent would give at the mapped values.
#'
#' @param distrib A reparametrized distribution, from [reparametrize()].
#' @param n The number of draws, a single non-negative whole number.
#' @param theta A named list of the new parameters, on the new parameter scale.
#'   A component of length `n` gives one draw per setting.
#'
#' @return A numeric vector of `n` draws.
#'
#' @seealso [reparametrize()];
#'   [distrib_quantile.ReparamContinuousDistrib()], which an inverse-transform
#'   parent uses.
#'
#' @examples
#' d <- reparametrize(
#'   gaussian1_distrib(),
#'   map = function(psi) list(mu = psi$mu, sigma = sqrt(psi$sigma2)),
#'   params = c("mu", "sigma2"),
#'   bounds = list(mu = c(-Inf, Inf), sigma2 = c(0, Inf)),
#'   links = list(mu = linkfunctions7::identity_link(),
#'                sigma2 = linkfunctions7::log_link())
#' )
#'
#' # The same draws the parent gives at the mapped parameters.
#' set.seed(1); a <- distrib_rng(d, 5, list(mu = 1, sigma2 = 4))
#' set.seed(1); b <- distrib_rng(gaussian1_distrib(), 5,
#'                               list(mu = 1, sigma = 2))
#' all.equal(a, b)
#'
#' @keywords internal
reparam_rng <- function(distrib, n, theta) {
  distrib_rng(distrib@parent_distrib, n, reparam_theta(distrib, theta))
}
S7::method(distrib_rng, ReparamContinuousDistrib) <- reparam_rng
S7::method(distrib_rng, ReparamDiscreteDistrib) <- reparam_rng

#' @title Gradient of a Reparametrized Distribution
#' @name distrib_gradient.ReparamContinuousDistrib
#' @aliases distrib_gradient.ReparamDiscreteDistrib
#'
#' @description
#' The parent's score carried by the Jacobian of the map,
#' \deqn{\ell^{(a)} = \sum_i \ell^{(i)}\,
#'       \frac{\partial\theta_i}{\partial\psi_a}.}
#' It is exact whenever the parent's score is, and the map's partials are
#' whatever [reparam_tables()] supplies: a hand-written table where the family
#' has one, and one stencil per partial otherwise.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density, \eqn{\theta} the parent's parameters and
#' \eqn{\psi} the new ones.
#'
#' @param distrib A reparametrized distribution, from [reparametrize()].
#' @param y The response, a numeric vector.
#' @param theta A named list of the new parameters, on the new parameter scale.
#' @param scale Either `"parameter"` (the default) or `"link"`. Handled by the
#'   generic, which applies the chain rule onto the link scale after this
#'   method has returned; the method itself always answers on the parameter
#'   scale.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of numeric vectors, one per new parameter, each the
#'   length of `y` recycled against `theta`.
#'
#' @seealso [chain_derivatives()] for the identity;
#'   [distrib_hessian.ReparamContinuousDistrib()] for the second order;
#'   [reparametrize()].
#'
#' @examples
#' d <- reparametrize(
#'   gaussian1_distrib(),
#'   map = function(psi) list(mu = psi$mu, sigma = sqrt(psi$sigma2)),
#'   params = c("mu", "sigma2"),
#'   bounds = list(mu = c(-Inf, Inf), sigma2 = c(0, Inf)),
#'   links = list(mu = linkfunctions7::identity_link(),
#'                sigma2 = linkfunctions7::log_link())
#' )
#' th <- list(mu = 1, sigma2 = 4)
#'
#' # The same numbers as the family written out by hand, to 1e-12: the map's
#' # partials are differenced here and written out there.
#' a <- distrib_gradient(d, c(0, 1, 2), th)
#' b <- distrib_gradient(gaussian2_distrib(), c(0, 1, 2), th)
#' max(abs(unlist(a[names(b)]) - unlist(b)))
#'
#' @keywords internal
reparam_gradient <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  reparam_chain(distrib, y, theta, 1L)
}
S7::method(distrib_gradient, ReparamContinuousDistrib) <- reparam_gradient
S7::method(distrib_gradient, ReparamDiscreteDistrib) <- reparam_gradient

#' @title Observed Hessian of a Reparametrized Distribution
#' @name distrib_hessian.ReparamContinuousDistrib
#' @aliases distrib_hessian.ReparamDiscreteDistrib
#'
#' @description
#' The second-order chain rule, which keeps the term in the parent's score as
#' well as the one in its Hessian:
#' \deqn{\ell^{(ab)} = \sum_{i,j} \ell^{(ij)}
#'       \frac{\partial\theta_i}{\partial\psi_a}
#'       \frac{\partial\theta_j}{\partial\psi_b}
#'       + \sum_i \ell^{(i)}
#'       \frac{\partial^2\theta_i}{\partial\psi_a \partial\psi_b}.}
#' The second sum is what a first-order chain does not have, and it is why the
#' map's second partials are needed as well as its first.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density, \eqn{\theta} the parent's parameters and
#' \eqn{\psi} the new ones.
#'
#' @param distrib A reparametrized distribution, from [reparametrize()].
#' @param y The response, a numeric vector.
#' @param theta A named list of the new parameters, on the new parameter scale.
#' @param scale Either `"parameter"` (the default) or `"link"`, handled by the
#'   generic after this method has returned.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of numeric vectors keyed as [hess_names()], each the
#'   length of `y` recycled against `theta`. Note the **order**: this route
#'   enumerates as [hess_names()] does, diagonal first, while a family written
#'   out by hand may enumerate lexicographically. Compare two such results by
#'   name and not by position.
#'
#' @seealso [distrib_gradient.ReparamContinuousDistrib()] for the first order;
#'   [distrib_expected_hessian.ReparamContinuousDistrib()], where the second
#'   sum drops; [chain_derivatives()].
#'
#' @examples
#' d <- reparametrize(
#'   gaussian1_distrib(),
#'   map = function(psi) list(mu = psi$mu, sigma = sqrt(psi$sigma2)),
#'   params = c("mu", "sigma2"),
#'   bounds = list(mu = c(-Inf, Inf), sigma2 = c(0, Inf)),
#'   links = list(mu = linkfunctions7::identity_link(),
#'                sigma2 = linkfunctions7::log_link())
#' )
#' th <- list(mu = 1, sigma2 = 4)
#'
#' a <- distrib_hessian(d, c(0, 1, 2), th)
#' b <- distrib_hessian(gaussian2_distrib(), c(0, 1, 2), th)
#'
#' # The two enumerate their components in different orders, so match by name.
#' rbind(reparametrized = names(a), written_out = names(b))
#' max(abs(unlist(a[names(b)]) - unlist(b)))
#'
#' @keywords internal
reparam_hessian <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  reparam_chain(distrib, y, theta, 2L)[hess_names(distrib@params)]
}
S7::method(distrib_hessian, ReparamContinuousDistrib) <- reparam_hessian
S7::method(distrib_hessian, ReparamDiscreteDistrib) <- reparam_hessian

#' @title Expected Information of a Reparametrized Distribution
#' @name distrib_expected_hessian.ReparamContinuousDistrib
#' @aliases distrib_expected_hessian.ReparamDiscreteDistrib
#'
#' @description
#' The same chain rule as the observed Hessian, with the term in the parent's
#' score dropped: expectation is linear and the map deterministic, and the
#' score has mean zero, so
#' \deqn{E[\ell^{(ab)}] = \sum_{i,j} E[\ell^{(ij)}]
#'       \frac{\partial\theta_i}{\partial\psi_a}
#'       \frac{\partial\theta_j}{\partial\psi_b}.}
#' A parent with an exact expected information therefore gives an exact one
#' here, and the map's second partials are not read at all.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density, \eqn{\theta} the parent's parameters and
#' \eqn{\psi} the new ones.
#'
#' @param distrib A reparametrized distribution, from [reparametrize()].
#' @param y The response, a numeric vector. Its length is what the result is
#'   recycled to; the values are read only where the parent's own expected
#'   information reads them.
#' @param theta A named list of the new parameters, on the new parameter scale.
#' @param scale Either `"parameter"` (the default) or `"link"`, handled by the
#'   generic after this method has returned.
#' @param ... Passed to the parent's method, which is where `approx` and
#'   `nsim` are read for a family that approximates its expected information.
#'
#' @return A named list of numeric vectors keyed as [hess_names()], each the
#'   length of `y` recycled against `theta`, in that enumeration's order.
#'
#' @seealso [distrib_hessian.ReparamContinuousDistrib()], which keeps the score
#'   term; [chain_derivatives()]; [reparametrize()].
#'
#' @examples
#' d <- reparametrize(
#'   gaussian1_distrib(),
#'   map = function(psi) list(mu = psi$mu, sigma = sqrt(psi$sigma2)),
#'   params = c("mu", "sigma2"),
#'   bounds = list(mu = c(-Inf, Inf), sigma2 = c(0, Inf)),
#'   links = list(mu = linkfunctions7::identity_link(),
#'                sigma2 = linkfunctions7::log_link())
#' )
#' th <- list(mu = 1, sigma2 = 4)
#'
#' a <- distrib_expected_hessian(d, c(0, 1, 2), th)
#' b <- distrib_expected_hessian(gaussian2_distrib(), c(0, 1, 2), th)
#'
#' # Matched by name, the two enumerating in different orders.
#' max(abs(unlist(a[names(b)]) - unlist(b)))
#'
#' @keywords internal
reparam_expected_hessian <- function(distrib, y, theta,
                                     scale = c("parameter", "link"),
                                     approx = c("bartlett", "integrate", "mc", "opg"),
                                     nsim = 10000, ...) {
  reparam_chain(distrib, y, theta, 2L, expected = TRUE)[hess_names(distrib@params)]
}
S7::method(distrib_expected_hessian, ReparamContinuousDistrib) <- reparam_expected_hessian
S7::method(distrib_expected_hessian, ReparamDiscreteDistrib) <- reparam_expected_hessian

#' @title Third-Order Derivatives of a Reparametrized Distribution
#' @name distrib_deriv3.ReparamContinuousDistrib
#' @aliases distrib_deriv3.ReparamDiscreteDistrib
#'
#' @description
#' The partition sum at order three, observed or expected. Every order the
#' parent has in closed form survives in closed form, so a family reparametrized
#' from one of the 42 that carry analytic third derivatives keeps them.
#'
#' @details
#' The sum reads the map's partials up to order three, which is where a
#' hand-written table begins to pay against the stencil: the numerical route's
#' accuracy fades with the order, and by the third it is several digits short of
#' the parent's own.
#'
#' @param distrib A reparametrized distribution, from [reparametrize()].
#' @param y The response, a numeric vector.
#' @param theta A named list of the new parameters, on the new parameter scale.
#' @param expected Should the expected derivatives be carried? A single logical,
#'   `FALSE` by default. Under expectation the term in the parent's score drops,
#'   the score having mean zero.
#' @param scale Either `"parameter"` (the default) or `"link"`, handled by the
#'   generic after this method has returned.
#' @param approx The strategy the parent uses for an expected derivative it has
#'   no closed form for, one of `"integrate"`, `"bartlett"`, `"mc"` or `"opg"`.
#'   Passed straight through; read only when `expected` is `TRUE` and only by a
#'   parent that approximates.
#' @param nsim The number of draws for `approx = "mc"`, 10000 by default.
#'   Passed through.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of numeric vectors keyed as
#'   [`deriv_names(distrib@params, 3)`][deriv_names], each the length of `y`
#'   recycled against `theta`.
#'
#' @seealso [chain_derivatives()] for the identity;
#'   [distrib_hessian.ReparamContinuousDistrib()] for the order below;
#'   [reparam_tables()], whose accuracy this inherits; [reparametrize()].
#'
#' @examples
#' d <- reparametrize(
#'   gaussian1_distrib(),
#'   map = function(psi) list(mu = psi$mu, sigma = sqrt(psi$sigma2)),
#'   params = c("mu", "sigma2"),
#'   bounds = list(mu = c(-Inf, Inf), sigma2 = c(0, Inf)),
#'   links = list(mu = linkfunctions7::identity_link(),
#'                sigma2 = linkfunctions7::log_link())
#' )
#' th <- list(mu = 1, sigma2 = 4)
#'
#' # Four components for a two-parameter family.
#' names(distrib_deriv3(d, c(0, 1, 2), th))
#'
#' @keywords internal
reparam_deriv3 <- function(distrib, y, theta, expected = FALSE,
                           scale = c("parameter", "link"),
                           approx = c("integrate", "bartlett", "mc", "opg"),
                           nsim = 10000, ...) {
  reparam_chain(distrib, y, theta, 3L, expected = expected)
}
S7::method(distrib_deriv3, ReparamContinuousDistrib) <- reparam_deriv3
S7::method(distrib_deriv3, ReparamDiscreteDistrib) <- reparam_deriv3

#' @title Fourth-Order Derivatives of a Reparametrized Distribution
#' @name distrib_deriv4.ReparamContinuousDistrib
#' @aliases distrib_deriv4.ReparamDiscreteDistrib
#'
#' @description
#' The partition sum at order four, observed or expected. This is the highest
#' order the toolkit carries, and it is what a marginal criterion's exact
#' gradient reads, so a reparametrized family used inside a penalized fit
#' reaches it through this method.
#'
#' @param distrib A reparametrized distribution, from [reparametrize()].
#' @param y The response, a numeric vector.
#' @param theta A named list of the new parameters, on the new parameter scale.
#' @param expected Should the expected derivatives be carried? A single logical,
#'   `FALSE` by default.
#' @param scale Either `"parameter"` (the default) or `"link"`, handled by the
#'   generic after this method has returned.
#' @param approx The strategy the parent uses for an expected derivative it has
#'   no closed form for, one of `"integrate"`, `"bartlett"`, `"mc"` or `"opg"`.
#'   Passed straight through.
#' @param nsim The number of draws for `approx = "mc"`, 10000 by default.
#'   Passed through.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of numeric vectors keyed as
#'   [`deriv_names(distrib@params, 4)`][deriv_names], each the length of `y`
#'   recycled against `theta`.
#'
#' @seealso [chain_derivatives()] for the identity;
#'   [distrib_deriv3.ReparamContinuousDistrib()] for the order below;
#'   [reparam_tables()], whose accuracy this inherits; [reparametrize()].
#'
#' @examples
#' d <- reparametrize(
#'   gaussian1_distrib(),
#'   map = function(psi) list(mu = psi$mu, sigma = sqrt(psi$sigma2)),
#'   params = c("mu", "sigma2"),
#'   bounds = list(mu = c(-Inf, Inf), sigma2 = c(0, Inf)),
#'   links = list(mu = linkfunctions7::identity_link(),
#'                sigma2 = linkfunctions7::log_link())
#' )
#' th <- list(mu = 1, sigma2 = 4)
#'
#' # Five components for a two-parameter family.
#' names(distrib_deriv4(d, c(0, 1, 2), th))
#'
#' @keywords internal
reparam_deriv4 <- function(distrib, y, theta, expected = FALSE,
                           scale = c("parameter", "link"),
                           approx = c("integrate", "bartlett", "mc", "opg"),
                           nsim = 10000, ...) {
  reparam_chain(distrib, y, theta, 4L, expected = expected)
}
S7::method(distrib_deriv4, ReparamContinuousDistrib) <- reparam_deriv4
S7::method(distrib_deriv4, ReparamDiscreteDistrib) <- reparam_deriv4

#' @title Response Derivative of a Reparametrized Distribution
#' @name distrib_grad_y.ReparamContinuousDistrib
#'
#' @description
#' The parent's, unchanged. A reparametrization touches the parameters and
#' leaves the response alone, so \eqn{\partial\log f/\partial y} is the
#' parent's at the mapped parameters and no chain rule enters.
#'
#' @details
#' Only the continuous class registers this. A discrete family has no response
#' derivative, and the base class refuses one, so there is nothing for a
#' discrete reparametrization to delegate.
#'
#' @param distrib A reparametrized distribution, from [reparametrize()].
#' @param y The response, a numeric vector.
#' @param theta A named list of the new parameters, on the new parameter scale.
#' @param ... Passed to the parent's method.
#'
#' @return A numeric vector, the length of `y` recycled against `theta`.
#'
#' @seealso [distrib_hess_y.ReparamContinuousDistrib()] for the second order;
#'   [distrib_gradient.ReparamContinuousDistrib()], where the map does enter;
#'   [reparametrize()].
#'
#' @examples
#' d <- reparametrize(
#'   gaussian1_distrib(),
#'   map = function(psi) list(mu = psi$mu, sigma = sqrt(psi$sigma2)),
#'   params = c("mu", "sigma2"),
#'   bounds = list(mu = c(-Inf, Inf), sigma2 = c(0, Inf)),
#'   links = list(mu = linkfunctions7::identity_link(),
#'                sigma2 = linkfunctions7::log_link())
#' )
#'
#' # Identical to the parent's at the mapped parameters.
#' all.equal(distrib_grad_y(d, c(0, 1, 2), list(mu = 1, sigma2 = 4)),
#'           distrib_grad_y(gaussian1_distrib(), c(0, 1, 2),
#'                          list(mu = 1, sigma = 2)))
#'
#' @keywords internal
S7::method(distrib_grad_y, ReparamContinuousDistrib) <- function(distrib, y, theta, ...) {
  distrib_grad_y(distrib@parent_distrib, y, reparam_theta(distrib, theta), ...)
}

#' @title Second Response Derivative of a Reparametrized Distribution
#' @name distrib_hess_y.ReparamContinuousDistrib
#'
#' @description
#' The parent's, unchanged, for the same reason as the first: a
#' reparametrization moves the parameters and leaves the response where it was.
#' Only the continuous class registers it.
#'
#' @param distrib A reparametrized distribution, from [reparametrize()].
#' @param y The response, a numeric vector.
#' @param theta A named list of the new parameters, on the new parameter scale.
#' @param ... Passed to the parent's method.
#'
#' @return A numeric vector, the length of `y` recycled against `theta`.
#'
#' @seealso [distrib_grad_y.ReparamContinuousDistrib()] for the first order;
#'   [reparametrize()].
#'
#' @examples
#' d <- reparametrize(
#'   gaussian1_distrib(),
#'   map = function(psi) list(mu = psi$mu, sigma = sqrt(psi$sigma2)),
#'   params = c("mu", "sigma2"),
#'   bounds = list(mu = c(-Inf, Inf), sigma2 = c(0, Inf)),
#'   links = list(mu = linkfunctions7::identity_link(),
#'                sigma2 = linkfunctions7::log_link())
#' )
#'
#' # A Gaussian's is -1/sigma^2, whatever the parametrization.
#' distrib_hess_y(d, c(0, 1, 2), list(mu = 1, sigma2 = 4))
#'
#' @keywords internal
S7::method(distrib_hess_y, ReparamContinuousDistrib) <- function(distrib, y, theta, ...) {
  distrib_hess_y(distrib@parent_distrib, y, reparam_theta(distrib, theta), ...)
}

#' @title Atoms of a Reparametrized Distribution
#' @name distrib_atoms.ReparamContinuousDistrib
#' @aliases distrib_atoms.ReparamDiscreteDistrib
#'
#' @description
#' The parent's, unchanged: a reparametrization does not move mass. Delegating
#' matters because [check_distrib()] consults this to decide which of its
#' checks apply, and a reparametrized mixed family would otherwise be validated
#' as though it were purely continuous.
#'
#' @param distrib A reparametrized distribution, from [reparametrize()].
#' @param theta A named list of the new parameters, on the new parameter scale.
#'
#' @return A named list with components `y` (the atom locations) and `p`
#'   (their probabilities), both numeric vectors of the same length, possibly
#'   of length zero.
#'
#' @seealso [distrib_atoms()] for the generic and what consults it;
#'   [check_distrib()]; [reparametrize()].
#'
#' @examples
#' # A reparametrized zero-adjusted gamma keeps its atom at zero.
#' za <- zero_adjusted(gamma2_distrib())
#' d <- reparametrize(
#'   za,
#'   map = function(psi) list(mu = psi$mu, sigma2 = psi$s, za = psi$p0),
#'   params = c("mu", "s", "p0"),
#'   bounds = list(mu = c(0, Inf), s = c(0, Inf), p0 = c(0, 1)),
#'   links = list(mu = linkfunctions7::log_link(),
#'                s = linkfunctions7::log_link(),
#'                p0 = linkfunctions7::logit_link())
#' )
#' distrib_atoms(d, list(mu = 2, s = 1, p0 = 0.3))
#'
#' @keywords internal
reparam_atoms <- function(distrib, theta) {
  distrib_atoms(distrib@parent_distrib, reparam_theta(distrib, theta))
}
S7::method(distrib_atoms, ReparamContinuousDistrib) <- reparam_atoms
S7::method(distrib_atoms, ReparamDiscreteDistrib) <- reparam_atoms


#' Describe a Trial Parameter Value
#'
#' @description
#' Renders the point a map was probed at as `(a = 1, b = 2)`, so that a failure
#' in [reparametrize()]'s construction check names the point instead of leaving
#' the caller to guess it.
#'
#' @param probe A named list of parameter values, each of length 1.
#'
#' @return A single string.
#'
#' @seealso [reparametrize()], the one caller.
#'
#' @keywords internal
.describe_probe <- function(probe) {
  paste0("(", paste(names(probe), unlist(probe), sep = " = ",
                    collapse = ", "), ")")
}

#' Write a Distribution in Different Coordinates
#'
#' @description
#' Returns the same law as `distrib`, parametrized by quantities of the
#' caller's choosing.
#'
#' @details
#' A reparametrization is not a link. A link changes the scale a parameter is
#' *modeled* on and leaves the parameter what it was; here the parameter
#' *is* the new quantity, and that is what the estimate, the standard
#' error and the confidence interval describe.
#'
#' **The map is written in ordinary R.** It takes a named list of the new
#' parameters and returns a named list of the parent's, and nothing in it
#' mentions derivatives:
#'
#' \preformatted{
#' function(psi) list(mu = psi$mean / gamma(1 + 1 / psi$sigma),
#'                    sigma = psi$sigma)
#' }
#'
#' **Where the map's derivatives come from.** Two routes, and `map_derivs`
#' chooses between them. Supply it, and the tables are hand-written closed
#' forms, exact at every order with no chain rule transcribed; the shipped
#' second parametrizations all do, and their formulas live in
#' [reparam_map_derivs()]. Leave it `NULL`, and each partial the chain needs
#' comes from one central stencil of [numericals7::fd_derivative()] applied to
#' the analytic map, a single stencil per order and never a chain of
#' differences. The stencil route is good to about \eqn{10^{-8}} at first
#' order and fades with the order, so a family fitted in earnest is worth a
#' table.
#'
#' **What is exact and what is inherited.** The derivatives of the
#' log-density are carried by
#' \deqn{\ell^{(I)}(\psi) = \sum_{\pi} \sum_{i_1 \dots i_{|\pi|}}
#'       \ell^{(i_1 \dots i_{|\pi|})}(\theta)
#'       \prod_{B \in \pi} \frac{\partial^{|B|}\theta_{i_B}}{\partial \psi_B},}
#' the sum over set partitions, so every order the parent has in closed form
#' survives in closed form. Expectation being linear and the map
#' deterministic, the expected derivatives obey the same formula, and the term
#' carrying the second derivative of the map drops because the score has mean
#' zero: a parent with an exact expected information gives an exact one here.
#' The density, the distribution function, the quantile function, the generator
#' and the response derivatives are the parent's at the mapped parameters.
#'
#' **Cost.** A parameter that varies by observation needs one run of the
#' map per observation; scalar parameters need one run in total, which is the
#' case a fit is in.
#'
#' @param distrib The distribution to rewrite.
#' @param map A function of a named list of the new parameters, returning a
#'   named list of the parent's.
#' @param params A character vector naming the new parameters.
#' @param bounds A named list of length-two numeric vectors, the open interval
#'   each new parameter lives in.
#' @param links A named list of \pkg{linkfunctions7} links, one per new
#'   parameter.
#' @param map_derivs An optional function returning, for each parent
#'   parameter, the non-zero partial derivatives of the map with respect to
#'   the new parameters to fourth order, keyed by the sorted tuple of
#'   new-parameter positions ("1", "1,2", "2,2,3,3", ...); a missing key is
#'   an exact zero. The shipped second parametrizations supply hand-written
#'   tables (see [reparam_map_derivs()]); when `NULL`, each
#'   needed partial comes from one finite-difference stencil on the map.
#' @param interpretation An optional named character vector describing each new
#'   parameter; defaults to the parameter names.
#' @param name An optional name for the result; defaults to the parent's with
#'   the new parameters appended.
#'
#' @return A distribution object of class
#'   [ReparamContinuousDistrib()] or
#'   [ReparamDiscreteDistrib()].
#'
#' @seealso [fixed()], [transformation()]
#'
#' @examples
#' # a gaussian in its mean and variance, obtained rather than written
#' d <- reparametrize(
#'   gaussian1_distrib(),
#'   map = function(psi) list(mu = psi$mu, sigma = sqrt(psi$sigma2)),
#'   params = c("mu", "sigma2"),
#'   bounds = list(mu = c(-Inf, Inf), sigma2 = c(0, Inf)),
#'   links = list(mu = linkfunctions7::identity_link(),
#'                sigma2 = linkfunctions7::log_link())
#' )
#' theta <- list(mu = 1, sigma2 = 4)
#' distrib_pdf(d, c(0, 1, 2), theta)
#'
#' # and it agrees with the family written out by hand
#' distrib_pdf(gaussian2_distrib(), c(0, 1, 2), theta)
#'
#' @export
reparametrize <- function(distrib, map, params, bounds, links,
                          map_derivs = NULL, interpretation = NULL,
                          name = NULL) {
  if (!S7::S7_inherits(distrib, continuous_distrib) &&
    !S7::S7_inherits(distrib, discrete_distrib)) {
    stop("'distrib' must inherit from 'continuous_distrib' or 'discrete_distrib'.",
      call. = FALSE
    )
  }
  if (!is.function(map)) {
    stop("'map' must be a function of the new parameters.", call. = FALSE)
  }
  if (!is.character(params) || !length(params) || anyDuplicated(params)) {
    stop("'params' must be distinct names for the new parameters.", call. = FALSE)
  }
  for (arg in c("bounds", "links")) {
    v <- get(arg)
    if (!is.list(v) || !setequal(names(v), params)) {
      stop(sprintf("'%s' must be a named list with one entry per new parameter.",
                   arg), call. = FALSE)
    }
  }
  if (is.null(interpretation)) {
    interpretation <- stats::setNames(params, params)
  }

  # The map has to return the parent's parameters, and the cheapest place to
  # find out is here rather than several frames down inside a derivative.
  # A trial value INSIDE each parameter's own interval. Taking the midpoint of
  # a half-line gives infinity and falling back to zero puts the probe outside
  # the domain, which for a map like sqrt((nu-2)/nu) produces a NaN whose
  # warning names neither the parameter nor the reason.
  probe <- stats::setNames(lapply(params, function(k) {
    b <- bounds[[k]]
    if (is.finite(b[1L]) && is.finite(b[2L])) return(mean(b))
    if (is.finite(b[1L])) return(b[1L] + 1)
    if (is.finite(b[2L])) return(b[2L] - 1)
    0
  }), params)
  got <- withCallingHandlers(
    tryCatch(map(probe), error = function(e) {
      stop(sprintf("'map' failed at %s: %s", .describe_probe(probe),
                   conditionMessage(e)), call. = FALSE)
    }),
    warning = function(w) {
      stop(sprintf("'map' warned at %s: %s", .describe_probe(probe),
                   conditionMessage(w)), call. = FALSE)
    }
  )
  if (any(!vapply(got, function(v) all(is.finite(unlist(v))), logical(1)))) {
    stop(sprintf("'map' returned a non-finite value at %s.",
                 .describe_probe(probe)), call. = FALSE)
  }
  if (!is.list(got) || is.null(names(got))) {
    stop("'map' must return a NAMED list of the parent's parameters.",
      call. = FALSE
    )
  }
  missing <- setdiff(distrib@params, names(got))
  if (length(missing)) {
    stop(sprintf(
      "'map' did not return '%s'. It must return every parameter of '%s': %s.",
      missing[1L], distrib@distrib_name,
      paste(distrib@params, collapse = ", ")
    ), call. = FALSE)
  }

  if (is.null(name)) {
    name <- sprintf("%s [%s]", distrib@distrib_name,
                    paste(params, collapse = ", "))
  }

  cls <- if (S7::S7_inherits(distrib, discrete_distrib)) {
    ReparamDiscreteDistrib
  } else {
    ReparamContinuousDistrib
  }

  cls(
    distrib_name = name,
    dimension = "univariate",
    bounds = distrib@bounds,
    params = params,
    params_interpretation = interpretation[params],
    n_params = length(params),
    params_bounds = bounds[params],
    link_params = links[params],
    params_smooth = stats::setNames(rep(TRUE, length(params)), params),
    parent_distrib = distrib,
    reparam_map = map,
    reparam_derivs = if (is.null(map_derivs)) {
      reparam_stencil_derivs(map, params, distrib@params)
    } else {
      map_derivs
    }
  )
}


#' @title Moments of a Reparametrized Distribution
#' @name mean.ReparamContinuousDistrib
#' @aliases mean.ReparamDiscreteDistrib variance.ReparamDiscreteDistrib
#'   skewness.ReparamDiscreteDistrib kurtosis.ReparamDiscreteDistrib
#'
#' @description
#' `mean()`, [variance()], [skewness()] and [kurtosis()] all delegate to the
#' parent at the mapped parameters. A reparametrization does not change the
#' law, so it does not change a moment; what changes is the coordinates the
#' moment is computed from. Delegating keeps the parent's closed forms, where
#' falling through to the `distrib` default would run a quadrature.
#'
#' @param x A reparametrized distribution, from [reparametrize()].
#' @param theta A named list of the new parameters, on the new parameter scale.
#' @param ... Passed to the parent's method, and from there to [moment()] if
#'   the parent has no closed form.
#'
#' @return A numeric vector: means for `mean()`, variances for [variance()],
#'   and the standardized third and fourth moments for [skewness()] and
#'   [kurtosis()], the excess in the last case. `NaN` or `Inf` wherever the
#'   parent's own method gives one.
#'
#' @seealso [moment()] for the numerical route the parent may take;
#'   [mean.distrib()] for the default this replaces; [reparametrize()].
#'
#' @examples
#' d <- reparametrize(
#'   gaussian1_distrib(),
#'   map = function(psi) list(mu = psi$mu, sigma = sqrt(psi$sigma2)),
#'   params = c("mu", "sigma2"),
#'   bounds = list(mu = c(-Inf, Inf), sigma2 = c(0, Inf)),
#'   links = list(mu = linkfunctions7::identity_link(),
#'                sigma2 = linkfunctions7::log_link())
#' )
#' th <- list(mu = 1, sigma2 = 4)
#'
#' # The parent's closed forms, in the new coordinates.
#' c(mean(d, th), variance(d, th), skewness(d, th), kurtosis(d, th))
#'
#' @keywords internal
reparam_mean <- function(x, theta, ...) {
  mean(x@parent_distrib, reparam_theta(x, theta), ...)
}
S7::method(mean, ReparamContinuousDistrib) <- reparam_mean
S7::method(mean, ReparamDiscreteDistrib) <- reparam_mean

#' @rdname mean.ReparamContinuousDistrib
#' @name variance.ReparamContinuousDistrib
#' @keywords internal
reparam_variance <- function(x, theta, ...) {
  variance(x@parent_distrib, reparam_theta(x, theta), ...)
}
S7::method(variance, ReparamContinuousDistrib) <- reparam_variance
S7::method(variance, ReparamDiscreteDistrib) <- reparam_variance

#' @rdname mean.ReparamContinuousDistrib
#' @name skewness.ReparamContinuousDistrib
#' @keywords internal
reparam_skewness <- function(x, theta, ...) {
  skewness(x@parent_distrib, reparam_theta(x, theta), ...)
}
S7::method(skewness, ReparamContinuousDistrib) <- reparam_skewness
S7::method(skewness, ReparamDiscreteDistrib) <- reparam_skewness

#' @rdname mean.ReparamContinuousDistrib
#' @name kurtosis.ReparamContinuousDistrib
#' @keywords internal
reparam_kurtosis <- function(x, theta, ...) {
  kurtosis(x@parent_distrib, reparam_theta(x, theta), ...)
}
S7::method(kurtosis, ReparamContinuousDistrib) <- reparam_kurtosis
S7::method(kurtosis, ReparamDiscreteDistrib) <- reparam_kurtosis
