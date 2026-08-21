#' @include distrib.R generics.R moments.R
NULL

#' @title S7 Classes for a Reparametrized Distribution
#' @name ReparamContinuousDistrib
#'
#' @description
#' The wrapper \code{\link{reparametrize}} returns: the same law as its parent,
#' written in different coordinates. There is one class per kind of parent, so
#' that a continuous parent keeps the defaults registered on
#' \code{\link{continuous_distrib}} and a discrete one those of
#' \code{\link{discrete_distrib}}.
#'
#' @inheritParams distrib
#' @param parent_distrib The distribution being rewritten.
#' @param reparam_map The map from the new parameters to the parent's.
#' @param reparam_derivs The function returning the map's keyed partial
#'   tables, as \code{\link{reparam_tables}} consumes them.
#' @return An object of the corresponding class.
#' @seealso \code{\link{reparametrize}}
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
#' \code{TRUE} for either of the two wrapper classes.
#'
#' @param distrib An object inheriting from class \code{"distrib"}.
#'
#' @return A single logical.
#'
#' @seealso \code{\link{reparametrize}}
#'
#' @keywords internal
is_reparam <- function(distrib) {
  S7::S7_inherits(distrib, ReparamContinuousDistrib) ||
    S7::S7_inherits(distrib, ReparamDiscreteDistrib)
}

#' The Parent's Parameters at the New Ones
#'
#' @description
#' Runs the map on plain numbers, which is what every probability function
#' needs before delegating to the parent.
#'
#' @param distrib A reparametrized distribution.
#' @param theta A named list of the new parameters.
#'
#' @return A named list of the parent's parameters.
#'
#' @seealso \code{\link{reparametrize}}
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
#' Returns, for every parent parameter, the partial derivatives of the map
#' component with respect to the new parameters, keyed by the sorted tuple
#' of new-parameter positions ("1", "1,2", "2,2,3,3", ...). A missing key
#' is an exact zero.
#'
#' @details
#' When the family supplies \code{map_derivs}, the tables are its
#' hand-written closed forms; the shipped second parametrizations all do,
#' and the formulas live in \code{\link{reparam_map_derivs}}. Otherwise
#' each needed partial comes from one finite-difference stencil of
#' \code{\link[numericals7]{fd_derivative}} applied to the analytic map --
#' a single stencil per order, never a chain of differences, at the
#' accuracy that construction carries (about 1e-8 at first order, fading
#' with the order). Exact tables are therefore the recommendation for any
#' family fitted in earnest.
#'
#' @param distrib A reparametrized distribution.
#' @param theta A named list of the new parameters, already aligned.
#'
#' @return A list over parent parameters of keyed partial tables.
#'
#' @seealso \code{\link{reparametrize}}
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
#' The numerical fallback behind \code{\link{reparam_tables}}: every
#' partial of the map is one central stencil of
#' \code{\link[numericals7]{fd_derivative}} in each direction, applied to
#' the analytic map, iterated across DISTINCT directions (a cross-variable
#' composition, which is not the forbidden same-variable nesting).
#'
#' @param map The map function.
#' @param params The new parameter names.
#' @param parent_params The parent parameter names.
#'
#' @return A function usable as \code{reparam_derivs}.
#'
#' @seealso \code{\link{reparametrize}}
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

#' The Chain Rule of Any Order Through a Reparametrization
#'
#' @description
#' Carries a derivative of the parent's log-density into the new coordinates.
#'
#' @details
#' With \eqn{\theta = h(\psi)}, the derivative of order \eqn{|I|} is
#' \deqn{\ell^{(I)}(\psi) = \sum_{\pi} \sum_{i_1 \dots i_{|\pi|}}
#'       \ell^{(i_1 \dots i_{|\pi|})}(\theta)
#'       \prod_{B \in \pi} \frac{\partial^{|B|}\theta_{i_B}}{\partial \psi_B}}
#' the outer sum running over the set partitions of the \strong{positions} of
#' \eqn{I} and the inner one over the assignment of a parent parameter to each
#' block. This is Faa di Bruno with a dense Jacobian, and it is the same
#' partition enumeration the wrappers of \code{\link{zero_inflated}} and the
#' rest already use.
#'
#' Blocks index positions rather than variables, which is what makes a repeated
#' index carry its multiplicity without a factorial correction.
#'
#' Expectation is linear and \eqn{h} is deterministic, so the expected
#' derivatives obey the same formula with the parent's expected derivatives in
#' place of the observed ones. The first-order term drops, the score having
#' mean zero.
#'
#' @param distrib A reparametrized distribution.
#' @param y The response.
#' @param theta A named list of the new parameters.
#' @param order The derivative order, 1 to 4.
#' @param expected Logical; if \code{TRUE}, carries the expected derivatives.
#'
#' @return A named list of component vectors, keyed as
#'   \code{\link{deriv_names}} keys them.
#'
#' @seealso \code{\link{reparametrize}}
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
#' Carries the parent's derivatives into new coordinates, given the jets of the
#' map. Separated from \code{\link{reparam_chain}} so that a family written in
#' its own right, rather than obtained through \code{\link{reparametrize}},
#' can use the same machinery instead of a second copy of it.
#'
#' @param parent The distribution whose derivatives are being carried.
#' @param y The response.
#' @param th_par The parent's parameters, as plain numbers.
#' @param maps The keyed partial tables of the map, as
#'   \code{\link{reparam_tables}} returns them.
#' @param new_params The names of the new parameters.
#' @param order The derivative order, 1 to 4.
#' @param expected Logical; if \code{TRUE}, carries the expected derivatives.
#'
#' @return A named list of component vectors.
#'
#' @seealso \code{\link{reparametrize}}
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
#' The partition sum of \code{\link{chain_derivatives}}, taking the inner
#' derivatives as tables rather than fetching them from a distribution.
#'
#' @details
#' Separated so that a family whose own derivatives are easiest to write in
#' coordinates it does not expose can carry them into the ones it does,
#' without a second copy of the enumeration. \code{\link{enet_distrib}} is
#' the case: its derivatives are compact in the two rates \eqn{(a, c)} and
#' the family is parametrized by \eqn{(\lambda, \alpha)}, which is a
#' bilinear map away.
#'
#' @param D A list of length \code{order}; \code{D[[k]]} is the inner
#'   derivative table of order \code{k}, keyed as
#'   \code{\link{deriv_names}(inner_params, k)}.
#' @param inner_params The names of the inner coordinates.
#' @param maps Per inner coordinate, a keyed table of the map's partials in
#'   the outer coordinates, as \code{\link{reparam_tables}} returns them.
#' @param new_params The names of the outer parameters.
#' @param order The derivative order, 1 to 4.
#' @param n The length to recycle a constant component to.
#'
#' @return A named list of component vectors.
#'
#' @seealso \code{\link{chain_derivatives}}
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
#' @description The parent's density at the mapped parameters.
#' @param distrib A reparametrized distribution.
#' @param y A numeric vector of observations.
#' @param theta A named list of the new parameters.
#' @param log Logical; if \code{TRUE}, returns the log-density.
#' @return A numeric vector.
#' @seealso \code{\link{reparametrize}}
#' @keywords internal
reparam_pdf <- function(distrib, y, theta, log = FALSE, ...) {
  distrib_pdf(distrib@parent_distrib, y, reparam_theta(distrib, theta), log = log)
}
S7::method(distrib_pdf, ReparamContinuousDistrib) <- reparam_pdf
S7::method(distrib_pdf, ReparamDiscreteDistrib) <- reparam_pdf

#' @title Distribution Function of a Reparametrized Distribution
#' @name distrib_cdf.ReparamContinuousDistrib
#' @description The parent's distribution function at the mapped parameters.
#' @param distrib A reparametrized distribution.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of the new parameters.
#' @param lower.tail Logical; if \code{TRUE}, probabilities are \eqn{P(Y \le q)}.
#' @param log.p Logical; if \code{TRUE}, returns log-probabilities.
#' @param ... Passed to the parent.
#' @return A numeric vector.
#' @seealso \code{\link{reparametrize}}
#' @keywords internal
reparam_cdf <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE, ...) {
  distrib_cdf(distrib@parent_distrib, q, reparam_theta(distrib, theta),
              lower.tail = lower.tail, log.p = log.p, ...)
}
S7::method(distrib_cdf, ReparamContinuousDistrib) <- reparam_cdf
S7::method(distrib_cdf, ReparamDiscreteDistrib) <- reparam_cdf

#' @title Quantile Function of a Reparametrized Distribution
#' @name distrib_quantile.ReparamContinuousDistrib
#' @description The parent's quantile function at the mapped parameters.
#' @param distrib A reparametrized distribution.
#' @param p A numeric vector of probabilities.
#' @param theta A named list of the new parameters.
#' @param lower.tail Logical; if \code{TRUE}, probabilities are \eqn{P(Y \le q)}.
#' @param log.p Logical; if \code{TRUE}, \code{p} is a log-probability.
#' @param ... Passed to the parent.
#' @return A numeric vector.
#' @seealso \code{\link{reparametrize}}
#' @keywords internal
reparam_quantile <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE, ...) {
  distrib_quantile(distrib@parent_distrib, p, reparam_theta(distrib, theta),
                   lower.tail = lower.tail, log.p = log.p, ...)
}
S7::method(distrib_quantile, ReparamContinuousDistrib) <- reparam_quantile
S7::method(distrib_quantile, ReparamDiscreteDistrib) <- reparam_quantile

#' @title Random Generation from a Reparametrized Distribution
#' @name distrib_rng.ReparamContinuousDistrib
#' @description The parent's generator at the mapped parameters.
#' @param distrib A reparametrized distribution.
#' @param n The number of draws.
#' @param theta A named list of the new parameters.
#' @return A numeric vector.
#' @seealso \code{\link{reparametrize}}
#' @keywords internal
reparam_rng <- function(distrib, n, theta) {
  distrib_rng(distrib@parent_distrib, n, reparam_theta(distrib, theta))
}
S7::method(distrib_rng, ReparamContinuousDistrib) <- reparam_rng
S7::method(distrib_rng, ReparamDiscreteDistrib) <- reparam_rng

#' @title Gradient of a Reparametrized Distribution
#' @name distrib_gradient.ReparamContinuousDistrib
#' @description
#' The parent's score carried by the Jacobian of the map,
#' \eqn{\ell^{(a)} = \sum_i \ell^{(i)} \partial\theta_i/\partial\psi_a}.
#' @param distrib A reparametrized distribution.
#' @param y The response.
#' @param theta A named list of the new parameters.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param ... Unused.
#' @return A named list of first derivatives.
#' @seealso \code{\link{reparametrize}}
#' @keywords internal
reparam_gradient <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  reparam_chain(distrib, y, theta, 1L)
}
S7::method(distrib_gradient, ReparamContinuousDistrib) <- reparam_gradient
S7::method(distrib_gradient, ReparamDiscreteDistrib) <- reparam_gradient

#' @title Observed Hessian of a Reparametrized Distribution
#' @name distrib_hessian.ReparamContinuousDistrib
#' @description
#' The second-order chain rule, which keeps the term in the parent's score and
#' the second derivative of the map:
#' \eqn{\ell^{(ab)} = \sum_{ij}\ell^{(ij)}h^i_a h^j_b + \sum_i \ell^{(i)} h^i_{ab}}.
#' @param distrib A reparametrized distribution.
#' @param y The response.
#' @param theta A named list of the new parameters.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param ... Unused.
#' @return A named list of second derivatives.
#' @seealso \code{\link{reparametrize}}
#' @keywords internal
reparam_hessian <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  reparam_chain(distrib, y, theta, 2L)[hess_names(distrib@params)]
}
S7::method(distrib_hessian, ReparamContinuousDistrib) <- reparam_hessian
S7::method(distrib_hessian, ReparamDiscreteDistrib) <- reparam_hessian

#' @title Expected Hessian of a Reparametrized Distribution
#' @name distrib_expected_hessian.ReparamContinuousDistrib
#' @description
#' The same congruence applied to the parent's expected information,
#' \eqn{\mathbb{E}[\ell^{(ab)}] = \sum_{ij}\mathbb{E}[\ell^{(ij)}]h^i_a h^j_b}.
#' The second-derivative term of the map drops out because the score has mean
#' zero, so a parent with an exact expected information gives an exact one
#' here.
#' @param distrib A reparametrized distribution.
#' @param y The response.
#' @param theta A named list of the new parameters.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param approx Passed to the parent when it has no closed form.
#' @param nsim Passed to the parent.
#' @param ... Unused.
#' @return A named list of expected second derivatives.
#' @seealso \code{\link{reparametrize}}
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
#' @description The partition sum at order three, observed or expected.
#' @param distrib A reparametrized distribution.
#' @param y The response.
#' @param theta A named list of the new parameters.
#' @param expected Logical; if \code{TRUE}, carries the expected derivatives.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param approx Passed to the parent.
#' @param nsim Passed to the parent.
#' @param ... Unused.
#' @return A named list of third-derivative components.
#' @seealso \code{\link{reparametrize}}
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
#' @description The partition sum at order four, observed or expected.
#' @param distrib A reparametrized distribution.
#' @param y The response.
#' @param theta A named list of the new parameters.
#' @param expected Logical; if \code{TRUE}, carries the expected derivatives.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param approx Passed to the parent.
#' @param nsim Passed to the parent.
#' @param ... Unused.
#' @return A named list of fourth-derivative components.
#' @seealso \code{\link{reparametrize}}
#' @keywords internal
reparam_deriv4 <- function(distrib, y, theta, expected = FALSE,
                           scale = c("parameter", "link"),
                           approx = c("integrate", "bartlett", "mc", "opg"),
                           nsim = 10000, ...) {
  reparam_chain(distrib, y, theta, 4L, expected = expected)
}
S7::method(distrib_deriv4, ReparamContinuousDistrib) <- reparam_deriv4
S7::method(distrib_deriv4, ReparamDiscreteDistrib) <- reparam_deriv4

#' @title Response Derivatives of a Reparametrized Distribution
#' @name distrib_grad_y.ReparamContinuousDistrib
#' @description
#' The parent's, unchanged: a reparametrization touches the parameters and not
#' the response.
#' @param distrib A reparametrized distribution.
#' @param y The response.
#' @param theta A named list of the new parameters.
#' @param ... Passed to the parent.
#' @return A numeric vector.
#' @seealso \code{\link{reparametrize}}
#' @keywords internal
S7::method(distrib_grad_y, ReparamContinuousDistrib) <- function(distrib, y, theta, ...) {
  distrib_grad_y(distrib@parent_distrib, y, reparam_theta(distrib, theta), ...)
}

#' @title Second Response Derivative of a Reparametrized Distribution
#' @name distrib_hess_y.ReparamContinuousDistrib
#' @description The parent's, unchanged.
#' @param distrib A reparametrized distribution.
#' @param y The response.
#' @param theta A named list of the new parameters.
#' @param ... Passed to the parent.
#' @return A numeric vector.
#' @seealso \code{\link{reparametrize}}
#' @keywords internal
S7::method(distrib_hess_y, ReparamContinuousDistrib) <- function(distrib, y, theta, ...) {
  distrib_hess_y(distrib@parent_distrib, y, reparam_theta(distrib, theta), ...)
}

#' @title Atoms of a Reparametrized Distribution
#' @name distrib_atoms.ReparamContinuousDistrib
#' @description The parent's, unchanged: a reparametrization does not move mass.
#' @param distrib A reparametrized distribution.
#' @param theta A named list of the new parameters.
#' @return A numeric vector of atom locations, possibly empty.
#' @seealso \code{\link{reparametrize}}
#' @keywords internal
reparam_atoms <- function(distrib, theta) {
  distrib_atoms(distrib@parent_distrib, reparam_theta(distrib, theta))
}
S7::method(distrib_atoms, ReparamContinuousDistrib) <- reparam_atoms
S7::method(distrib_atoms, ReparamDiscreteDistrib) <- reparam_atoms


#' Describe a Trial Parameter Value
#'
#' @description
#' Renders the point a map was probed at, so that a failure names it instead of
#' leaving the caller to guess.
#'
#' @param probe A named list of parameter values.
#'
#' @return A single string.
#'
#' @seealso \code{\link{reparametrize}}
#'
#' @keywords internal
.describe_probe <- function(probe) {
  paste0("(", paste(names(probe), unlist(probe), sep = " = ",
                    collapse = ", "), ")")
}

#' Write a Distribution in Different Coordinates
#'
#' @description
#' Returns the same law as \code{distrib}, parametrized by quantities of the
#' caller's choosing.
#'
#' @details
#' A reparametrization is not a link. A link changes the scale a parameter is
#' \emph{modeled} on and leaves the parameter what it was; here the parameter
#' \emph{is} the new quantity, and that is what the estimate, the standard
#' error and the confidence interval describe.
#'
#' \strong{The map is written in ordinary R.} It takes a named list of the new
#' parameters and returns a named list of the parent's, and nothing in it
#' mentions derivatives:
#'
#' \preformatted{
#' function(psi) list(mu = psi$mean / gamma(1 + 1 / psi$sigma),
#'                    sigma = psi$sigma)
#' }
#'
#' The derivatives of the map come from running that same expression on
#' \strong{jets} -- values carrying every partial derivative to fourth order --
#' so they are exact at every order with no chain rule transcribed. The
#' arithmetic operators and the mathematical functions dispatch on them. A map
#' that branches on the value of a parameter is rejected rather than
#' approximated, a comparison having no derivative to carry.
#'
#' \strong{What is exact and what is inherited.} The derivatives of the
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
#' \strong{Cost.} A parameter that varies by observation needs one run of the
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
#'   tables (see \code{\link{reparam_map_derivs}}); when \code{NULL}, each
#'   needed partial comes from one finite-difference stencil on the map.
#' @param interpretation An optional named character vector describing each new
#'   parameter; defaults to the parameter names.
#' @param name An optional name for the result; defaults to the parent's with
#'   the new parameters appended.
#'
#' @return A distribution object of class
#'   \code{\link{ReparamContinuousDistrib}} or
#'   \code{\link{ReparamDiscreteDistrib}}.
#'
#' @seealso \code{\link{fixed}}, \code{\link{transformation}}
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
#' @description
#' The parent's, at the mapped parameters. A reparametrization does not change
#' the law, so it does not change a moment: what changes is the coordinates the
#' moment is computed from. Delegating keeps the parent's closed forms instead
#' of falling through to a quadrature.
#' @param x A reparametrized distribution.
#' @param theta A named list of the new parameters.
#' @param ... Passed to the parent.
#' @return A numeric vector.
#' @seealso \code{\link{reparametrize}}
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
