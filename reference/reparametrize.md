# Write a Distribution in Different Coordinates

Returns the same law as `distrib`, parametrized by quantities of the
caller's choosing.

## Usage

``` r
reparametrize(
  distrib,
  map,
  params,
  bounds,
  links,
  map_derivs = NULL,
  interpretation = NULL,
  name = NULL
)
```

## Arguments

- distrib:

  The distribution to rewrite.

- map:

  A function of a named list of the new parameters, returning a named
  list of the parent's.

- params:

  A character vector naming the new parameters.

- bounds:

  A named list of length-two numeric vectors, the open interval each new
  parameter lives in.

- links:

  A named list of linkfunctions7 links, one per new parameter.

- map_derivs:

  An optional function returning, for each parent parameter, the
  non-zero partial derivatives of the map with respect to the new
  parameters to fourth order, keyed by the sorted tuple of new-parameter
  positions ("1", "1,2", "2,2,3,3", ...); a missing key is an exact
  zero. The shipped second parametrizations supply hand-written tables
  (see
  [`reparam_map_derivs`](https://statmodels7.github.io/distributions7/reference/reparam_map_derivs.md));
  when `NULL`, each needed partial comes from one finite-difference
  stencil on the map.

- interpretation:

  An optional named character vector describing each new parameter;
  defaults to the parameter names.

- name:

  An optional name for the result; defaults to the parent's with the new
  parameters appended.

## Value

A distribution object of class
[`ReparamContinuousDistrib`](https://statmodels7.github.io/distributions7/reference/ReparamContinuousDistrib.md)
or
[`ReparamDiscreteDistrib`](https://statmodels7.github.io/distributions7/reference/ReparamContinuousDistrib.md).

## Details

A reparametrization is not a link. A link changes the scale a parameter
is *modeled* on and leaves the parameter what it was; here the parameter
*is* the new quantity, and that is what the estimate, the standard error
and the confidence interval describe.

**The map is written in ordinary R.** It takes a named list of the new
parameters and returns a named list of the parent's, and nothing in it
mentions derivatives:


    function(psi) list(mu = psi$mean / gamma(1 + 1 / psi$sigma),
                       sigma = psi$sigma)

The derivatives of the map come from running that same expression on
**jets** – values carrying every partial derivative to fourth order – so
they are exact at every order with no chain rule transcribed. The
arithmetic operators and the mathematical functions dispatch on them. A
map that branches on the value of a parameter is refused rather than
approximated, a comparison having no derivative to carry.

**What is exact and what is inherited.** The derivatives of the
log-density are carried by \$\$\ell^{(I)}(\psi) = \sum\_{\pi} \sum\_{i_1
\dots i\_{\|\pi\|}} \ell^{(i_1 \dots i\_{\|\pi\|})}(\theta) \prod\_{B
\in \pi} \frac{\partial^{\|B\|}\theta\_{i_B}}{\partial \psi_B},\$\$ the
sum over set partitions, so every order the parent has in closed form
survives in closed form. Expectation being linear and the map
deterministic, the expected derivatives obey the same formula, and the
term carrying the second derivative of the map drops because the score
has mean zero: a parent with an exact expected information gives an
exact one here. The density, the distribution function, the quantile
function, the generator and the response derivatives are the parent's at
the mapped parameters.

**Cost.** A parameter that varies by observation needs one run of the
map per observation; scalar parameters need one run in total, which is
the case a fit is in.

## See also

[`fixed`](https://statmodels7.github.io/distributions7/reference/fixed.md),
[`transformation`](https://statmodels7.github.io/distributions7/reference/transformation.md)

## Examples

``` r
# a gaussian in its mean and variance, obtained rather than written
d <- reparametrize(
  gaussian1_distrib(),
  map = function(psi) list(mu = psi$mu, sigma = sqrt(psi$sigma2)),
  params = c("mu", "sigma2"),
  bounds = list(mu = c(-Inf, Inf), sigma2 = c(0, Inf)),
  links = list(mu = linkfunctions7::identity_link(),
               sigma2 = linkfunctions7::log_link())
)
theta <- list(mu = 1, sigma2 = 4)
distrib_pdf(d, c(0, 1, 2), theta)
#> [1] 0.1760327 0.1994711 0.1760327

# and it agrees with the family written out by hand
distrib_pdf(gaussian2_distrib(), c(0, 1, 2), theta)
#> [1] 0.1760327 0.1994711 0.1760327
```
