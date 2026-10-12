# Write a Distribution in Different Coordinates

Returns the same law as `distrib`, parametrized by quantities of the
caller's choosing.

The function is meant for trying a parametrization the package does not
ship. Every derivative it returns is assembled from the parent's
derivatives through the map, one order at a time, so it costs more than
a family written for that parametrization; a parametrization used in a
real analysis is better written as a distribution of its own, with its
derivatives in compiled code, which is how the families of this package
are written.

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
  name = NULL,
  params_smooth = NULL
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

  An optional function of two arguments, the named list of new
  parameters and `order`, an integer from 1 to 4, returning for each
  parent parameter the non-zero partial derivatives of the map with
  respect to the new parameters up to order `order`, keyed by the sorted
  tuple of new-parameter positions ("1", "1,2", "2,2,3,3", ...); a
  missing key is an exact zero. Each caller passes the highest order
  that it reads, so a gradient passes order one and a fourth derivative
  passes order four. The shipped second parametrizations supply
  hand-written tables (see
  [`reparam_map_derivs()`](https://statmodels7.github.io/distributions7/reference/reparam_map_derivs.md));
  when `NULL`, each needed partial comes from one finite-difference
  stencil on the map.

- interpretation:

  An optional named character vector describing each new parameter;
  defaults to the parameter names.

- name:

  An optional name for the result; defaults to the parent's with the new
  parameters appended.

- params_smooth:

  An optional named logical vector, one entry per new parameter, saying
  which of them the log-density is differentiable in. The map is
  arbitrary, so the default cannot read which new parameter reaches a
  kink of the parent: it is `TRUE` for every new parameter when every
  parameter of the parent is smooth, and `FALSE` for every one
  otherwise.

## Value

A distribution object of class
[`ReparamContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/ReparamContinuousDistrib.md)
or
[`ReparamDiscreteDistrib()`](https://statmodels7.github.io/distributions7/reference/ReparamContinuousDistrib.md).

## Details

A link changes the scale on which a parameter is modeled and leaves the
parameter itself unchanged. In a reparametrization the parameter is the
new quantity, and the estimate, the standard error and the confidence
interval describe that quantity.

**The map is written in ordinary R.** It takes a named list of the new
parameters and returns a named list of the parent's, and nothing in it
mentions derivatives:


    function(psi) list(mu = psi$mean / gamma(1 + 1 / psi$sigma),
                       sigma = psi$sigma)

**Where the map's derivatives come from.** Two routes, and `map_derivs`
chooses between them. Supply it, and the tables are hand-written closed
forms, exact at every order with no chain rule transcribed; the shipped
second parametrizations all do, and their formulas live in
[`reparam_map_derivs()`](https://statmodels7.github.io/distributions7/reference/reparam_map_derivs.md).
Leave it `NULL`, and each partial the chain needs comes from one central
stencil of
[`numericals7::fd_derivative()`](https://statmodels7.github.io/numericals7/reference/fd_derivative.html)
applied to the analytic map, a single stencil per order and never a
chain of differences. The stencil route is accurate to about
\\10^{-11}\\ at first order and loses accuracy quickly with the order
(the third and fourth tables can be wrong in their leading digits), so a
family used for anything beyond a first exploration needs `map_derivs`.

**What is exact and what is inherited.** The derivatives of the
log-density are carried by \$\$\ell^{(I)}(\psi) = \sum\_{\pi} \sum\_{i_1
\dots i\_{\|\pi\|}} \ell^{(i_1 \dots i\_{\|\pi\|})}(\theta) \prod\_{B
\in \pi} \frac{\partial^{\|B\|}\theta\_{i_B}}{\partial \psi_B},\$\$ the
sum over set partitions, so every order the parent has in closed form
survives in closed form. Expectation being linear and the map
deterministic, the expected derivatives obey the same formula, and the
term carrying the parent's score (multiplied by a derivative of the map
of the order being formed) drops because the score has mean zero: a
parent with an exact expected information gives an exact one here. The
density, the distribution function, the quantile function, the generator
and the response derivatives are the parent's at the mapped parameters.

**Cost.** The density, the distribution function, the quantile function
and the generator call the map once per evaluation, on whole vectors.
When `map_derivs` is `NULL`, the derivatives call the map at the points
of the stencil, several times for each observation (more at higher
orders), and separately for each observation when a parameter varies by
observation. A supplied `map_derivs` is called once per evaluation.

## See also

[`fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.md),
[`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md)

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
