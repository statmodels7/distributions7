# Third-Order Derivatives of a Reparametrized Distribution

The partition sum at order three, observed or expected. Every order the
parent has in closed form survives in closed form, so a family
reparametrized from one of the 42 that carry analytic third derivatives
keeps them.

## Usage

``` r
reparam_deriv3(
  distrib,
  y,
  theta,
  expected = FALSE,
  scale = c("parameter", "link"),
  approx = c("integrate", "bartlett", "mc", "opg"),
  nsim = 10000,
  ...
)
```

## Arguments

- distrib:

  A reparametrized distribution, from
  [`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md).

- y:

  The response, a numeric vector.

- theta:

  A named list of the new parameters, on the new parameter scale.

- expected:

  Should the expected derivatives be carried? A single logical, `FALSE`
  by default. Under expectation the term in the parent's score drops,
  the score having mean zero.

- scale:

  Either `"parameter"` (the default) or `"link"`, handled by the generic
  after this method has returned.

- approx:

  The strategy the parent uses for an expected derivative it has no
  closed form for, one of `"integrate"`, `"bartlett"`, `"mc"` or
  `"opg"`. Passed straight through; read only when `expected` is `TRUE`
  and only by a parent that approximates.

- nsim:

  The number of draws for `approx = "mc"`, 10000 by default. Passed
  through.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors keyed as
[`deriv_names(distrib@params, 3)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md),
each the length of `y` recycled against `theta`.

## Details

The sum reads the map's partials up to order three, which is where a
hand-written table begins to pay against the stencil: the numerical
route's accuracy fades with the order, and by the third it is several
digits short of the parent's own.

## See also

[`chain_derivatives()`](https://statmodels7.github.io/distributions7/reference/chain_derivatives.md)
for the identity;
[`distrib_hessian.ReparamContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.ReparamContinuousDistrib.md)
for the order below;
[`reparam_tables()`](https://statmodels7.github.io/distributions7/reference/reparam_tables.md),
whose accuracy this inherits;
[`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md).

## Examples

``` r
d <- reparametrize(
  gaussian1_distrib(),
  map = function(psi) list(mu = psi$mu, sigma = sqrt(psi$sigma2)),
  params = c("mu", "sigma2"),
  bounds = list(mu = c(-Inf, Inf), sigma2 = c(0, Inf)),
  links = list(mu = linkfunctions7::identity_link(),
               sigma2 = linkfunctions7::log_link())
)
th <- list(mu = 1, sigma2 = 4)

# Four components for a two-parameter family.
names(distrib_deriv3(d, c(0, 1, 2), th))
#> [1] "mu_mu_mu"             "mu_mu_sigma2"         "mu_sigma2_sigma2"    
#> [4] "sigma2_sigma2_sigma2"
```
